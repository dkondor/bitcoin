/*
 * dump2.cpp -- extract the whole Bitcoin blockchain from the database of
 * 	the Satoshi client, dump everything as TSV files
 *  note: all output file names have to be specified as command line args,
 * 	look for the arguments like -DUMP_* below
 *
 * Copyright 2013-2018 Kondor Dániel <kondor.dani@gmail.com>
 *
 * Distributed under the MIT software license, see the accompanying
 * file COPYING or http://www.opensource.org/licenses/mit-license.php.
 *
 */

//~ #include "bitcoinrpc.h"
#include "net.h"
#include "init.h"
#include "util.h"
#include "chain.h"
#include "base58.h"
#include "validation.h"
#include "chainparams.h"
#include "script/standard.h"
#include "utilstrencodings.h"

#include "read_table.h"

#include <string>
#include <string.h>
#include <unordered_map>


struct txoutput { // store transaction outputs (for the sake of remembering them when they are spent later)
	txoutput():txid(-1),addr(-1),value(0) {}
	txoutput(int64_t addr_, int64_t val, int64_t txid_):txid(txid_),addr(addr_),value(val) {  }
	int64_t txid; // transaction ID
	int64_t addr; // address ID
	int64_t value; // output value (in Satoshis, must be positive, but not checked)
};

struct txoutid { //ezek szerepelnek a bemenetnél
	txoutid() {}
	txoutid(uint256 id, uint32_t n) { txid = id; nout = n; }
	uint256 txid;
	uint32_t nout;
};


// hash function for transaction outputs
struct txouthash {
	uint64_t operator()( const txoutid &p) const {
		uint64_t r = p.txid.GetCheapHash();
		uint64_t n = p.nout;
		r ^= n + 0xc6a4a7935bd1e995UL + (r<<6) + (r>>2);
		return r;
	}
};

struct txouteq {
	bool operator()( const txoutid& p1, const txoutid& p2 ) const {
		if(p1.txid != p2.txid) return false;
		if(p1.nout != p2.nout) return false;
		return true;
	}
};

typedef std::unordered_map<txoutid,txoutput,txouthash,txouteq> txoutmap;

static char* path_combine(const char* str1, const char* str2) {
	int n1 = strlen(str1);
	int n2 = strlen(str2);
	char* ret = (char*)malloc(sizeof(char)*(n1+n2+2));
	if(!ret) return 0;
	sprintf(ret,"%s/%s",str1,str2);
	return ret;
}

static FILE* open_out_zip(const char* fn) {
	const char gzip[] = "/bin/gzip -c > ";
	char* cmd = (char*)malloc( sizeof(char) * (strlen(fn) + strlen(gzip) + 8));
	if(!cmd) return 0;
	sprintf(cmd,"%s%s.gz",gzip,fn);
	FILE* ret = popen(cmd,"w");
	free(cmd);
	return ret;
}

static FILE* open_in_zip(const char* fn) {
	const char gzip[] = "/bin/gzip -cd ";
	char* cmd = (char*)malloc( sizeof(char) * (strlen(fn) + strlen(gzip) + 8));
	if(!cmd) return 0;
	sprintf(cmd,"%s%s",gzip,fn);
	FILE* ret = popen(cmd,"r");
	free(cmd);
	return ret;
}

/* alternative way to store addresses in a hashtable: store the original boost::variant values
 * hash function -- use the lowest 64 bits for all, do not care about type
 * (hash collisions can be resolved by the equal function) */
class addr_v_hash : public boost::static_visitor<size_t> {
public:
	size_t operator()(const CKeyID& id) const { return id.GetCheapHash(); }
    size_t operator()(const CScriptID& id) const { return id.GetCheapHash(); }
    size_t operator()(const WitnessV0KeyHash& id) const { return id.GetCheapHash(); }
    size_t operator()(const WitnessV0ScriptHash& id) const { return id.GetCheapHash(); }
    size_t operator()(const WitnessUnknown& id) const { return ReadLE64(id.program); }
    size_t operator()(const CNoDestination& no) const { return 0UL; }
    size_t operator()(const CTxDestination& id) const { return boost::apply_visitor(*this,id); }
};
/* compare for equality -- only objects of the same type will compare equal, taken from the boost::variant tutorial */
class addr_v_eq : public boost::static_visitor<bool> {
public:
	template <typename T, typename U>
    bool operator()( const T &, const U & ) const { return false; }
    template <typename T> /* all potential types have equality defined */
    bool operator()( const T & lhs, const T & rhs ) const { return lhs == rhs; }
};
typedef std::unordered_map<CTxDestination,int64_t,addr_v_hash,addr_v_eq> addr_v_map;

int dumpblocks()
{
	char* ftxout = 0; FILE* stxout = 0; // output file for transaction outputs
	char* ftxin = 0; FILE* stxin = 0; // output file for transaction inputs
	char* ftx = 0; FILE* stx = 0; // output file for transaction overview (mapping to block and number of inputs / outputs)
	char* ftxh = 0; FILE* stxh = 0; // output file for transaction hashes (mapping from txIDs used in all other outputs)
	char* fbh = 0; FILE* sbh = 0; // output file for block hashes (mapping from blockIDs used in all other outputs)
	char* fmissing = 0; FILE* smissing = 0; // output file for missing transaction inputs (should be empty, anything in it is an error)
	char* ftmultiple = 0; FILE* stmultiple = 0; // output file for transaction outputs with multiple addresses (multisign) -- the txout and txin file will only include the first address
	char* fnonstandard = 0; FILE* snonstandard = 0; // output file for nonstandard transaction outputs
	char* faddresses = 0; FILE* saddresses = 0; // output file for address ID mapping to address strings (IDs are used in all other outputs)
	char* funspent = 0; FILE* sunspent = 0; // output file with unspent transaction outputs (can be read later)
	char* outdir = 0; // output directory, if given, just use default filenames
	int64_t bmax = chainActive.Height() + 1; // last block to write out
	int64_t bmin = 0; // first block to write out
	int64_t txmin = 0; // first transaction ID to use
	int64_t T1 = 10000; // output progress after this many blocks each
	int64_t Tn = T1;
	int err = 0;

	bool out_zip = gArgs.GetBoolArg("-DUMP_zip", false);

	int64_t txs = 0; // counter for transactions processed (used as txIDs as well)

	txoutmap outmap1;
	addr_v_map addr_map;
	
	int64_t b = 0; // current block index
	int64_t naddr = 0; // number of addresses encountered in total
	
	const auto& cp = Params();
	const auto& par = cp.GetConsensus();
	
	
	//get output file names or output dir name
	ftxout = gArgs.GetArg("-DUMP_txout"); // output file for transaction outputs
	ftxin = gArgs.GetArg("-DUMP_txin"); // output file for transaction inputs
	ftx = gArgs.GetArg("-DUMP_tx"); // output file for transaction overview (mapping to block and number of inputs / outputs)
	ftxh = gArgs.GetArg("-DUMP_txh"); // output file for transaction hashes (mapping from txIDs used in all other outputs)
	fbh = gArgs.GetArg("-DUMP_bh"); // output file for block hashes (mapping from blockIDs used in all other outputs)
	fmissing = gArgs.GetArg("-DUMP_missing"); // output file for missing transaction inputs (should be empty, anything in it is an error)
	ftmultiple = gArgs.GetArg("-DUMP_multiple"); // output file for transaction outputs with multiple addresses (multisign) -- the txout and txin file will only include the first address
	fnonstandard = gArgs.GetArg("-DUMP_nonstandard"); // output file for nonstandard transaction outputs
	faddresses = gArgs.GetArg("-DUMP_addresses"); // output file for address ID mapping to address strings (IDs are used in all other outputs)
	funspent = gArgs.GetArg("-DUMP_unspent");
	
	outdir = gArgs.GetArg("-DUMP_outdir");
	if(outdir) {
		if(!ftxout) ftxout = path_combine(outdir,"txout.dat");
		if(!ftxin) ftxin = path_combine(outdir,"txin.dat");
		if(!ftx) ftx = path_combine(outdir,"tx.dat");
		if(!ftxh) ftxh = path_combine(outdir,"txh.dat");
		if(!fbh) fbh = path_combine(outdir,"bh.dat");
		if(!fmissing) fmissing = path_combine(outdir,"missing.dat");
		if(!ftmultiple) ftmultiple = path_combine(outdir,"multiple.dat");
		if(!fnonstandard) fnonstandard = path_combine(outdir,"nonstandard.dat");
		if(!faddresses) faddresses = path_combine(outdir,"addresses.dat");
		if(!funspent) funspent = path_combine(outdir,"txout_unspent.dat");
		
		if(! (ftxout && ftxin && ftx && ftxh && fbh && fmissing && ftmultiple && fnonstandard && faddresses && funspent) ) {
			fprintf(stderr,"dumpblocks(): Error allocating memory!\n");
			goto dumpblocks_end;
		}
	}
	
	bmax = gArgs.GetArg("-DUMP_bmax",bmax);
	bmin = gArgs.GetArg("-DUMP_bmin",0);
	txmin = gArgs.GetArg("-Dump_txmin",0);
	if(txmin) b = bmin;
	
	if( ftxout == 0 || ftxin == 0 || ftx == 0 || ftxh == 0 || fbh == 0
	   || fmissing == 0 || fnonstandard == 0 || faddresses == 0 || ftmultiple == 0 ) {
		fprintf(stderr,"dumpblocks(): error: not all filenames given!\n");
		goto dumpblocks_end;
	}

	if(out_zip) {
		stxout = open_out_zip(ftxout);
		stxin = open_out_zip(ftxin);
		stx = open_out_zip(ftx);
		stxh = open_out_zip(ftxh);
		sbh = open_out_zip(fbh);
		smissing = open_out_zip(fmissing);
		snonstandard = open_out_zip(fnonstandard);
		saddresses = open_out_zip(faddresses);
		stmultiple = open_out_zip(ftmultiple);
		if(funspent) sunspent = open_out_zip(funspent);
	}
	else {
		stxout = fopen(ftxout,"w");
		stxin = fopen(ftxin,"w");
		stx = fopen(ftx,"w");
		stxh = fopen(ftxh,"w");
		sbh = fopen(fbh,"w");
		smissing = fopen(fmissing,"w");
		snonstandard = fopen(fnonstandard,"w");
		saddresses = fopen(faddresses,"w");
		stmultiple = fopen(ftmultiple,"w");
		if(funspent) sunspent = fopen(funspent,"w");
	}

	if( ! ( stxout && stxin && stx && stxh && sbh && smissing && snonstandard && saddresses && stmultiple )) {
		fprintf(stderr,"dumpblocks(): error opening output files!");
		goto dumpblocks_end;
	}
	if(funspent && !sunspent) {
		fprintf(stderr,"dumpblocks(): error opening output files!");
		goto dumpblocks_end;
	}
	
	
	/* optionally read addresses and unspent transaction outputs from previous run */
	{
		bool load_zip = gArgs.GetBoolArg("-DUMP_load_zip",false);
		char* addr_load_fn = gArgs.GetArg("-DUMP_load_addr");
		char* outmap_load_fn = gArgs.GetArg("-DUMP_load_unspent");
		
		if(addr_load_fn) {
			FILE* faddrload = 0;
			if(load_zip) faddrload = open_in_zip(addr_load_fn);
			else faddrload = fopen(addr_load_fn,"r");
			if(!faddrload) {
				fprintf(stderr,"dumpblocks(): Error opening input file %s!\n",addr_load_fn);
				goto dumpblocks_end;
			}
			read_table2 rt(faddrload);
			rt.set_fn(addr_load_fn);
			while(rt.read_line()) {
				int64_t addrid;
				const char* addrstr;
				size_t addrlen;
				if(! (rt.read_int64(addrid) && rt.read_string(&addrstr,&addrlen) ) ) break;
				if(addrid != naddr) {
					fprintf(stderr,"dumpblocks(): Address IDs in input file %s not sorted on line %lu!\n",addr_load_fn,rt.get_line());
					err = 1; break;
				}
				
				CTxDestination addr2 = DecodeDestination(std::string(addrstr,addrlen));
				if(boost::get<CNoDestination>(&addr2)) {
					fprintf(stderr,"dumpblocks(): invalid address in input file %s line %lu!\n",addr_load_fn,rt.get_line());
					err = 1;
					break;
				}
				addr_map.insert(std::make_pair(addr2,naddr));
				naddr++;
			}
			if(!err) if(rt.get_last_error() != T_EOF) {
				fprintf(stderr,"dumpblocks(): ");
				rt.write_error(stderr);
				err = 1;
			}
			if(load_zip) pclose(faddrload);
			else fclose(faddrload);
			free(addr_load_fn);
			if(err) goto dumpblocks_end;
		}
		
		if(outmap_load_fn) {
			FILE* foutmapload = 0;
			if(load_zip) foutmapload = open_in_zip(outmap_load_fn);
			else foutmapload = fopen(outmap_load_fn,"r");
			if(!foutmapload) {
				fprintf(stderr,"dumpblocks(): Error opening input file %s!\n",outmap_load_fn);
				goto dumpblocks_end;
			}
			read_table2 rt(foutmapload);
			rt.set_fn(outmap_load_fn);
			while(rt.read_line()) {
				txoutid txoid;
				txoutput txout;
				const char* txhashstr;
				size_t txhashlen;
				if(!rt.read_string(&txhashstr,&txhashlen)) break;
				if(!rt.read(txoid.nout,txout.txid,txout.addr,txout.value)) break;
				if(txhashlen != 64) {
					fprintf(stderr,"dumpblocks(): Invalid transaction hash value in input file %s, linr %lu: %.*s!\n",outmap_load_fn,rt.get_line(),(int)txhashlen,txhashstr);
					err = 1; break;
				}
				for(int i=0;i<64;i++) if(HexDigit(txhashstr[i]) == -1) {
					fprintf(stderr,"dumpblocks(): Invalid transaction hash value in input file %s, linr %lu: %.*s!\n",outmap_load_fn,rt.get_line(),(int)txhashlen,txhashstr);
					err = 1; break;
				}
				if(err) break;
				txoid.txid.SetHex(txhashstr);
				if(outmap1.insert(std::make_pair(txoid,txout)).second == false) {
					fprintf(stderr,"dumpblocks(): transaction output appears multiple times in input file %s (line %lu)!\n",outmap_load_fn,rt.get_line());
					err = 1; break;
				}
			}
			if(rt.get_last_error() != T_EOF) {
				fprintf(stderr,"dumpblocks(): ");
				rt.write_error(stderr);
				err = 1;
			}
			if(load_zip) pclose(foutmapload);
			else fclose(foutmapload);
			free(outmap_load_fn);
			if(err) goto dumpblocks_end;
		}
	}
	
	for(;b<bmax;b++) {
		if(b > chainActive.Height()) break; // already at the end of blockchain
		CBlock block;
		CBlockIndex* pblockindex = chainActive[b];
		if(!pblockindex) {
			fprintf(stderr,"Error finding block %ld!\n",b);
			goto dumpblocks_end; 
		}
		
		if(!ReadBlockFromDisk(block, pblockindex, par)) {
			fprintf(stderr,"Error reading block %ld (hash: %s)!\n",b,pblockindex->GetBlockHash().GetHex().c_str());
			goto dumpblocks_end;
		}
		if(b < bmin) { // just count the number of transactions in old blocks
			txmin += block.vtx.size();
			txs += block.vtx.size();
			continue;
		}
		
		unsigned int btxs = 0;
		std::string bhash = pblockindex->GetBlockHash().GetHex();
		int64_t btime = pblockindex->GetBlockTime();
		
		// process all transactions in this block
		for(const auto& ptx : block.vtx) {
			const CTransaction& tx = *ptx;
			std::string txh = tx.GetHash().GetHex();

			int txin = 0; // total transaction inputs; note: we do not count newly minted Bitcoins (coinbase transactions)
			int txout = 0; // total transaction outputs

			// process transaction inputs
			if(!tx.IsCoinBase()) {
				for(uint32_t i=0;i<tx.vin.size();i++) {
					const CTxIn& ctxin = tx.vin[i];

					// find previous transaction which the current one references
					uint256 outhash = ctxin.prevout.hash;
					uint32_t ivout = ctxin.prevout.n;

					txoutid outid1(outhash,ivout);
					txoutmap::iterator out1 = outmap1.find(outid1);
					
					if(out1 == outmap1.end()) {
						// missing previous output, write to separate file (this is an error)
						err = fprintf(smissing,"%ld\t%s\t%u\n",txs,outhash.ToString().c_str(),i); if(err < 0) break;
					}
					
					else {
						const txoutput out2 = out1->second;
						// output: txID, i, prev_txID, prev_i, addr_id, value (Satoshis)
						err = fprintf(stxin,"%ld\t%u\t%ld\t%u\t%ld\t%ld\n",txs,i,out2.txid,ivout,out2.addr,out2.value); if(err < 0) break;
						
						// delete previous output from the map of tx outputs
						outmap1.erase(out1);
					}
					txin++;
				} //for(i) -- process all tx inputs
			}
			
			if(err < 0) break;

			//kimenetek feldolgozása
			for(uint32_t i=0;i<tx.vout.size();i++) {
				const CTxOut& ctxout = tx.vout[i]; 
				int64_t value = ctxout.nValue;
				std::vector<CTxDestination> addresses;
				txnouttype type;
				int nRequired;
				if (!ExtractDestinations(ctxout.scriptPubKey, type, addresses, nRequired)) {
					type = TX_NONSTANDARD;
				}
				
				if(addresses.size() > 1) {
					// a transaction output can include multiple addresses (multisign)
					// write these to a separate file
					if(stmultiple) {
						for(const auto& it : addresses) {
							int64_t addrid;
							auto it2 = addr_map.find(it);
							if(it2 == addr_map.end()) {
								addr_map.insert(std::make_pair(it,naddr));
								addrid = naddr;
								std::string addr_str = EncodeDestination(it);
								err = fprintf(saddresses,"%ld\t%s\n",naddr,addr_str.c_str()); if(err < 0) break;
								naddr++;
							}
							else addrid = it2->second;
							err = fprintf(stmultiple,"%ld\t%u\t%ld\n",txs,i,addrid); if(err < 0) break;
						}
					}
				}
				if(err < 0) break;
				
				int64_t addrid1;
				if(type == TX_NONSTANDARD || addresses.size() == 0) {
					err = fprintf(snonstandard,"%ld\t%u\n",txs,i); if(err < 0) break;
					// use special addr_id for these, but include in normal output
					// so spending can be linked later from the txin file
					addrid1 = -1;
				}
				else {
					// write out only the first output address as the addressee in the main output file
					auto it1 = addr_map.find(addresses[0]);
					if(it1 == addr_map.end()) {
						addr_map.insert(std::make_pair(addresses[0],naddr));
						addrid1 = naddr;
						std::string addr_str = EncodeDestination(addresses[0]);
						err = fprintf(saddresses,"%ld\t%s\n",naddr,addr_str.c_str()); if(err < 0) break;
						naddr++;
					}
					else addrid1 = it1->second;
				}
				
				// output: txID, i, addr_id, sum (in Satoshis)
				err = fprintf(stxout,"%ld\t%u\t%ld\t%ld\n",txs,i,addrid1,value); if(err < 0) break;
				txout++;

				// save the current transaction output to be found later when it is spent
				txoutput out1(addrid1,value,txs);
				txoutid id1(tx.GetHash(),i);

				outmap1[id1] = out1;
			}//for(i) -- process all outputs
			if(err < 0) break;

			// output for each transaction once: txID, blockID, number of inputs and outputs
			err = fprintf(stx,"%ld\t%ld\t%d\t%d\n",txs,b,txin,txout); if(err < 0) break;

			// separate output (txhash file): txID, tx hash
			err = fprintf(stxh,"%ld\t%s\n",txs,txh.c_str()); if(err < 0) break;
			txs++;
			btxs++;
		} //for(tx) -- process all transactions in the current block

		if(pblockindex->nTx != btxs) {
			fprintf(stderr,"missing transactions for block %ld (hash: %s)!\n",b,bhash.c_str());
			break;
		}

		// print out block properties
		err = fprintf(sbh,"%ld\t%s\t%ld\t%u\n",b,bhash.c_str(),btime,btxs);
		//~ blocks++;
		if(err < 0) { fprintf(stderr,"dumpblocks(): error writing output!\n"); break; }
		
		if(b == Tn) {
			fprintf(stderr,"%ld blocks processed\n",b);
			Tn += T1;
			fflush(stderr);
		}
		if(ShutdownRequested()) break;
	} //for(b) -- process all blocks
	
	if(sunspent && err >= 0) {
		if(ShutdownRequested()) fprintf(stderr,"dumpblocks(): exit request, saving unspent transaction inputs");
		for(const auto& p : outmap1) {
			err = fprintf(sunspent,"%s\t%u\t%ld\t%ld\t%ld\n",p.first.txid.GetHex().c_str(),p.first.nout,
				p.second.txid,p.second.addr,p.second.value);
			if(err < 0) { fprintf(stderr,"dumpblocks(): error writing output!\n"); break; }
		}
	}
	
	
dumpblocks_end:
	if(out_zip) {
		if(sbh) pclose(sbh);
		if(stx) pclose(stx);
		if(stxh) pclose(stxh);
		if(stxin) pclose(stxin);
		if(stxout) pclose(stxout);
		if(smissing) pclose(smissing);
		if(snonstandard) pclose(snonstandard);
		if(stmultiple) pclose(stmultiple);
		if(sunspent) pclose(sunspent);
	}
	else {
		if(sbh) fclose(sbh);
		if(stx) fclose(stx);
		if(stxh) fclose(stxh);
		if(stxin) fclose(stxin);
		if(stxout) fclose(stxout);
		if(smissing) fclose(smissing); 
		if(snonstandard) fclose(snonstandard);
		if(stmultiple) fclose(stmultiple);
		if(sunspent) fclose(sunspent);
	}
	
	
	if(outdir) free(outdir);
	if(fbh) free(fbh);
	if(ftx) free(ftx);
	if(ftxh) free(ftxh);
	if(ftxin) free(ftxin);
	if(ftxout) free(ftxout);
	if(fmissing) free(fmissing);
	if(fnonstandard) free(fnonstandard);
	if(ftmultiple) free(ftmultiple);
	if(funspent) free(funspent);
	
	fprintf(stderr,"dumpblocks(): min blockID: %ld, max blockID: %ld, min txID: %ld, max txID %ld\n",bmin,b,txmin,txs);

	return 0;
}


