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

#include <string>
#include <string.h>
#include <unordered_map>


struct txoutput { // store transaction outputs (for the sake of remembering them when they are spent later)
	txoutput():addr(-1),txid(-1),value(0) {}
	txoutput(int64_t addr_, int64_t val, int64_t txid_):addr(addr_),txid(txid_),value(val) {  }
	int64_t txid; // transaction ID
	int64_t addr; // address ID
	int64_t value; // output value (in Satoshis, must be positive, but not checked)
};

struct txoutid { //ezek szerepelnek a bemenetnél
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

char* path_combine(const char* str1, const char* str2) {
	int n1 = strlen(str1);
	int n2 = strlen(str2);
	char* ret = (char*)malloc(sizeof(char)*(n1+n2+2));
	if(!ret) return 0;
	sprintf(ret,"%s/%s",str1,str2);
	return ret;
}

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
	char* outdir = 0; // output directory, if given, just use default filenames
	int64_t bmax = chainActive.Height() + 1; // last block to write out
	int64_t T1 = 10000; // output progress after this many blocks each
	int64_t Tn = T1;

	int64_t txs = 0; // counter for transactions processed (used as txIDs as well)

	txoutmap outmap1;
	std::unordered_map<std::string,int64_t> addr_map;
	
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
	
	outdir = gArgs.GetArg("-DUMP_outdir");
	if(outdir) {
		if(!ftxout) {
			ftxout = path_combine(outdir,"txout.dat");
			if(!ftxout) { fprintf(stderr,"dumpblocks(): Error allocating memory!\n"); goto dumpblocks_end; }
		}
		if(!ftxin) {
			ftxin = path_combine(outdir,"txin.dat");
			if(!ftxin) { fprintf(stderr,"dumpblocks(): Error allocating memory!\n"); goto dumpblocks_end; }
		}
		if(!ftx) {
			ftx = path_combine(outdir,"tx.dat");
			if(!ftx) { fprintf(stderr,"dumpblocks(): Error allocating memory!\n"); goto dumpblocks_end; }
		}
		if(!ftxh) {
			ftxh = path_combine(outdir,"txh.dat");
			if(!ftxh) { fprintf(stderr,"dumpblocks(): Error allocating memory!\n"); goto dumpblocks_end; }
		}
		if(!fbh) {
			fbh = path_combine(outdir,"bh.dat");
			if(!fbh) { fprintf(stderr,"dumpblocks(): Error allocating memory!\n"); goto dumpblocks_end; }
		}
		if(!fmissing) {
			fmissing = path_combine(outdir,"missing.dat");
			if(!fmissing) { fprintf(stderr,"dumpblocks(): Error allocating memory!\n"); goto dumpblocks_end; }
		}
		if(!ftmultiple) {
			ftmultiple = path_combine(outdir,"multiple.dat");
			if(!ftmultiple) { fprintf(stderr,"dumpblocks(): Error allocating memory!\n"); goto dumpblocks_end; }
		}
		if(!fnonstandard) {
			fnonstandard = path_combine(outdir,"nonstandard.dat");
			if(!fnonstandard) { fprintf(stderr,"dumpblocks(): Error allocating memory!\n"); goto dumpblocks_end; }
		}
		if(!faddresses) {
			faddresses = path_combine(outdir,"addresses.dat");
			if(!faddresses) { fprintf(stderr,"dumpblocks(): Error allocating memory!\n"); goto dumpblocks_end; }
		}
	}
	
	bmax = gArgs.GetArg("-DUMP_bmax",bmax);
	
	
	if( ftxout == 0 || ftxin == 0 || ftx == 0 || ftxh == 0 || fbh == 0
	   || fmissing == 0 || fnonstandard == 0 || faddresses == 0 || ftmultiple == 0 ) {
		fprintf(stderr,"dumpblocks(): error: not all filenames given!\n");
		goto dumpblocks_end;
	}

	stxout = fopen(ftxout,"w");
	stxin = fopen(ftxin,"w");
	stx = fopen(ftx,"w");
	stxh = fopen(ftxh,"w");
	sbh = fopen(fbh,"w");
	smissing = fopen(fmissing,"w");
	snonstandard = fopen(fnonstandard,"w");
	saddresses = fopen(faddresses,"w");
	stmultiple = fopen(ftmultiple,"w");

	if( ! ( stxout && stxin && stx && stxh && sbh && smissing && snonstandard && saddresses && stmultiple )) {
		fprintf(stderr,"dumpblocks(): error opening output files!");
		goto dumpblocks_end;
	}
	
	for(;b<bmax;b++) {
		if(b > chainActive.Height()) break; // already at the end of blockchain
		CBlock block;
		CBlockIndex* pblockindex = chainActive[b];
		if(!pblockindex) {
			fprintf(stderr,"Error finding block %ld!\n",b);
			goto dumpblocks_end; 
		}
		
		std::string bhash = pblockindex->GetBlockHash().GetHex();
		int64_t btime = pblockindex->GetBlockTime();
		if(!ReadBlockFromDisk(block, pblockindex, par)) {
			fprintf(stderr,"Error reading block %ld (hash: %s)!\n",b,bhash.c_str());
			goto dumpblocks_end;
		}
		
		unsigned int btxs = 0;

		// process all transactions in this block
		for(unsigned int i1=0;i1<block.vtx.size();i1++) {
			const CTransaction& tx = *(block.vtx[i1]);

			std::string txh = tx.GetHash().GetHex();
		    //~ unsigned int txtime = btime; // note: no transaction timestamps stored in the blockchain, so we are using the block timestamps

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
						fprintf(smissing,"%ld\t%s\t%u\n",txs,outhash.ToString().c_str(),i);
					}
					else {
						const txoutput& out2 = out1->second;
						// output: txID, i, prev_txID, prev_i, addr_id, value (Satoshis)
						fprintf(stxin,"%ld\t%u\t%ld\t%u\t%ld\t%ld\n",txs,i,out2.txid,ivout,out2.addr,out2.value);
					}
					
					// delete previous output from the map of tx outputs
					outmap1.erase(out1);
					txin++;
				} //for(i) -- process all tx inputs
			}

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
							const std::string& s = CBitcoinAddress(it).ToString();
							int64_t addrid;
							auto it2 = addr_map.find(s);
							if(it2 == addr_map.end()) {
								addr_map.insert(std::make_pair(s,naddr));
								addrid = naddr;
								fprintf(saddresses,"%ld\t%s\n",naddr,s.c_str());
								naddr++;
							}
							else addrid = it2->second;
							fprintf(stmultiple,"%ld\t%u\t%ld\n",txs,i,addrid);
						}
					}
				}
				
				int64_t addrid1;
				if(type == TX_NONSTANDARD || addresses.size() == 0) {
					fprintf(snonstandard,"%ld\t%u\n",txs,i);
					// use special addr_id for these, but include in normal output
					// so spending can be linked later from the txin file
					addrid1 = -1;
				}
				else {
					// write out only the first output address as the addressee in the main output file
					const std::string& addr = CBitcoinAddress(addresses[0]).ToString();
					auto it1 = addr_map.find(addr);
					if(it1 == addr_map.end()) {
						addr_map.insert(std::make_pair(addr,naddr));
						addrid1 = naddr;
						fprintf(saddresses,"%ld\t%s\n",naddr,addr.c_str());
						naddr++;
					}
					else addrid1 = it1->second;
				}
				
				// output: txID, i, addr_id, sum (in Satoshis)
				fprintf(stxout,"%ld\t%u\t%ld\t%ld\n",txs,i,addrid1,value);
				txout++;

				// save the current transaction output to be found later when it is spent
				txoutput out1(addrid1,value,txs);
				txoutid id1(tx.GetHash(),i);

				outmap1[id1] = out1;
			}//for(i) -- process all outputs
			//~ if(err) {
				//~ break;
			//~ }

			// output for each transaction once: txID, blockID, number of inputs and outputs
			fprintf(stx,"%ld\t%ld\t%d\t%d\n",txs,b,txin,txout);

			// separate output (txhash file): txID, tx hash
			fprintf(stxh,"%ld\t%s\n",txs,txh.c_str());
			txs++;
			btxs++;
		} //for(tx) -- process all transactions in the current block

		if(pblockindex->nTx != btxs) {
			fprintf(stderr,"missing transactions for block %ld (hash: %s)!\n",b,bhash.c_str());
			break;
		}

		// print out block properties
		fprintf(sbh,"%ld\t%s\t%ld\t%u\n",b,bhash.c_str(),btime,btxs);
		//~ blocks++;
		
		
		if(b == Tn) {
			fprintf(stderr,"%ld blocks processed\n",b);
			Tn += T1;
			fflush(stderr);
		}
		if(ShutdownRequested()) break;
	} //for(b) -- process all blocks
	
	
dumpblocks_end:
	if(sbh) fclose(sbh); if(fbh) free(fbh);
	if(stx) fclose(stx); if(ftx) free(ftx);
	if(stxh) fclose(stxh); if(ftxh) free(ftxh);
	if(stxin) fclose(stxin); if(ftxin) free(ftxin);
	if(stxout) fclose(stxout); if(ftxout) free(ftxout);
	if(smissing) fclose(smissing); if(fmissing) free(fmissing);
	if(snonstandard) fclose(snonstandard); if(fnonstandard) free(fnonstandard);
	if(stmultiple) fclose(stmultiple); if(ftmultiple) free(ftmultiple);
	if(outdir) free(outdir);

	fprintf(stderr,"dumpblocks(): Processed %ld blocks and %ld transactions\n",b,txs);

	return 0;
}


