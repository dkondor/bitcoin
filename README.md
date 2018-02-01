bitcoind-dump-tsv
-----------------

Bitcoin Core (Satoshi) client with minimal modifications to output the whole
blockchain as TSV files.

This is achieved in a special "dump" mode, activated with the -DUMP command
line argument, where the blockchain is read from the database and written to
files specified as further argument, after which the program will exit.

All other functionality is unchanged, should be compatible with an unmodified
Bitcoin Core client in any way (including using the same data directory).
Modifications are based on version 0.15.1.

For original build instructions, see https://github.com/bitcoin/bitcoin

or

https://bitcoincore.org


Usage (dump mode)
-----------------

src/bitcoind -DUMP -DUMP_outdir=/path/to/output

The following files will be created or overwritten in the output dir:

	- txout.dat -- output file for transaction outputs; format: txID, output_seq, addrID, sum
	- txin.dat -- output file for transaction inputs; format: txID, input_seq, prev_txID, prev_output_seq, prev_i, addrID, sum
	- tx.dat -- output file for transaction overview (mapping to block and number of inputs / outputs); format: txID, blockID, n_inputs, n_outputs
	- txh.dat -- output file for transaction hashes (mapping from txIDs used in all other outputs); format: txID, hash
	- bh.dat -- output file for block hashes (mapping from blockIDs used in all other outputs); format: blockID, hash, block_timestamp, n_txs
	- missing.dat -- output file for missing transaction inputs (should be empty, anything in it is an error); format: txID, input_seq
	- multiple.dat -- output file for transaction outputs with multiple addresses (multisign);
		the txout and txin file will only include the first address; format: txID, output_seq, addrID
	- nonstandard.dat -- output file for nonstandard transaction outputs; format: txID, output_seq
	- addresses.dat -- output file for address ID mapping to address strings (IDs are used in all other outputs); format: addrID, address


All output files use numeric IDs to refer to transactions, blocks and addresses (these are all counters starting from 0). These are mapped to
hashes in the three files txh.dat, bh.dat and addresses.dat respectively. A special value of -1 for txID means a bug in the processing (should not
happen). A special value of -1 for an addrID means that the address could not be decoded. This is not necessarily and error, there are certain
nonstandard transactions where this can happen. The flow of bitcoins can still be followed in these cases as all transaction inputs are linked
to the corresponding previous transactions outputs in the txin.dat.

All sums are in Satoshis (1e-9 BTC).

Transaction inputs and outputs include a sequence number (input_seq and output_seq respectively), which identifies the input / output. These are counters
starting from 0 for each transaction. I'm not sure if these will correspond to the same used by other Bitcoin clients, but can be used to map inputs to
previous outputs. The txin.dat file includes this information, the prev_txID and prev_output_seq columns include this information.

Mining rewards (coinbase transaction) can be identified by having zero inputs. For all other transactions, the sum of inputs should be greater than the
sum of outputs, but this is not checked explicitely during processing.



What is Bitcoin?
----------------

Bitcoin is an experimental digital currency that enables instant payments to
anyone, anywhere in the world. Bitcoin uses peer-to-peer technology to operate
with no central authority: managing transactions and issuing money are carried
out collectively by the network. Bitcoin Core is the name of open source
software which enables the use of this currency.

For more information, as well as an immediately useable, binary version of
the Bitcoin Core software, see https://bitcoin.org/en/download, or read the
[original whitepaper](https://bitcoincore.org/bitcoin.pdf).

License
-------

Bitcoin Core is released under the terms of the MIT license. See [COPYING](COPYING) for more
information or see https://opensource.org/licenses/MIT.

