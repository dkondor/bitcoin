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

See http://www.vo.elte.hu/bitcoin for research overview and downloadable processed data.


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

Development Process
-------------------

The `master` branch is regularly built and tested, but is not guaranteed to be
completely stable. [Tags](https://github.com/bitcoin/bitcoin/tags) are created
regularly to indicate new official, stable release versions of Bitcoin Core.

The contribution workflow is described in [CONTRIBUTING.md](CONTRIBUTING.md).

The developer [mailing list](https://lists.linuxfoundation.org/mailman/listinfo/bitcoin-dev)
should be used to discuss complicated or controversial changes before working
on a patch set.

Developer IRC can be found on Freenode at #bitcoin-core-dev.

Testing
-------

Testing and code review is the bottleneck for development; we get more pull
requests than we can review and test on short notice. Please be patient and help out by testing
other people's pull requests, and remember this is a security-critical project where any mistake might cost people
lots of money.

### Automated Testing

Developers are strongly encouraged to write [unit tests](src/test/README.md) for new code, and to
submit new unit tests for old code. Unit tests can be compiled and run
(assuming they weren't disabled in configure) with: `make check`. Further details on running
and extending unit tests can be found in [/src/test/README.md](/src/test/README.md).

There are also [regression and integration tests](/test), written
in Python, that are run automatically on the build server.
These tests can be run (if the [test dependencies](/test) are installed) with: `test/functional/test_runner.py`

The Travis CI system makes sure that every pull request is built for Windows, Linux, and OS X, and that unit/sanity tests are run automatically.

### Manual Quality Assurance (QA) Testing

Changes should be tested by somebody other than the developer who wrote the
code. This is especially important for large or high-risk changes. It is useful
to add a test plan to the pull request description if testing the changes is
not straightforward.

Translations
------------

Changes to translations as well as new translations can be submitted to
[Bitcoin Core's Transifex page](https://www.transifex.com/projects/p/bitcoin/).

Translations are periodically pulled from Transifex and merged into the git repository. See the
[translation process](doc/translation_process.md) for details on how this works.

**Important**: We do not accept translation changes as GitHub pull requests because the next
pull from Transifex would automatically overwrite them again.

Translators should also subscribe to the [mailing list](https://groups.google.com/forum/#!forum/bitcoin-translators).
