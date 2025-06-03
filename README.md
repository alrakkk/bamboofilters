# bamboofilters for bioinformatics

A high-performance approximate membership query data structure for k-mer analysis in genomic data. Bamboo filters implement insertion, lookup, and deletion operations with constant-time cost, and support dynamic resizing through expanding and compressing the hash table.

## Features

- **High throughput** for insertions, lookups, and deletions
- **Low false positive rates** (configurable based on memory usage)
- **Memory efficient** storage of k-mers
- **Dynamic resizing** capability without rebuilding the entire structure
- **Optimized for genomic data** with support for variable k-mer sizes

## Installation

```bash
# Clone the repository
git clone https://github.com/alrakkk/bamboofilters.git
cd bamboofilters

# Create build directory
mkdir build
cd build

# Configure and build
cmake ..
make
```

## E. coli Benchmark

The Bamboo Filter includes tools to benchmark performance using the E. coli reference genome.

### Download E. coli Genome

```bash
# From the build directory
make download_ecoli
```

## Usage

The main benchmark tool is `bamboofilter_ecoli` which supports various operations:

```bash
# Basic usage with default parameters
./tools/bamboofilter_ecoli --fasta-path /path/to/ecoli.fasta --auto-test --false-positives


# Test with different k-mer sizes
./tools/bamboofilter_ecoli --fasta-path /path/to/ecoli.fasta --kmer-size 15
./tools/bamboofilter_ecoli --fasta-path /path/to/ecoli.fasta --kmer-size 25

# Test false positive rates
./tools/bamboofilter_ecoli --fasta-path /path/to/ecoli.fasta --false-positives --kmer-size 25

# Test deletion performance
./tools/bamboofilter_ecoli --fasta-path /path/to/ecoli.fasta --action delete --kmer-size 25

# Control maximum number of k-mers to process
./tools/bamboofilter_ecoli --fasta-path /path/to/ecoli.fasta --kmer-size 25 --max-kmers 50000
```

## Understanding Benchmark Results

The benchmark tool outputs detailed performance metrics to `bamboo_filter_results.csv` with the following columns:

- `kmer_size`: The size of k-mers being tested
- `num_elements`: Number of k-mers inserted
- `insert_time_ms`: Time taken to insert all k-mers (ms)
- `insert_throughput`: Insert operations per second
- `lookup_time_ms`: Time taken for lookups (ms)
- `lookup_throughput`: Lookup operations per second
- `false_positive_rate`: Percentage of false positives (lower is better)
- `peak_memory_mb`: Maximum memory usage (MB)
- `expand_memory_mb`: Memory growth during expansion
- `bits_per_element`: Memory efficiency metric

## Comparative Evaluation

Bamboo Filters are compared against these existing data structures:

| Algorithm | Description                                                                                                                                                                                                                                                                                               |
| :-------: | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
|    DBF    | D. Guo, J. Wu, H. Chen, Y. Yuan, and X. Luo, "The Dynamic Bloom Filters," IEEE Transactions on Knowledge and Data Engineering, vol. 22, no. 1, pp. 120–133, 2009. Implementation: https://github.com/CGCL-codes/DCF/tree/master/src_DBF                                                                   |
|    EBF    | Y. Wu, J. He, S. Yan, J. Wu, T. Yang, O. Ruas, G. Zhang, and B. Cui, "Elastic Bloom Filter: Deletable and ExpandableFilter Using Elastic Fingerprints," IEEE Transactions on Computers, vol. 1, no. 1, pp. 1–8, 2021. Implementation: https://github.com/Dustin-He/ElasticBloomFilter/blob/master/bloom.h |
|    SBF    | K. Xie, Y. Min, D. Zhang, J. Wen, and G. Xie, "A Scalable Bloom Filter for Membership Queries," in Proceedings of the Global Communications Conference. IEEE, 2007, pp. 543–547. Implementation: https://github.com/Dustin-He/ElasticBloomFilter/blob/master/ScalableBF.h                                 |
|    DCF    | H. Chen, L. Liao, H. Jin, and J. Wu, "The Dynamic Cuckoo Filter," in Proceedings of International Conference on Network Protocols. IEEE, 2017, pp. 1–10. Implementation: https://github.com/CGCL-codes/DCF                                                                                                |
|   SuRF    | H. Zhang, H. Lim, V. Leis, D. G. Andersen, M. Kaminsky, K. Keeton, and A. Pavlo, "SuRF: Practical Range Query Filtering with Fast Succinct Tries," in Proceedings of International Conference on Management of Data. ACM, 2018, pp. 323–336. Implementation: https://github.com/efficient/SuRF            |
|   E2CF    | S. Yu, S. Wu,H. Chen, and H. Jin, "The entry-extensible cuckoo filter," in Proceedings of International Federation for Information Processing. Springer, 2020, pp. 373–385. Implementation: https://github.com/CGCL-codes/E2CF                                                                            |
|    TCF    | J. Apple, "Stretching your data with taffy filters," arXiv preprint arXiv:2109.01947v4, 2022. Implementation: https://github.com/jbapple/libfilter                                                                                                                                                        |

## Performance Results

Our benchmarks show that Bamboo Filters outperform existing alternatives:

- Highest insertion throughput (5-10M operations/second)
- Superior positive lookup throughput
- Efficient deletion operations (~5.4M operations/second)
- Minimal performance degradation after multiple resizing operations
- Lower false positive rates for genomic data than comparable structures
- Memory efficiency comparable to state-of-the-art AMQ data structures
