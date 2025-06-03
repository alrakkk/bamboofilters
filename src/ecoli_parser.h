/* ecoli_parser.h */
#pragma once
#include <cstdint>
#include <string>
#include <vector>

// Read FASTA / .fna -> returns sequence length, sequence in seq_out.
size_t read_fasta_sequence(const std::string& path, std::string& seq_out);

// Hash every k-mer (k ≥ 1) in the sequence into a 64-bit vector.
std::vector<uint64_t> hash_kmers(const std::string& sequence, uint32_t k);
