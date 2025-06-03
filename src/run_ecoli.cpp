#include "bamboofilter/bamboofilter.hpp"  
#include "src/ecoli_parser.h"
#include <iostream>
#include <string>
#include <vector>
#include <chrono>
#include <cstdlib>
#include <cmath>
#include <fstream>  
#include <sys/resource.h>  
#include <iomanip>  // For formatted output
#include <random> 


// Generate mutated k-mers that aren't in the original dataset
std::vector<std::string> generate_false_positive_tests(const std::string& sequence, int kmerSize, size_t numTests) {
    std::vector<std::string> test_kmers;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> pos_dist(0, sequence.length() - kmerSize);
    std::uniform_int_distribution<> base_dist(0, 3);
    
    const char bases[4] = {'A', 'C', 'G', 'T'};
    
    for (size_t i = 0; i < numTests; i++) {
        // Get a random k-mer from the sequence
        size_t pos = pos_dist(gen);
        std::string kmer = sequence.substr(pos, kmerSize);
        
        // Mutate one random position to make it different
        size_t mut_pos = gen() % kmerSize;
        char original = kmer[mut_pos];
        char replacement;
        do {
            replacement = bases[base_dist(gen)];
        } while (replacement == original);
        
        kmer[mut_pos] = replacement;
        test_kmers.push_back(kmer);
    }
    
    return test_kmers;
}

void run_benchmark(const std::string& sequence, int kmerSize, size_t maxKmers, bool measureFalsePositives) {
    size_t numKmers = sequence.length() - kmerSize + 1;
    size_t maxKmersToProcess = std::min(numKmers, maxKmers);
    
    // Pre-calculate capacity needed (~1.5x the number of elements)
    uint32_t initialCapacity = std::max(1024u, static_cast<uint32_t>(maxKmersToProcess * 1.5));
    std::cout << "Creating filter with initial capacity: " << initialCapacity << std::endl;
    
    // Use a higher split condition (8) to reduce expansion frequency
    BambooFilter filter(initialCapacity, 8);
    
    // Store k-mers for verification
    std::vector<std::string> inserted_kmers;
    if (measureFalsePositives) {
        inserted_kmers.reserve(maxKmersToProcess);
    }
    
    // Use smaller batches for processing
    const size_t BATCH_SIZE = 1000;
    auto startTime = std::chrono::high_resolution_clock::now();
    size_t peakMemory = 0;
    
    std::cout << "Processing up to " << maxKmersToProcess << " k-mers in batches of " << BATCH_SIZE << std::endl;
    
    // Insert phase
    size_t inserted = 0;
    for (size_t i = 0; i < maxKmersToProcess; i += BATCH_SIZE) {
        size_t endBatch = std::min(i + BATCH_SIZE, maxKmersToProcess);
        
        for (size_t j = i; j < endBatch; j++) {
            if (j + kmerSize <= sequence.length()) {
                std::string kmer = sequence.substr(j, kmerSize);
                filter.Insert(kmer.c_str());
                inserted++;
                
                if (measureFalsePositives && inserted_kmers.size() < 100000) {
                    inserted_kmers.push_back(kmer);
                }
                
                // Check memory periodically
                if (j % 1000 == 0) {
                    size_t currentMemory = get_memory_usage();
                    peakMemory = std::max(peakMemory, currentMemory);
                }
            }
        }
        
        if (i % 10000 == 0) {
            std::cout << "Processed batch " << i << " to " << endBatch 
                    << " (" << std::fixed << std::setprecision(2) << (100.0 * endBatch / maxKmersToProcess) << "%)" << std::endl;
        }
    }
    
    auto insertEndTime = std::chrono::high_resolution_clock::now();
    auto insertDuration = std::chrono::duration_cast<std::chrono::milliseconds>(insertEndTime - startTime).count();
    
    // Lookup phase - check for true positives
    auto lookupStartTime = std::chrono::high_resolution_clock::now();
    size_t lookupSuccess = 0;
    size_t lookupTests = std::min(inserted, size_t(10000));  // Limit lookup tests
    
    for (size_t i = 0; i < lookupTests; i++) {
        size_t idx = (i * 97) % inserted; // Use prime number to distribute tests
        std::string kmer = sequence.substr(idx, kmerSize);
        if (filter.Lookup(kmer.c_str())) {
            lookupSuccess++;
        }
    }
    
    auto lookupEndTime = std::chrono::high_resolution_clock::now();
    auto lookupDuration = std::chrono::duration_cast<std::chrono::milliseconds>(lookupEndTime - lookupStartTime).count();
    
    // False positive testing
    double falsePositiveRate = 0.0;
    if (measureFalsePositives) {
        auto fpStartTime = std::chrono::high_resolution_clock::now();
        
        // Generate k-mers that shouldn't be in the filter
        size_t fpTests = 10000;
        std::vector<std::string> fp_test_kmers = generate_false_positive_tests(sequence, kmerSize, fpTests);
        
        std::cout << "Peak memory: " << (peakMemory / (1024.0 * 1024.0)) << " MB\n";
        
        if (filter.last_expand_before_memory > 0) {
            std::cout << "Memory before expand: " << (filter.last_expand_before_memory / (1024.0 * 1024.0)) << " MB\n";
            std::cout << "Memory after expand: " << (filter.last_expand_after_memory / (1024.0 * 1024.0)) << " MB\n";
            std::cout << "Memory growth during expand: " << ((filter.last_expand_after_memory - filter.last_expand_before_memory) / (1024.0 * 1024.0)) << " MB\n";
        }
        std::cout << "Bits per element: " << (peakMemory * 8.0 / inserted) << "\n";

        size_t falsePositives = 0;
        for (const auto& kmer : fp_test_kmers) {
            if (filter.Lookup(kmer.c_str())) {
                falsePositives++;
            }
        }
        
        falsePositiveRate = static_cast<double>(falsePositives) / fpTests;
        auto fpEndTime = std::chrono::high_resolution_clock::now();
        auto fpDuration = std::chrono::duration_cast<std::chrono::milliseconds>(fpEndTime - fpStartTime).count();
        
        std::cout << "False positive test time: " << fpDuration << " ms" << std::endl;
    }
    
    // Output results
    std::cout << "\nResults for k=" << kmerSize << ", elements=" << inserted << ":\n";
    std::cout << "---------------------------------------------\n";
    std::cout << "False positive rate: " << (falsePositiveRate * 100.0) << "%\n";
    // [Other metrics]
    std::cout << "---------------------------------------------\n";
    std::cout << "Insert time: " << insertDuration << " ms\n";
    std::cout << "Insert throughput: " << (inserted * 1000.0 / insertDuration) << " k-mers/sec\n";
    std::cout << "Lookup time: " << lookupDuration << " ms\n";
    std::cout << "Lookup throughput: " << (lookupTests * 1000.0 / lookupDuration) << " k-mers/sec\n";
    
    if (measureFalsePositives) {
        std::cout << "False positive rate: " << (falsePositiveRate * 100.0) << "%\n";
    }
    
    std::cout << "Peak memory: " << (peakMemory / (1024.0 * 1024.0)) << " MB\n";
    std::cout << "Bits per element: " << (peakMemory * 8.0 / inserted) << "\n";
    std::cout << "---------------------------------------------\n";
    
    // Append to CSV file
    std::ofstream csv("bamboo_filter_results.csv", std::ios::app);
    if (csv.is_open()) {
        // Check if file is empty to write header
        csv.seekp(0, std::ios::end);
        if (csv.tellp() == 0) {
            csv << "kmer_size,num_elements,insert_time_ms,insert_throughput,lookup_time_ms,lookup_throughput,false_positive_rate,peak_memory_mb,expand_memory_mb,bits_per_element\n";
        }
        
        csv << kmerSize << ","
            << inserted << ","
            << insertDuration << ","
            << std::fixed << std::setprecision(2) << (inserted * 1000.0 / insertDuration) << ","
            << lookupDuration << ","
            << std::fixed << std::setprecision(2) << (lookupTests * 1000.0 / lookupDuration) << ","
            << std::fixed << std::setprecision(2) << (measureFalsePositives ? (falsePositiveRate * 100.0) : 0.0) << ","
            << std::fixed << std::setprecision(4) << (peakMemory / (1024.0 * 1024.0)) << ","
            << std::fixed << std::setprecision(4) << (filter.last_expand_after_memory > 0 ? ((filter.last_expand_after_memory - filter.last_expand_before_memory) / (1024.0 * 1024.0)) : 0.0) << ","
            << std::fixed << std::setprecision(4) << (peakMemory * 8.0 / inserted) << "\n";
        
        csv.close();
    }
}

// Function for insert-only action
void run_insert_only(const std::string& sequence, int kmerSize, size_t maxKmers) {
    
    size_t numKmers = sequence.length() - kmerSize + 1;
    size_t maxKmersToProcess = std::min(numKmers, maxKmers);
    
    // Pre-calculate capacity needed (~1.5x the number of elements)
    uint32_t initialCapacity = std::max(1024u, static_cast<uint32_t>(maxKmersToProcess * 1.5));
    std::cout << "Creating filter with initial capacity: " << initialCapacity << std::endl;
    
    BambooFilter filter(initialCapacity, 8);
    
    // Use smaller batches for processing
    const size_t BATCH_SIZE = 1000;
    auto startTime = std::chrono::high_resolution_clock::now();
    size_t peakMemory = 0;
    size_t initialMemory = get_memory_usage();
    
    std::cout << "Processing up to " << maxKmersToProcess << " k-mers in batches of " << BATCH_SIZE << std::endl;
    std::cout << "Initial memory usage: " << (initialMemory / (1024.0 * 1024.0)) << " MB" << std::endl;
    
    // Insert phase
    size_t inserted = 0;
    for (size_t i = 0; i < maxKmersToProcess; i += BATCH_SIZE) {
        size_t endBatch = std::min(i + BATCH_SIZE, maxKmersToProcess);
        
        for (size_t j = i; j < endBatch; j++) {
            if (j + kmerSize <= sequence.length()) {
                std::string kmer = sequence.substr(j, kmerSize);
                filter.Insert(kmer.c_str());
                inserted++;
                
                // Check memory periodically
                if (j % 1000 == 0) {
                    size_t currentMemory = get_memory_usage();
                    peakMemory = std::max(peakMemory, currentMemory);
                }
            }
        }
        
        if (i % 10000 == 0) {
            std::cout << "Processed batch " << i << " to " << endBatch 
                    << " (" << std::fixed << std::setprecision(2) << (100.0 * endBatch / maxKmersToProcess) << "%)" << std::endl;
        }
    }
    
    auto insertEndTime = std::chrono::high_resolution_clock::now();
    auto insertDuration = std::chrono::duration_cast<std::chrono::milliseconds>(insertEndTime - startTime).count();

    std::cout << "Peak memory: " << (peakMemory / (1024.0 * 1024.0)) << " MB\n";
    if (filter.last_expand_before_memory > 0) {
        std::cout << "Memory before expand: " << (filter.last_expand_before_memory / (1024.0 * 1024.0)) << " MB\n";
        std::cout << "Memory after expand: " << (filter.last_expand_after_memory / (1024.0 * 1024.0)) << " MB\n";
        std::cout << "Memory growth during expand: " << ((filter.last_expand_after_memory - filter.last_expand_before_memory) / (1024.0 * 1024.0)) << " MB\n";
    }
    std::cout << "Memory growth: " << ((peakMemory - initialMemory) / (1024.0 * 1024.0)) << " MB\n";
    
    // Output results
    std::cout << "\nResults for k=" << kmerSize << ", elements=" << inserted << ":\n";
    std::cout << "---------------------------------------------\n";
    std::cout << "Insert time: " << insertDuration << " ms\n";
    std::cout << "Insert throughput: " << (inserted * 1000.0 / insertDuration) << " k-mers/sec\n";
    std::cout << "Initial memory: " << (initialMemory / (1024.0 * 1024.0)) << " MB\n";
    std::cout << "Peak memory: " << (peakMemory / (1024.0 * 1024.0)) << " MB\n";
    std::cout << "Memory growth: " << ((peakMemory - initialMemory) / (1024.0 * 1024.0)) << " MB\n";
    std::cout << "Bits per element: " << (peakMemory * 8.0 / inserted) << "\n";
    std::cout << "---------------------------------------------\n";
    
    // Append to CSV file
    std::ofstream csv("bamboo_filter_results.csv", std::ios::app);
    if (csv.is_open()) {
        // Check if file is empty to write header
        csv.seekp(0, std::ios::end);
        if (csv.tellp() == 0) {
            csv << "action,kmer_size,num_elements,insert_time_ms,insert_throughput,peak_memory_mb,bits_per_element\n";
        }
        
        csv << "insert,"
            << kmerSize << ","
            << inserted << ","
            << insertDuration << ","
            << (inserted * 1000.0 / insertDuration) << ","
            << (peakMemory / (1024.0 * 1024.0)) << ","
            << (peakMemory * 8.0 / inserted) << "\n";
        
        csv.close();
    }
}

// Function for lookup-only action
void run_lookup_only(const std::string& sequence, int kmerSize, size_t maxKmers) {
    size_t numKmers = sequence.length() - kmerSize + 1;
    size_t maxKmersToProcess = std::min(numKmers, maxKmers);
    
    // Pre-calculate capacity needed (~1.5x the number of elements)
    uint32_t initialCapacity = std::max(1024u, static_cast<uint32_t>(maxKmersToProcess * 1.5));
    std::cout << "Creating filter with initial capacity: " << initialCapacity << std::endl;
    
    BambooFilter filter(initialCapacity, 8);
    
    // Use smaller batches for processing
    const size_t BATCH_SIZE = 1000;
    
    std::cout << "Processing up to " << maxKmersToProcess << " k-mers in batches of " << BATCH_SIZE << std::endl;
    
    // Insert phase (without timing)
    size_t inserted = 0;
    for (size_t i = 0; i < maxKmersToProcess; i += BATCH_SIZE) {
        size_t endBatch = std::min(i + BATCH_SIZE, maxKmersToProcess);
        
        for (size_t j = i; j < endBatch; j++) {
            if (j + kmerSize <= sequence.length()) {
                std::string kmer = sequence.substr(j, kmerSize);
                filter.Insert(kmer.c_str());
                inserted++;
            }
        }
        
        if (i % 10000 == 0) {
            std::cout << "Inserted batch " << i << " to " << endBatch 
                    << " (" << std::fixed << std::setprecision(2) << (100.0 * endBatch / maxKmersToProcess) << "%)" << std::endl;
        }
    }
    
    std::cout << "Inserted " << inserted << " k-mers. Starting lookup tests..." << std::endl;
    
    // Lookup phase - check for true positives
    auto lookupStartTime = std::chrono::high_resolution_clock::now();
    size_t lookupSuccess = 0;
    size_t lookupTests = std::min(inserted, size_t(100000));  // More lookup tests for lookup-focused test
    
    for (size_t i = 0; i < lookupTests; i++) {
        size_t idx = (i * 97) % inserted; // Use prime number to distribute tests
        std::string kmer = sequence.substr(idx, kmerSize);
        if (filter.Lookup(kmer.c_str())) {
            lookupSuccess++;
        }
        
        if (i % 10000 == 0 && i > 0) {
            std::cout << "Performed " << i << " lookups (" 
                    << std::fixed << std::setprecision(2) << (100.0 * i / lookupTests) << "%)" << std::endl;
        }
    }
    
    auto lookupEndTime = std::chrono::high_resolution_clock::now();
    auto lookupDuration = std::chrono::duration_cast<std::chrono::milliseconds>(lookupEndTime - lookupStartTime).count();
    
    // Generate k-mers that shouldn't be in the filter
    size_t fpTests = 10000;
    std::vector<std::string> fp_test_kmers = generate_false_positive_tests(sequence, kmerSize, fpTests);
    
    auto fpStartTime = std::chrono::high_resolution_clock::now();
    size_t falsePositives = 0;
    for (const auto& kmer : fp_test_kmers) {
        if (filter.Lookup(kmer.c_str())) {
            falsePositives++;
        }
    }
    
    double falsePositiveRate = static_cast<double>(falsePositives) / fpTests;
    auto fpEndTime = std::chrono::high_resolution_clock::now();
    auto fpDuration = std::chrono::duration_cast<std::chrono::milliseconds>(fpEndTime - fpStartTime).count();
    
    // Output results
    std::cout << "\nResults for k=" << kmerSize << ", elements=" << inserted << ":\n";
    std::cout << "---------------------------------------------\n";
    std::cout << "Lookup time (for " << lookupTests << " lookups): " << lookupDuration << " ms\n";
    std::cout << "Lookup throughput: " << (lookupTests * 1000.0 / lookupDuration) << " lookups/sec\n";
    std::cout << "False positive rate: " << (falsePositiveRate * 100.0) << "%\n";
    std::cout << "False positive test time: " << fpDuration << " ms\n";
    std::cout << "---------------------------------------------\n";
    
    // Append to CSV file
    std::ofstream csv("bamboo_filter_results.csv", std::ios::app);
    if (csv.is_open()) {
        // Check if file is empty to write header
        csv.seekp(0, std::ios::end);
        if (csv.tellp() == 0) {
            csv << "action,kmer_size,num_elements,lookup_time_ms,lookup_throughput,false_positive_rate\n";
        }
        
        csv << "lookup,"
            << kmerSize << ","
            << inserted << ","
            << lookupDuration << ","
            << (lookupTests * 1000.0 / lookupDuration) << ","
            << (lookupSuccess * 100.0 / lookupTests) << ","
            << (falsePositiveRate * 100.0) << "\n";
        
        csv.close();
    }
}

// Function for delete-only action
void run_delete_only(const std::string& sequence, int kmerSize, size_t maxKmers) {
    size_t numKmers = sequence.length() - kmerSize + 1;
    size_t maxKmersToProcess = std::min(numKmers, maxKmers);
    
    // Pre-calculate capacity needed (~1.5x the number of elements)
    uint32_t initialCapacity = std::max(1024u, static_cast<uint32_t>(maxKmersToProcess * 1.5));
    std::cout << "Creating filter with initial capacity: " << initialCapacity << std::endl;
    
    BambooFilter filter(initialCapacity, 8);
    
    // Store all k-mers for verification
    std::vector<std::string> kmers;
    kmers.reserve(maxKmersToProcess);
    
    // Use smaller batches for processing
    const size_t BATCH_SIZE = 1000;
    
    std::cout << "Processing up to " << maxKmersToProcess << " k-mers in batches of " << BATCH_SIZE << std::endl;
    
    // Insert phase (without timing)
    size_t inserted = 0;
    for (size_t i = 0; i < maxKmersToProcess; i += BATCH_SIZE) {
        size_t endBatch = std::min(i + BATCH_SIZE, maxKmersToProcess);
        
        for (size_t j = i; j < endBatch; j++) {
            if (j + kmerSize <= sequence.length()) {
                std::string kmer = sequence.substr(j, kmerSize);
                filter.Insert(kmer.c_str());
                kmers.push_back(kmer);
                inserted++;
            }
        }
        
        if (i % 10000 == 0) {
            std::cout << "Inserted batch " << i << " to " << endBatch 
                    << " (" << std::fixed << std::setprecision(2) << (100.0 * endBatch / maxKmersToProcess) << "%)" << std::endl;
        }
    }
    
    std::cout << "Inserted " << inserted << " k-mers. Starting delete tests..." << std::endl;
    
    // Delete phase
    auto deleteStartTime = std::chrono::high_resolution_clock::now();
    size_t deleteSuccess = 0;
    
    for (size_t i = 0; i < kmers.size(); i++) {
        const std::string& kmer = kmers[i];
        filter.Delete(kmer.c_str());
        
        // Verify deletion every 100 items to avoid slowdown from excessive verification
        if (i % 100 == 0) {
            if (!filter.Lookup(kmer.c_str())) {
                deleteSuccess++;
            }
        }
        
        if (i % 10000 == 0 && i > 0) {
            std::cout << "Deleted " << i << " k-mers (" 
                    << std::fixed << std::setprecision(2) << (100.0 * i / kmers.size()) << "%)" << std::endl;
        }
    }
    
    auto deleteEndTime = std::chrono::high_resolution_clock::now();
    auto deleteDuration = std::chrono::duration_cast<std::chrono::milliseconds>(deleteEndTime - deleteStartTime).count();
    
    // Calculate delete success rate
    size_t verificationCount = kmers.size() / 100 + 1;
    double deleteSuccessRate = (double)deleteSuccess / verificationCount * 100.0;
    
    // Output results
    std::cout << "\nResults for k=" << kmerSize << ", elements=" << inserted << ":\n";
    std::cout << "---------------------------------------------\n";
    std::cout << "Delete time: " << deleteDuration << " ms\n";
    std::cout << "Delete throughput: " << (kmers.size() * 1000.0 / deleteDuration) << " deletes/sec\n";
    std::cout << "Delete success rate: " << deleteSuccessRate << "%\n";
    std::cout << "---------------------------------------------\n";
    
    // Append to CSV file
    std::ofstream csv("bamboo_filter_results.csv", std::ios::app);
    if (csv.is_open()) {
        // Check if file is empty to write header
        csv.seekp(0, std::ios::end);
        if (csv.tellp() == 0) {
            csv << "action,kmer_size,num_elements,delete_time_ms,delete_throughput,delete_success_rate\n";
        }
        
        csv << "delete,"
            << kmerSize << ","
            << inserted << ","
            << deleteDuration << ","
            << (kmers.size() * 1000.0 / deleteDuration) << ","
            << deleteSuccessRate << "\n";
        
        csv.close();
    }
}

int main(int argc, char* argv[]) {
    std::string fastaPath;
    int kmerSize = 25;
    std::string action = "benchmark";
    size_t maxProcessKmers = 10000000; 
    bool limitReadSize = false;
    size_t maxReadSize = 0;
    bool autoTest = false;
    bool measureFalsePositives = false;
    
    // Parse command line args
    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "--fasta-path" && i + 1 < argc) {
            fastaPath = argv[++i];
        } else if (arg == "--kmer-size" && i + 1 < argc) {
            kmerSize = std::stoi(argv[++i]);
        } else if (arg == "--action" && i + 1 < argc) {
            action = argv[++i];
        } else if (arg == "--max-kmers" && i + 1 < argc) {
            maxProcessKmers = std::stoul(argv[++i]);
        } else if (arg == "--limit-read" && i + 1 < argc) {
            limitReadSize = true;
            maxReadSize = std::stoul(argv[++i]);
        } else if (arg == "--auto-test") {
            autoTest = true;
        } else if (arg == "--false-positives") {
            measureFalsePositives = true;
        }
    }
    
    if (fastaPath.empty()) {
        std::cerr << "Error: FASTA path is required\n";
        std::cerr << "Usage: " << argv[0] << " --fasta-path <path> [options]\n";
        std::cerr << "Options:\n";
        std::cerr << "  --kmer-size <k>            : k-mer size (default: 25)\n";
        std::cerr << "  --action <action>          : benchmark, insert, lookup, delete (default: benchmark)\n";
        std::cerr << "  --max-kmers <count>        : Maximum k-mers to process (default: 10000000)\n";
        std::cerr << "  --limit-read <size>        : Limit genome read size in base pairs\n";
        std::cerr << "  --auto-test                : Run tests with different k-mer sizes (15, 20, 25, 30)\n";
        std::cerr << "  --false-positives          : Measure false positive rate\n";
        return 1;
    }
    
    try {
        // Read the genome file
        std::string sequence;
        std::ifstream in(fastaPath);
        if (!in) throw std::runtime_error("Cannot open FASTA: " + fastaPath);
        
        std::string line;
        size_t charsRead = 0;
        
        std::cout << "Reading genome file..." << std::endl;
        
        while (std::getline(in, line)) {
            if (line.empty() || line[0] == '>') continue;
            
            for (char c : line) {
                if (limitReadSize && charsRead >= maxReadSize) break;
                
                switch (std::toupper(c)) {
                    case 'A': case 'C': case 'G': case 'T': 
                        sequence.push_back(std::toupper(c)); 
                        charsRead++;
                        break;
                    default: break;
                }
            }
            
            if (limitReadSize && charsRead >= maxReadSize) break;
            
            // Show progress for large files
            if (sequence.length() % 1000000 == 0) {
                std::cout << "Read " << (sequence.length() / 1000000) << " million base pairs...\r" << std::flush;
            }
        }
        
        std::cout << "Read " << sequence.length() << " base pairs from genome" << std::endl;
        
        if (autoTest) {
        // Run benchmarks with various k-mer sizes
        std::vector<int> kSizes = {15, 20, 25, 30};
        for (int k : kSizes) {
            std::cout << "\n=== Running benchmark with k=" << k << " ===\n" << std::endl;
            
            if (action == "insert") {
                run_insert_only(sequence, k, maxProcessKmers);
            } else if (action == "lookup") {
                run_lookup_only(sequence, k, maxProcessKmers);
            } else if (action == "delete") {
                run_delete_only(sequence, k, maxProcessKmers);
            } else {
                // Default is benchmark
                run_benchmark(sequence, k, maxProcessKmers, measureFalsePositives);
            }
        }
} else {
    // Run a single test with specified parameters
    std::cout << "k=" << kmerSize << " | kmers=" << (sequence.length() - kmerSize + 1) << std::endl;
    
    if (action == "insert") {
        run_insert_only(sequence, kmerSize, maxProcessKmers);
    } else if (action == "lookup") {
        run_lookup_only(sequence, kmerSize, maxProcessKmers);
    } else if (action == "delete") {
        run_delete_only(sequence, kmerSize, maxProcessKmers);
    } else {
        // Default is benchmark
        run_benchmark(sequence, kmerSize, maxProcessKmers, measureFalsePositives);
    }
}
        
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}