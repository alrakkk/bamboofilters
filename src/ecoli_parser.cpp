/* src/ecoli_parser.cpp */
#include "ecoli_parser.h"
#include <cctype>
#include <fstream>
#include <stdexcept>

// ---------- FASTA reader -------------------------------------------------
size_t read_fasta_sequence(const std::string& path, std::string& seq)
{
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot open FASTA: " + path);
    seq.clear();
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '>') continue;
        for (char c : line) {
            switch (std::toupper(c)) {
                case 'A': case 'C': case 'G': case 'T': seq.push_back(std::toupper(c)); break;
                default: break;
            }
        }
    }
    return seq.size();
}

// ---------- MurmurHash3-x64 (low 64 bits) -------------------------------
static inline uint64_t rotl(uint64_t x, int r){return (x<<r)|(x>>(64-r));}
static inline uint64_t fmix(uint64_t k){
    k^=k>>33; k*=0xff51afd7ed558ccdULL;
    k^=k>>33; k*=0xc4ceb9fe1a85ec53ULL;
    k^=k>>33; return k;
}
static uint64_t murmur64low(const void* key,int len,uint32_t seed=42){
    const uint8_t* d=(const uint8_t*)key; int n=len/16;
    uint64_t h1=seed,h2=seed,c1=0x87c37b91114253d5ULL,c2=0x4cf5ad432745937fULL;
    const uint64_t* blocks=(const uint64_t*)d;
    for(int i=0;i<n;++i){
        uint64_t k1=blocks[i*2],k2=blocks[i*2+1];
        k1*=c1; k1=rotl(k1,31); k1*=c2; h1^=k1;
        h1=rotl(h1,27); h1+=h2; h1=h1*5+0x52dce729;
        k2*=c2; k2=rotl(k2,33); k2*=c1; h2^=k2;
        h2=rotl(h2,31); h2+=h1; h2=h2*5+0x38495ab5;
    }
    const uint8_t* tail=d+n*16; uint64_t k1=0,k2=0;
    switch(len&15){
        case15:k2^=uint64_t(tail[14])<<48;
        case14:k2^=uint64_t(tail[13])<<40;
        case13:k2^=uint64_t(tail[12])<<32;
        case12:k2^=uint64_t(tail[11])<<24;
        case11:k2^=uint64_t(tail[10])<<16;
        case10:k2^=uint64_t(tail[9])<<8;
        case9 :k2^=uint64_t(tail[8]);
               k2*=c2;k2=rotl(k2,33);k2*=c1;h2^=k2;
        case8 :k1^=uint64_t(tail[7])<<56;
        case7 :k1^=uint64_t(tail[6])<<48;
        case6 :k1^=uint64_t(tail[5])<<40;
        case5 :k1^=uint64_t(tail[4])<<32;
        case4 :k1^=uint64_t(tail[3])<<24;
        case3 :k1^=uint64_t(tail[2])<<16;
        case2 :k1^=uint64_t(tail[1])<<8;
        case1 :k1^=uint64_t(tail[0]);
               k1*=c1;k1=rotl(k1,31);k1*=c2;h1^=k1;
    }
    h1^=len; h2^=len; h1+=h2; h2+=h1;
    h1=fmix(h1); h2=fmix(h2); h1+=h2;
    return h1;
}

// ---------- fast 2-bit rolling hash (k ≤ 32) ----------------------------
static inline uint64_t enc(char nt){
    return (nt=='C'||nt=='c')?1ULL:(nt=='G'||nt=='g')?2ULL:(nt=='T'||nt=='t')?3ULL:0ULL;
}

std::vector<uint64_t> hash_kmers(const std::string& s, uint32_t k)
{
    const size_t n=s.size();
    if(k==0||n<k) return {};
    std::vector<uint64_t> out; out.reserve(n-k+1);

    if(k<=32){                               // 2-bit packed
        uint64_t mask=(k==32)?~0ULL:((1ULL<<(2*k))-1);
        uint64_t code=0;
        for(uint32_t i=0;i<k;++i) code=(code<<2)|enc(s[i]);
        out.push_back(code);
        for(size_t i=k;i<n;++i){
            code=((code<<2)|enc(s[i]))&mask;
            out.push_back(code);
        }
    }else{                                   // large k → Murmur
        for(size_t i=0;i+k<=n;++i)
            out.push_back(murmur64low(&s[i],k));
    }
    return out;
}
