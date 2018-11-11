#ifndef KFS_H_
#define KFS_H_

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif /* getline() support */

#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#include <string>
#include <vector>
#include <deque>
#include <map>
#include <unordered_map>
#include <algorithm>
#include <iostream>
#include <exception>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <getopt.h>
#include <system_error>
#include <chrono>

#include "absl/container/flat_hash_map.h"

#define KFS_LINE_MAX 268435456

namespace kfs
{
  class KFS
  {
  private:
    int _k;
    std::string _in_fn;
    FILE* _in_stream;
    absl::flat_hash_map<long, int> _mer_counts;

  public:
    void parse_input(void);
    void reread_input_and_query_table(void);
    void query_fasta_sequence(char* sequence);
    void process_fasta_record(char* header, char* sequence);
    void debug_mer_counts(void);
    long mer_to_key(std::string &s);
    std::string key_to_mer(long key);

    void initialize_command_line_options(int argc, char** argv);
    void print_usage(FILE* os);
    void print_version(FILE* os);
    std::string client_kfs_opt_string(void);
    struct option* client_kfs_long_options(void);
    std::string client_kfs_name(void);
    std::string client_kfs_version(void);
    std::string client_kfs_authors(void);
    std::string client_kfs_usage(void);

    int k() const { return _k; }
    void k(const int& k) { _k = k; }
    std::string in_fn() const { return _in_fn; }
    void in_fn(const std::string& ifn) { _in_fn = ifn; }

    inline absl::flat_hash_map<long, int>& mer_counts(void) { return _mer_counts; }
    inline int mer_count(const long& k) { return _mer_counts[k]; }
    inline void mer_counts(const absl::flat_hash_map<long, int>& mc) { _mer_counts = mc; }
    inline void set_mer_count(const long& k, const int& v) { _mer_counts[k] = v; }
    inline void erase_mer_count(const long& k) { _mer_counts.erase(k); }
    inline void increment_mer_count(const long& k) { _mer_counts[k]++; }

    std::map<unsigned char, int> create_fmap(void) {
      std::map<unsigned char, int> m;
      m['A'] = 0;
      m['C'] = 1;
      m['T'] = 2;
      m['G'] = 3;
      return m;
    } 
    std::map<unsigned char, int> fmap = create_fmap();

    std::map<int, unsigned char> create_rmap(void) {
      std::map<int, unsigned char> m;
      m[0] = 'A';
      m[1] = 'C';
      m[2] = 'T';
      m[3] = 'G';
      return m;
    }
    std::map<int, unsigned char> rmap = create_rmap();

    inline static void reverse_complement_string(std::string &s) {
      static unsigned char base_complement_map[256] = {
          0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,  13,  14,  15,
         16,  17,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28,  29,  30,  31,
         32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  43,  44,  45,  46,  47,
         48,  49,  50,  51,  52,  53,  54,  55,  56,  57,  58,  59,  60,  61,  62,  63,
         64, 'T', 'V', 'G', 'H', 'E', 'F', 'C', 'D', 'I', 'J', 'M', 'L', 'K', 'N', 'O',
        'P', 'Q', 'Y', 'S', 'A', 'A', 'B', 'W', 'X', 'R', 'Z',  91,  92,  93,  94,  95,
         64, 't', 'v', 'g', 'h', 'e', 'f', 'c', 'd', 'i', 'j', 'm', 'l', 'k', 'n', 'o',
        'p', 'q', 'y', 's', 'a', 'a', 'b', 'w', 'x', 'r', 'z', 123, 124, 125, 126, 127,
        128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143,
        144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159,
        160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175,
        176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191,
        192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207,
        208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223,
        224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239,
        240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255
      };
      std::reverse(s.begin(), s.end());
      for (auto i = s.begin(); i != s.end(); ++i) {
        *i = base_complement_map[(int) *i];
      }
    }

    FILE* in_stream(void) { return _in_stream; }
    void in_stream(FILE** isp) { _in_stream = *isp; }
    void initialize_in_stream(void) {
      FILE* in_fp = NULL;
      in_fp = this->in_fn().empty() ? stdin : std::fopen(this->in_fn().c_str(), "r");
      if (!in_fp) {
        std::fprintf(stderr, "Error: Input file handle could not be created\n");
        std::exit(ENODATA); /* No message is available on the STREAM head read queue (POSIX.1) */
      }
      this->in_stream(&in_fp);
    }
    void close_in_stream(void) {
      std::fclose(this->in_stream());
    }

    KFS();
    ~KFS();
  };

  KFS::KFS() {
    k(-1);
    in_fn("");
  }
  
  KFS::~KFS() {
  }
}

#endif // KFS_H_
