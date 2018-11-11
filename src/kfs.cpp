#include "./kfs.hpp"

int
main(int argc, char** argv)
{
  kfs::KFS kfs;
  kfs.initialize_command_line_options(argc, argv);
  kfs.parse_input();
  kfs.reread_input_and_query_table();
}

void
kfs::KFS::parse_input(void)
{
  char* buf = NULL;
  size_t buf_len = 0;
  ssize_t buf_read = 0;
  char header_str[LINE_MAX] = {0};
  char* sequence_str = NULL;
  ssize_t sequence_read = 0;
  char* sequence_intermediate_buf = NULL;
  ssize_t line_count = 0;
  long record_count = 0;

  // initialize timing
  std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

  this->initialize_in_stream();

  sequence_str = (char*) malloc(KFS_LINE_MAX);
  if (!sequence_str) {
    std::fprintf(stderr, "Error: Could not allocate memory for sequence buffer\n");
    std::exit(ENOMEM);
  }

  sequence_intermediate_buf = (char*) malloc(KFS_LINE_MAX);
  if (!sequence_intermediate_buf) {
    std::fprintf(stderr, "Error: Could not allocate memory for sequence intermediate buffer\n");
    std::exit(ENOMEM);
  }

  while ((buf_read = getline(&buf, &buf_len, this->in_stream())) != EOF) {
    if (buf[0] == '>') {
      if ((strlen(header_str) > 0) && (strlen(sequence_str) > 0)) {
        this->process_fasta_record(header_str, sequence_str);
        sequence_read = 0;
        record_count++;
      }
      // read in next header
      std::sscanf(buf, ">%[^\n\t]\n", header_str);
    }
    else {
      std::sscanf(buf, "%s\n", sequence_intermediate_buf);
      std::memcpy(sequence_str + sequence_read, sequence_intermediate_buf, buf_read);
      sequence_read += (buf_read - 1);
      if (sequence_read >= KFS_LINE_MAX) {
        std::fprintf(stderr, "Error: Input FASTA record is longer than KFS_LINE_MAX buffer\n");
        std::exit(EOVERFLOW);
      }
    }
    line_count += 1;
  }

  // process final record
  if ((strlen(header_str) > 0) && (strlen(sequence_str) > 0)) {
    this->process_fasta_record(header_str, sequence_str);
    sequence_read = 0;
    record_count++;
  }

  // cleanup
  free(buf);
  free(sequence_str);
  free(sequence_intermediate_buf);
  this->close_in_stream();

  // conclude timing
  std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
  std::fprintf(stderr, "-finished %ld seqs.\n", record_count);
  std::fprintf(stderr, "%lld seconds to load the array.\n", static_cast<long long int>(time_span.count()));
}

void
kfs::KFS::process_fasta_record(char* header, char* sequence)
{
  std::string seq(sequence);
  std::string n("N");

  std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
  std::deque<char> window(seq.begin(), seq.begin() + this->k());
  // walk over all windows across sequence
  for (size_t i = this->k(); i <= seq.length(); ++i) {
    std::string mer_f(window.begin(), window.end());
    window.pop_front();
    window.push_back(seq[i]);
    std::size_t n_found = mer_f.find(n);
    if (n_found != std::string::npos) {
      continue;
    }
    std::string mer_r(mer_f);
    reverse_complement_string(mer_r);
    bool mer_f_is_canonical = (mer_f.compare(mer_r) <= 0);
    long mer_key = this->mer_to_key((mer_f_is_canonical) ? mer_f : mer_r);
    if (!mer_count(mer_key)) {
      set_mer_count(mer_key, 1);
    }
    else {
      increment_mer_count(mer_key);
    }
  }
}

void
kfs::KFS::reread_input_and_query_table(void)
{
  char* buf = NULL;
  size_t buf_len = 0;
  ssize_t buf_read = 0;
  char header_str[LINE_MAX] = {0};
  char* sequence_str = NULL;
  ssize_t sequence_read = 0;
  char* sequence_intermediate_buf = NULL;
  ssize_t line_count = 0;
  long record_count = 0;

  // initialize timing
  std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

  this->initialize_in_stream();

  sequence_str = (char*) malloc(KFS_LINE_MAX);
  if (!sequence_str) {
    std::fprintf(stderr, "Error: Could not allocate memory for sequence buffer\n");
    std::exit(ENOMEM);
  }

  sequence_intermediate_buf = (char*) malloc(KFS_LINE_MAX);
  if (!sequence_intermediate_buf) {
    std::fprintf(stderr, "Error: Could not allocate memory for sequence intermediate buffer\n");
    std::exit(ENOMEM);
  }

  while ((buf_read = getline(&buf, &buf_len, this->in_stream())) != EOF) {
    if (buf[0] == '>') {
      if ((strlen(header_str) > 0) && (strlen(sequence_str) > 0)) {
        this->query_fasta_sequence(sequence_str);
        sequence_read = 0;
        record_count++;
      }
      // read in next header
      std::sscanf(buf, ">%[^\n\t]\n", header_str);
    }
    else {
      std::sscanf(buf, "%s\n", sequence_intermediate_buf);
      std::memcpy(sequence_str + sequence_read, sequence_intermediate_buf, buf_read);
      sequence_read += (buf_read - 1);
      if (sequence_read >= KFS_LINE_MAX) {
        std::fprintf(stderr, "Error: Input FASTA record is longer than KFS_LINE_MAX buffer\n");
        std::exit(EOVERFLOW);
      }
    }
    line_count += 1;
  }

  // process final record
  if ((strlen(header_str) > 0) && (strlen(sequence_str) > 0)) {
    this->query_fasta_sequence(sequence_str);
    sequence_read = 0;
    record_count++;
  }

  // cleanup
  free(buf);
  free(sequence_str);
  free(sequence_intermediate_buf);
  this->close_in_stream();

  // conclude timing
  std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
  std::fprintf(stderr, "-looked at %ld seqs.\n", record_count);
  std::fprintf(stderr, "%lld seconds to lookup the kmers.\n", static_cast<long long int>(time_span.count()));
}

void
kfs::KFS::query_fasta_sequence(char* sequence)
{
  std::string seq(sequence);
  std::string n("N");

  std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
  std::deque<char> window(seq.begin(), seq.begin() + this->k());
  // walk over all windows across sequence
  for (size_t i = this->k(); i <= seq.length(); ++i) {
    std::string mer_f(window.begin(), window.end());
    window.pop_front();
    window.push_back(seq[i]);
    std::size_t n_found = mer_f.find(n);
    if (n_found != std::string::npos) {
      continue;
    }
    std::string mer_r(mer_f);
    reverse_complement_string(mer_r);
    bool mer_f_is_canonical = (mer_f.compare(mer_r) <= 0);
    long mer_key = this->mer_to_key((mer_f_is_canonical) ? mer_f : mer_r);

    // count result
    int mer_count_result = mer_count(mer_key);

    // dumb hack to disable -Wunused-variable compile warning
    (void)mer_count_result;
  }
}

void
kfs::KFS::debug_mer_counts(void)
{
  std::vector<int> v;
  for(absl::flat_hash_map<long, int>::iterator it = mer_counts().begin(); it != mer_counts().end(); ++it) {
    v.push_back(it->first);
    std::fprintf(stderr, "%s\t%d\n", this->key_to_mer(it->first).c_str(), it->second);
  }
}

long
kfs::KFS::mer_to_key(std::string &s)
{
  long key = 0;
  for (int j = this->k() - 1; j > -1; --j) {
    key += fmap[s[j]] * 1 << (2*(this->k() - 1 - j));
  }
  return key;
}

std::string
kfs::KFS::key_to_mer(long key)
{
  std::string s;
  if (this->k() == 1) {
    s.append(1, rmap[key]);
  }
  else {
    for (int j = this->k() - 1; j > -1; --j) {
      int x = 0;
      while (key >= pow(4,j)) {
        key -= pow(4,j);
        x++;
      }
      s.append(1, rmap[x]);
    }
  }
  return s;
}

void
kfs::KFS::initialize_command_line_options(int argc, char** argv)
{
  int client_long_index;
  int client_opt = getopt_long(argc,
                               argv,
                               this->client_kfs_opt_string().c_str(),
                               this->client_kfs_long_options(),
                               &client_long_index);
  int k = -1;

  opterr = 0; /* disable error reporting by GNU getopt */

  while (client_opt != -1) {
    switch (client_opt) {
    case 'k':
      std::sscanf(optarg, "%d", &k);
      this->k(k);
      break;
    case 'i':
      this->in_fn(optarg);
      break;
    case 'h':
      this->print_usage(stdout);
      std::exit(EXIT_SUCCESS);
    case 'v':
      this->print_version(stdout);
      std::exit(EXIT_SUCCESS);
    case '?':
      this->print_usage(stdout);
      std::exit(EXIT_SUCCESS);
    default:
      break;
    }
    client_opt = getopt_long(argc,
                             argv,
                             this->client_kfs_opt_string().c_str(),
                             this->client_kfs_long_options(),
                             &client_long_index);
  }

  bool error_flagged = false;
  
  if (this->k() == -1) {
    std::fprintf(stderr, "Error: Specify k value\n");
    error_flagged = true;
  }

  if (this->in_fn().length() == 0) {
    std::fprintf(stderr, "Error: Specify input filename value\n");
    error_flagged = true;
  }

  if (error_flagged) {
    this->print_usage(stderr);
    std::exit(ENODATA);
  }
}

std::string
kfs::KFS::client_kfs_opt_string(void)
{
  std::string _s("k:i:hv?");
  return _s;
}

struct option*
kfs::KFS::client_kfs_long_options(void)
{
  static struct option _k = { "k",              required_argument,   NULL,    'k' };
  static struct option _i = { "i",              required_argument,   NULL,    'k' };
  static struct option _h = { "help",           no_argument,         NULL,    'h' };
  static struct option _v = { "version",        no_argument,         NULL,    'v' };
  static struct option _0 = { NULL,             no_argument,         NULL,     0  };
  std::vector<struct option> _s;
  _s.push_back(_k);
  _s.push_back(_i);
  _s.push_back(_h);
  _s.push_back(_v);
  _s.push_back(_0);
  return &_s[0];
}

std::string
kfs::KFS::client_kfs_name(void)
{
  std::string _s("kfs");
  return _s;
}

std::string
kfs::KFS::client_kfs_version(void)
{
  std::string _s("0.0.1");
  return _s;
}

std::string
kfs::KFS::client_kfs_authors(void)
{
  std::string _s("Alex Reynolds");
  return _s;
}

std::string
kfs::KFS::client_kfs_usage(void)
{
  std::string _s("\n"				\
		 "  Usage:\n"			\
		 "\n"						\
		 "  $ kfs [arguments] < input\n");
  return _s;
}

void
kfs::KFS::print_usage(FILE* os)
{
  this->print_version(os);
}

void
kfs::KFS::print_version(FILE* os)
{
  std::fprintf(os,
               "%s\n"     \
               "  version: %s\n"     \
               "  author:  %s\n",
               this->client_kfs_name().c_str(),
               this->client_kfs_version().c_str(),
               this->client_kfs_authors().c_str());
}
