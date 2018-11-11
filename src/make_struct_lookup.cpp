#include <fstream>
#include <string>
#include <map>
#include <iostream>
#include <chrono>
//build the structure. measure how much RAM it consumes.
//then measure how long it takes to lookup in the data structure

#define k 29

std::string rc(std::string seq){
  std::string rc;
  for (int i = seq.length()-1; i>=0; i--){
    if (seq[i] == 'A'){
      rc.push_back('T');
    } else if (seq[i] == 'C'){
      rc.push_back('G');
    } else if (seq[i] == 'G'){
      rc.push_back('C');
    } else if (seq[i] == 'T'){
      rc.push_back('A');
    }
  }
  return rc;
}

int main(int argc, char* argv[]){
  using namespace std::chrono;
  //initialize the data structure
  std::string thisline;
  std::map<std::string, int> kmer_map;
  std::string header;
  std::string seq;
  //open the fasta file
  std::ifstream inFile;
  inFile.open(argv[1]);

  //construct the kmer-lookup structure
  int i = 0;
  high_resolution_clock::time_point t1 = high_resolution_clock::now();
  while (getline(inFile,thisline)){
    if (thisline[0] == '>'){
      header = thisline.substr(1,thisline.size());
      //std::cout << header << '\n';
    } else {
      seq = thisline;
      //now add the kmers
      for (int j=0; j< thisline.size() - k + 1; j++){
        kmer_map[seq.substr(j,j+k)] = stoi(header);
        kmer_map[rc(seq.substr(j,j+k))] = stoi(header);
      }
      i++;
    }
  }
  std::cout << "  -finished " << i << " seqs.\n";
  inFile.close();
  high_resolution_clock::time_point t2 = high_resolution_clock::now();
  duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
  std::cout << time_span.count() << " seconds to load the array." << '\n';

  //now lookup the kmers
  inFile.open(argv[1]);
  t1 = high_resolution_clock::now();
  int lookup;
  while (getline(inFile,thisline)){
    if (thisline[0] != '>'){
      seq = thisline;
      //now lookup the kmers
      for (int j=0; j< thisline.size() - k + 1; j++){
        lookup = kmer_map[seq.substr(j,j+k)];
      }
    }
  }
  std::cout << "  - looked at " << i << " seqs.\n";
  inFile.close();
  t2 = high_resolution_clock::now();
  time_span = duration_cast<duration<double>>(t2 - t1);
  std::cout << time_span.count() << " seconds to lookup the kmers." << '\n';

}