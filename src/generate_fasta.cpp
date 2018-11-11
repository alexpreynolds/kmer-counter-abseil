//generate a fasta file to count k-mers
#include <iostream>
#include <random>

char gen_base(int q);

char gen_base(int q){
  if (q <= 30){
    return 'A';
  } else if ((q > 30) && (q <=60) ){
    return 'T';
  } else if ((q > 60) && (q <=80) ){
    return 'C';
  } else if (q > 80){
    return 'G';
  }
  return 'N';
}

int main() {
  unsigned seed = 1;
  std::default_random_engine generator (seed);
  std::poisson_distribution<int> poisson (59);
  std::geometric_distribution<int> geo (0.05);
  std::uniform_int_distribution<int> uniform (1,100);
  int printval;
  int i=0;
  while(i<2000000){
    if (i % 2 == 0){
      printval = poisson(generator);
    } else {
      printval = geo(generator) + 29;
    }
    if (printval >= 29){
      std::cout << '>' << i << '\n';
      //std::cout << printval << '\n';
      for (int j = 0; j < printval; j++){
        std::cout << gen_base(uniform(generator));
      }
      std::cout << '\n';
      i++;
    }
  }
  return 0;
}
