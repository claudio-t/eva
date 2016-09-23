#include <iostream>
#include <armadillo>
#include <array>


int main(int argc, char * argv[]) {
    
  arma::mat A = arma::randu<arma::mat>(4,5);
  arma::mat B = arma::randu<arma::mat>(4,5);
  
  std::cout << A*B.t() << std::endl;
  
  constexpr static int N = 2;
  std::array<double, N*N> arr;
  return 0;
    
}
