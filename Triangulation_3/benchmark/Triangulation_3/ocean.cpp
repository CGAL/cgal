#include <fstream>
#include <CGAL/Random.h>

int main()
{
  int N=100;
  std::cout.precision(17);
  CGAL::Random rng;

  for(int i = 0; i < N; ++i){
    double x = rng.get_double(-1.0, 1.0);
    double y = rng.get_double(-1.0, 1.0);

    for(int j = 0; j < N; j++){
      std::cout << x << " " << y  <<  " " << rng.get_double(-1.0, 1.0) << std::endl;
    }
  }

  return 0;
}
