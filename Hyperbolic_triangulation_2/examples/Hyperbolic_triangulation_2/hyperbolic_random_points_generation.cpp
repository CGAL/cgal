#include <fstream>
#include <string>
#include <sstream>

// CGAL headers

#include <CGAL/Hyperbolic_random_points_in_disc_2.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef K::Point_2 Point_2;
typedef K::FT FT;

using namespace std;

int main()
{
  typedef CGAL::Creator_uniform_2<FT, Point_2> Creator;
  
  FT r = 1;
  cout << "Please, enter radius of the disk:" << std::endl;
  cin >> r;
  assert(r > 0);
  
  FT eps = 0.001;
  cout << "Please, enter epsilon:" << std::endl;
  cin >> eps;
  assert(eps > 0);
  
  for( int nb = 10000, degree = 4; degree < 8; nb = nb*10, degree++) {
    
    vector<Point_2> pts;
    // disk Poincare is of radius 1 by default
    Hyperbolic_random_points_in_disc_2<K>(pts, 10*nb, 1, eps);
    
    string file_name;
    if (eps == 1e-3) {
      file_name = string("pointsEps3Nb");
    }
    
    if (eps == 1e-7) {
      file_name = string("pointsEps7Nb");
    }
    
    // file name
    stringstream convert;
    convert << degree;//add the value of Number to the characters in the stream
    file_name = file_name + convert.str();
    
    // output file
    ofstream output;
    output.open(file_name.c_str());
    
    // write radius to the file
    output << r << std::endl;
    
    // write eps to the file
    output << eps << std::endl;
    
    // write nb to the file
    output << nb << std::endl;
    
    // write random points to the file
    for(int i = 0; i < 10*nb ; i++) {
      output << pts[i] << std::endl;
    }
    
    output.close();
    
    std::cout << "file name: " << file_name << std::endl;
    std::cout << "Done." << std::endl;
  }
  
  return 0;
}
