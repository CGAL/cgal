#include <fstream>
#include <string>

// CGAL headers

#include <CGAL/Random.h>
#include <CGAL/point_generators_2.h>

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
  
  int nb = 0;
  cout << "Please, enter number of points:" << std::endl;
  cin >> nb;
  assert(nb >= 0);
  
  CGAL::Random_points_in_disc_2<Point_2, Creator> in_disk(r);
  
  vector<Point_2> pts(nb);
  
  // file name
  string file_name("points.cin");
  
  // output file
  ofstream output;
  output.open(file_name.c_str());
  
  // write radius to the file
  output << r << std::endl;
  
  // write random points to the file
  for(int i = 0; i < nb ; i++) {
    output << *in_disk << std::endl;
    in_disk++;
  }
  
  output.close();
  
  std::cout << "file name: " << file_name << std::endl;
  std::cout << "Done." << std::endl;
  
  return 0;
}
