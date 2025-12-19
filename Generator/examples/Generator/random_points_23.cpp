#include <CGAL/Simple_cartesian.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/point_generators_3.h>
#include <string>

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_2 Point_2;
typedef K::Point_3 Point_3;
typedef CGAL::Random_points_in_square_2<Point_2> Random_points_in_square_2;
typedef CGAL::Random_points_in_cube_3<Point_3> Random_points_in_cube_3;

int main(int argc, char* argv[])
{
  int n = (argc>1)? std::stoi(argv[1]): 10;
  int dim = (argc>2)? std::stoi(argv[2]): 2;

  std::cout.precision(17);

  if(dim == 2){
    Random_points_in_square_2 rpg(1.0);

    for(int i = 0; i < n;  ++i){
      std::cout << *rpg++ << "\n";
    }
  }else if(dim == 3){
    Random_points_in_cube_3 rpg(1.0);

    for(int i = 0; i < n;  ++i){
      std::cout << *rpg++ << "\n";
    }

  }else{
    std::cout << "dimension must be 2 or 3" << std::endl;
  }
  return 0;
}

