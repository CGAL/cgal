#include <fstream>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Set_movable_separability_2/Single_mold_translational_casting/is_pullout_direction.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef Kernel::Vector_2                                        Vector_2;
typedef Kernel::Point_2                                        Point_2;
typedef Kernel::Direction_2                               Direction_2;
typedef CGAL::Polygon_2<Kernel>                           Polygon_2;

namespace SMS = CGAL::Set_movable_separability_2;
namespace casting = SMS::Single_mold_translational_casting;

// The main program:
int main(int  argc, char* argv[])
{
  Polygon_2 polygon;

  const char* filename = (argc > 1) ? argv[1] : "polygon.dat";
  std::ifstream input_file(filename);
  if (! input_file.is_open()) {
    std::cerr << "Failed to open the " << filename << std::endl;
    return -1;
  }
  input_file >> polygon;
  input_file.close();

  // Example for is_pullout_direction_single_mold_translational_casting_2 that
  // accepts the edge
  size_t index(0);
  for (auto e_it = polygon.edges_begin(); e_it != polygon.edges_end(); ++e_it,
         ++index)
  {
    auto orientation = polygon.orientation();
    auto segment_outer_circle =
      SMS::internal::get_segment_outer_circle<Kernel>(*e_it, orientation);
    auto d = segment_outer_circle.first;
    d = d.perpendicular(CGAL::CLOCKWISE);
    auto res = casting::is_pullout_direction(polygon, e_it, d);
    std::cout << "The polygon is " << (res ? "" : "not ")
              << "castable using edge "
              << index << " in vartical translation (" << d << ")" << std::endl;

  }

  std::cout << "-----------------------------------"<< std::endl;

  // Example for is_pullout_direction_single_mold_translational_casting_2 that
  // do not accepts the edge
  {
    Vector_2 v (Point_2(0,0), Point_2(1,0));
    Direction_2 d(v);
    auto res = casting::is_pullout_direction(polygon, d);
    if (res != polygon.edges_end()) {
      std::cout << "The polygon is castable in direction d (" << d
                << ") using edge "<< *res << std::endl;

    }
    else {
      std::cout << "The polygon is not castable in direction d (" << d << ")"
                << std::endl;
    }
  }

  return 0;
}
