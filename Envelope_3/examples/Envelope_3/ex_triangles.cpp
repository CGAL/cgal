#include <CGAL/Cartesian.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Timer.h>

#include <CGAL/envelope_3.h>
#include <CGAL/Envelope_triangles_traits_3.h>


#include <iostream>
#include <vector>

typedef CGAL::Gmpq                                     NT;
typedef CGAL::Cartesian<NT>                            Kernel;
typedef Kernel::Point_3                                Point_3;

typedef CGAL::Envelope_triangles_traits_3<Kernel>      Traits_3;
typedef Traits_3::Xy_monotone_surface_3                Xy_monotone_surface_3;
typedef CGAL::Envelope_diagram_2<Traits_3>             Envelope_diagram_2;
                                                                
using std::cout;
using std::endl;

int main(int argc, char **argv)
{
  std::vector<Xy_monotone_surface_3> triangles;
  
  int triangles_num;
  std::cin >> triangles_num;
  std::cout << "number of triangles to read: " << triangles_num << std::endl;
  Point_3 a,b,c;
  for(int i=0; i<triangles_num; ++i)
  {
    std::cin >> a >> b >> c;
    Xy_monotone_surface_3 p_triangle(a,b,c);
    triangles.push_back(p_triangle);
  }
  cout << "deal with " << triangles.size() << " triangles" << endl;

  CGAL::Timer timer;
  Envelope_diagram_2 result;
  timer.start();
  CGAL::lower_envelope_3(triangles.begin(), triangles.end(), result);
  timer.stop();
  cout << "construct envelope took " << timer.time() << " seconds" << endl;

  cout << "after lower envelope construction result has " << endl <<
          result.number_of_vertices() << " vertices" << endl <<
          result.number_of_halfedges() << " halfedges (" <<
          result.number_of_edges() << " edges)" << endl <<
          result.number_of_faces() << " faces" << endl;

  return 0;
}
