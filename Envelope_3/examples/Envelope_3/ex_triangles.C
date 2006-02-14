#include <CGAL/Cartesian.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Timer.h>

#include <CGAL/Envelope_divide_and_conquer_3.h>
#include <CGAL/Envelope_triangles_traits_3.h>
#include <CGAL/Envelope_pm_dcel.h>

#include <CGAL/Arrangement_2.h>

#include <iostream>
#include <vector>

typedef CGAL::Gmpq                                     NT;
typedef CGAL::Cartesian<NT>                            Kernel;
typedef Kernel::Point_3                                Point_3;

typedef CGAL::Envelope_triangles_traits_3<Kernel>      Traits;
typedef Traits::Xy_monotone_surface_3                  Xy_monotone_surface_3;

typedef CGAL::Envelope_pm_dcel<Traits,
                               Xy_monotone_surface_3>  Dcel;

typedef CGAL::Arrangement_2<Traits, Dcel>              Arrangement_2;

typedef CGAL::Envelope_divide_and_conquer_3<Traits, Arrangement_2>
                                                       Envelope_alg;
                                                                
using std::cout;
using std::endl;

int main(int argc, char **argv)
{
  Traits traits;
  Envelope_alg envelope_algorithm;
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
  Arrangement_2 result;
  timer.start();
  envelope_algorithm.construct_lu_envelope(triangles.begin(),
                                           triangles.end(), result);
  timer.stop();
  cout << "construct envelope took " << timer.time() << " seconds" << endl;

  cout << "after lower envelope construction result has " << endl <<
          result.number_of_vertices() << " vertices" << endl <<
          result.number_of_halfedges() << " halfedges (" <<
          result.number_of_edges() << " edges)" << endl <<
          result.number_of_faces() << " faces" << endl;

  return 0;
}
