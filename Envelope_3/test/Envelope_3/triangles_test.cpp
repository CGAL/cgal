#include <CGAL/Exact_rational.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Timer.h>

#include <CGAL/Env_triangle_traits_3.h>
#include <CGAL/envelope_3.h>
#include "Envelope_triangles_test_3.h"

#include <iostream>
#include <cassert>
#include <list>
#include <set>
#include <vector>
#include <map>

// leda_rational, or Gmpq, or Quotient<MP_float>
typedef CGAL::Exact_rational                                     NT;
typedef CGAL::Cartesian<NT>                                      Kernel;

typedef Kernel::Point_3                                          Point_3;
typedef Kernel::Triangle_3                                       Triangle_3;

typedef CGAL::Env_triangle_traits_3<Kernel>                      Traits_3;

typedef Traits_3::Xy_monotone_surface_3                          Xy_monotone_surface_3;
typedef Traits_3::Surface_3                                      Surface_3;

typedef CGAL::Envelope_diagram_2<Traits_3>                       Envelope_diagram_2;
typedef CGAL::Envelope_triangles_test_3<Traits_3, 
                                        Envelope_diagram_2>      Envelope_test;

std::ostream & operator << (std::ostream  & os, Triangle_3* p_data)
{
  return os << *p_data;
}

using std::cout;
using std::endl;

bool test_one_file(std::ifstream& inp)
{
  Traits_3 traits;
  std::vector<Surface_3> triangles;
  
  int triangles_num = 0;
  inp >> triangles_num;
  std::cout << "number of triangles to read: " << triangles_num << std::endl;
  Point_3 a,b,c;
  for(int i=0; i<triangles_num; ++i)
  {
    inp >> a >> b >> c;
    Surface_3 p_triangle(a,b,c);
    triangles.push_back(p_triangle);
  }
  cout << "deal with " << triangles.size() << " triangles" << endl;

  // make xy-monotone surfaces
  std::vector<Surface_3>::iterator begin;
  std::vector<Xy_monotone_surface_3> xy_monotones;
  
  for(begin = triangles.begin(); begin != triangles.end(); ++begin)
    traits.make_xy_monotone_3_object()(*begin, true, std::back_inserter(xy_monotones));

  CGAL::Timer timer;
  Envelope_diagram_2 result, test_result;
  timer.start();
  CGAL::lower_envelope_xy_monotone_3(xy_monotones.begin(), xy_monotones.end(), result);
  timer.stop();
  cout << "construct envelope took " << timer.time() << " seconds" << endl;
  timer.reset();

  cout << "after lower envelope construction result has " << endl <<
          result.number_of_vertices() << " vertices" << endl <<
          result.number_of_halfedges() << " halfedges (" << result.number_of_edges() << " edges)" << endl <<
          result.number_of_faces() << " faces" << endl;

  Envelope_test env_test;
  timer.start();
  env_test.construct_lu_envelope(xy_monotones.begin(), xy_monotones.end(), test_result);
  timer.stop();
  cout << "construct test map took " << timer.time() << " seconds" << endl;
  timer.reset();

  cout << "finished computing test map. it has " << endl <<
          test_result.number_of_vertices() << " vertices" << endl <<
          test_result.number_of_halfedges() << " halfedges (" <<
          test_result.number_of_edges() << " edges)" << endl <<
          test_result.number_of_faces() << " faces" << endl;

  timer.start();
  bool test2 = env_test.compare_lu_envelopes_test2(test_result, result, /*only faces = */false);
  timer.stop();

  return test2;
}


int main(int argc, char **argv)
{
  if(argc < 2)
  {
    std::cerr<<"Missing input file\n";
    std::exit (-1);
  }

  int success = 0;
  for(int i=1; i<argc; ++i)
  {
    std::string str(argv[i]);
    if(str.empty())
      continue;

    std::ifstream inp(argv[i]);
    if(!inp.is_open())
    {
      std::cerr<<"Failed to open " <<str<<"\n";
      return (-1);
    }
    if (! test_one_file(inp))
    {
      inp.close();
      std::cout<<str<<": ERROR\n";
      success = -1;
    }
    else
    {
      std::cout<<str<<": succeeded\n";
    }
    inp.close();
  }
  
  return success;
}
