#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/corefinement_operations.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <fstream>
#include <sstream>

const char* cube_a =
"OFF 8 12 0\n\
0 0 0\n\
0 1 0\n\
1 1 0\n\
1 0 0\n\
1 0 1\n\
1 1 1\n\
0 1 1\n\
0 0 1\n\
3  0 1 2\n\
3  3 0 2\n\
3  4 2 5\n\
3  4 3 2\n\
3  1 6 5\n\
3  2 1 5\n\
3  0 6 1\n\
3  0 7 6\n\
3  4 5 6\n\
3  7 4 6\n\
3  3 4 7\n\
3  0 3 7\n";

const char* cube_b=
"OFF\n\
8 12 0\n\
-1 -1 -1\n\
-1 0 -1\n\
0 0 -1\n\
0 -1 -1\n\
0 -1 0\n\
0 0 0\n\
-1 0 0\n\
-1 -1 0\n\
3  0 1 2\n\
3  3 0 2\n\
3  4 2 5\n\
3  4 3 2\n\
3  1 6 5\n\
3  2 1 5\n\
3  0 6 1\n\
3  0 7 6\n\
3  4 5 6\n\
3  7 4 6\n\
3  3 4 7\n\
3  0 3 7\n";

const char* inside_a=
"OFF 4 4 0\n\
1 0.5 0.5\n\
0.5 1 0.5\n\
0.5 0 0.5\n\
0.5 0.5 1\n\
3 0 2 1\n\
3 0 1 3\n\
3 0 3 2\n\
3 1 2 3\n";

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel ;

int main()
{
  typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
  typedef CGAL::Polyhedron_corefinement<Polyhedron> Corefinement;

  std::stringstream fpa(cube_a);
  std::stringstream fpb(cube_b);
  Polyhedron pa;
  Polyhedron pb;

  fpa >> pa;
  fpb >> pb;

  assert( pa.size_of_vertices() == 8);
  assert( pb.size_of_vertices() == 8);
  
  {
  Corefinement coref;
  std::list<std::vector<Kernel::Point_3> > polylines;
  std::vector<std::pair<Polyhedron*, int> > result;
  coref( pa, pb, std::back_inserter(polylines), std::back_inserter(result), Corefinement::Intersection_tag );
  
  assert( polylines.size() == 1 );  
  assert( polylines.begin()->size() == 1 );  
  assert( result.size() == 0 );
  }

  pb.clear();
  std::stringstream fpc(inside_a);
  fpc >> pb;
  assert( pb.size_of_vertices() == 4);

  {
  Corefinement coref;
  std::list<std::vector<Kernel::Point_3> > polylines;
  std::vector<std::pair<Polyhedron*, int> > result;
  coref( pa, pb, std::back_inserter(polylines), std::back_inserter(result), Corefinement::Intersection_tag );
  
  assert( polylines.size() == 4 ); 
  assert( polylines.begin()->size() == 1 );  
  assert( result.size() == 1 );
  assert( result[0].first->size_of_vertices() == 4);
  }
}
