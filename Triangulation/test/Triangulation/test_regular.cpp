#include <CGAL/Epick_d.h>
#include <CGAL/Epeck_d.h>
#include <CGAL/point_generators_d.h>
#include <CGAL/Regular_triangulation.h>
#include <CGAL/algorithm.h>

#include <tilted_grid.h>
#include <vector>
#include <string>
#include <fstream>
#include <cstdlib>
#include <algorithm>

using namespace std;

template<typename K, typename RTri>
void test(const int d, const string & type, const int N)
{
  typedef typename K::Weighted_point_d Weighted_point;
  typedef typename K::Point_d Bare_point;

  typedef typename RTri::Full_cell_handle Full_cell_handle;
  typedef typename RTri::Face Face;
  typedef typename RTri::Finite_full_cell_const_iterator Finite_full_cell_const_iterator;
  typedef typename RTri::Finite_vertex_iterator Finite_vertex_iterator;

  // Instantiate a random point generator
  CGAL::Random rng(0);
  typedef CGAL::Random_points_in_cube_d<Bare_point> Random_points_iterator;
  Random_points_iterator rand_it(d, 1.0, rng);

  RTri rt(d);
  cerr << "\nBuilding Regular triangulation of (" << type << d << ") dimension with " << N << " points";
  assert(rt.empty());

  // Generate N random points
  std::vector<Weighted_point> points;
  for (int i = 0; i < N; ++i)
    points.push_back(Weighted_point(*rand_it++, rng.get_double(0., 10.)));

  rt.insert(points.begin(),  points.end());
  cerr << "\nChecking topology and geometry...";
  assert( rt.is_valid(true) );

  cerr << "\nTraversing finite full_cells... ";
  size_t nbfs(0), nbis(0);
  Finite_full_cell_const_iterator fsit = rt.finite_full_cells_begin();
  while( fsit != rt.finite_full_cells_end() )
    ++fsit, ++nbfs;
  cerr << nbfs << " + ";
  vector<Full_cell_handle> infinite_full_cells;
  rt.tds().incident_full_cells(rt.infinite_vertex(), back_inserter(infinite_full_cells));
  nbis = infinite_full_cells.size();
  cerr << nbis << " = " << (nbis+nbfs)
  << " = " << rt.number_of_full_cells();
  cerr << "\nThe triangulation has current dimension " << rt.current_dimension();
  CGAL_assertion( rt.number_of_full_cells() == nbis+nbfs);

  cerr << "\nTraversing finite vertices... ";
  size_t nbfv(0);
  Finite_vertex_iterator fvit = rt.finite_vertices_begin();
  while( fvit != rt.finite_vertices_end() )
    ++fvit, ++nbfv;
  cerr << nbfv <<endl;

  // Count convex hull vertices:
  if( rt.maximal_dimension() > 1 )
  {
    typedef vector<Face> Faces;
    Faces edges;
    back_insert_iterator<Faces> out(edges);
    rt.tds().incident_faces(rt.infinite_vertex(), 1, out);
    cout << "\nThere are " << edges.size() << " vertices on the convex hull.";
    edges.clear();
  }
  else // rt.maximal_dimension() == 1
  {
    typedef vector<Full_cell_handle> Cells;
    Cells cells;
    back_insert_iterator<Cells> out(cells);
    rt.tds().incident_full_cells(rt.infinite_vertex(), out);
    cout << "\nThere are " << cells.size() << " vertices on the convex hull.";
    cells.clear();
  }

  // CLEAR
  rt.clear();
  assert( -1 == rt.current_dimension() );
  assert( rt.empty() );
  assert( rt.is_valid() );
  //std::cerr << ((rt.is_valid(true)) ? "VALID!" : "NOT VALID :(") << std::endl;
}

template< int D >
void go(const int N)
{
  typedef CGAL::Epick_d<CGAL::Dimension_tag<D> > KI;
  typedef CGAL::Regular_triangulation<KI> Triangulation;
  test<KI, Triangulation>(D, "inexact static", N);

  typedef CGAL::Epick_d<CGAL::Dynamic_dimension_tag> KI_dyn;
  typedef CGAL::Regular_triangulation<KI_dyn> Triangulation_dyn;
  test<KI_dyn, Triangulation_dyn>(D, "inexact dynamic", N);

  typedef CGAL::Epeck_d<CGAL::Dimension_tag<D> > KE;
  typedef CGAL::Regular_triangulation<KE> TriangulationE;
  test<KE, TriangulationE>(D, "exact static", N);

  typedef CGAL::Epeck_d<CGAL::Dynamic_dimension_tag> KE_dyn;
  typedef CGAL::Regular_triangulation<KE_dyn> TriangulationE_dyn;
  test<KE_dyn, TriangulationE_dyn>(D, "exact dynamic", N);
}

void test_inserting_points_at_the_same_position()
{
  const int DIM = 5;
  typedef CGAL::Epick_d<CGAL::Dimension_tag<DIM> > FK;
  typedef CGAL::Regular_triangulation<FK> RTri;

  typedef RTri::Vertex_handle Vertex_handle;
  typedef FK::Weighted_point_d Weighted_point;
  typedef FK::Point_d Bare_point;

  RTri rt(DIM);

  cerr << "\nTesting insertion of points at the same position";
  assert(rt.empty());

  std::vector<double> pt;
  pt.push_back(1.2);
  pt.push_back(20.3);
  pt.push_back(10.);
  pt.push_back(8.);
  pt.push_back(7.1);

  //=======================

  // First point
  Vertex_handle vh = rt.insert(Weighted_point(Bare_point(pt.begin(), pt.end()), 0.));
  assert(rt.number_of_vertices() == 1);
  assert(rt.number_of_hidden_vertices() == 0);
  assert(vh != Vertex_handle());

  // Same point
  Vertex_handle vh2 = rt.insert(Weighted_point(Bare_point(pt.begin(), pt.end()), 0.));
  assert(rt.number_of_vertices() == 1);
  assert(rt.number_of_hidden_vertices() == 0);
  assert(vh == vh2);

  // Same position, bigger weight
  vh = rt.insert(Weighted_point(Bare_point(pt.begin(), pt.end()), 1.3));
  assert(rt.number_of_vertices() == 1);
  assert(rt.number_of_hidden_vertices() == 1);
  assert(vh != Vertex_handle());

  // Same point
  vh2 = rt.insert(Weighted_point(Bare_point(pt.begin(), pt.end()), 1.3));
  assert(rt.number_of_vertices() == 1);
  assert(rt.number_of_hidden_vertices() == 1);
  assert(vh == vh2);

  // Same position, lower weight
  vh = rt.insert(Weighted_point(Bare_point(pt.begin(), pt.end()), 1.1));
  assert(rt.number_of_vertices() == 1);
  assert(rt.number_of_hidden_vertices() == 2);
  assert(vh == Vertex_handle());

  //=======================

  // New position
  pt[3] = 0.78;
  vh = rt.insert(Weighted_point(Bare_point(pt.begin(), pt.end()), 0.));
  assert(rt.number_of_vertices() == 2);
  assert(rt.number_of_hidden_vertices() == 2);
  assert(vh != Vertex_handle());

  // Same point
  vh2 = rt.insert(Weighted_point(Bare_point(pt.begin(), pt.end()), 0.));
  assert(rt.number_of_vertices() == 2);
  assert(rt.number_of_hidden_vertices() == 2);
  assert(vh == vh2);

  // Same position, bigger weight
  vh = rt.insert(Weighted_point(Bare_point(pt.begin(), pt.end()), 1.3));
  assert(rt.number_of_vertices() == 2);
  assert(rt.number_of_hidden_vertices() == 3);
  assert(vh != Vertex_handle());

  // Same point
  vh2 = rt.insert(Weighted_point(Bare_point(pt.begin(), pt.end()), 1.3));
  assert(rt.number_of_vertices() == 2);
  assert(rt.number_of_hidden_vertices() == 3);
  assert(vh == vh2);

  // Same position, lower weight
  vh = rt.insert(Weighted_point(Bare_point(pt.begin(), pt.end()), 1.1));
  assert(rt.number_of_vertices() == 2);
  assert(rt.number_of_hidden_vertices() == 4);
  assert(vh == Vertex_handle());

  //=======================

  // New position
  pt[4] = 1.78;
  vh = rt.insert(Weighted_point(Bare_point(pt.begin(), pt.end()), 0.2));
  assert(rt.number_of_vertices() == 3);
  assert(rt.number_of_hidden_vertices() == 4);
  assert(vh != Vertex_handle());

  //=======================

  // New position
  pt[1] = 1.78;
  vh = rt.insert(Weighted_point(Bare_point(pt.begin(), pt.end()), 0.8));
  assert(rt.number_of_vertices() == 4);
  assert(rt.number_of_hidden_vertices() == 4);
  assert(vh != Vertex_handle());

  //=======================

  // New position
  pt[2] = 1.78;
  vh = rt.insert(Weighted_point(Bare_point(pt.begin(), pt.end()), 0.25));
  assert(rt.number_of_vertices() == 5);
  assert(rt.number_of_hidden_vertices() == 4);
  assert(vh != Vertex_handle());

  //=======================

  // New position
  pt[0] = 12.13;
  vh = rt.insert(Weighted_point(Bare_point(pt.begin(), pt.end()), 1.25));
  assert(rt.number_of_vertices() == 6);
  assert(rt.number_of_hidden_vertices() == 4);
  assert(vh != Vertex_handle());

  // Same position, bigger weight
  vh = rt.insert(Weighted_point(Bare_point(pt.begin(), pt.end()), 1.3));
  assert(rt.number_of_vertices() == 6);
  assert(rt.number_of_hidden_vertices() == 5);
  assert(vh != Vertex_handle());

  //=======================

  // New position
  pt[0] = 9.13;
  vh = rt.insert(Weighted_point(Bare_point(pt.begin(), pt.end()), 1.25));
  assert(rt.number_of_vertices() == 7);
  assert(rt.number_of_hidden_vertices() == 5);
  assert(vh != Vertex_handle());

  //=======================

  cerr << "\nChecking topology and geometry...";
  assert(rt.is_valid(true));

  rt.clear();
  assert(-1 == rt.current_dimension());
  assert(rt.empty());
  assert(rt.is_valid());
  //std::cerr << ((rt.is_valid(true)) ? "VALID!" : "NOT VALID :(") << std::endl;
}


int main(int argc, char **argv)
{
  srand(static_cast<unsigned int>(time(nullptr)));
  int N = 10;
  if( argc > 1 )
    N = atoi(argv[1]);

  test_inserting_points_at_the_same_position();

  //go<5>(N);
  go<4>(N);
  go<3>(N);
  go<2>(N);
  go<1>(N);

  cerr << endl;
  return 0;
}
