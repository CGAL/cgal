#include <CGAL/Hyperbolic_octagon_translation.h>
#include <CGAL/Periodic_4_hyperbolic_Delaunay_triangulation_2.h>
#include <CGAL/Periodic_4_hyperbolic_Delaunay_triangulation_traits_2.h>

#include <CGAL/point_generators_2.h>

typedef CGAL::Periodic_4_hyperbolic_Delaunay_triangulation_traits_2<>           Traits;
typedef Traits::FT                                                              NT;
typedef CGAL::Periodic_4_hyperbolic_Delaunay_triangulation_2<Traits>            Triangulation;

typedef Triangulation::Vertex_handle                                            Vertex_handle;
typedef Triangulation::Point                                                    Point;
typedef Triangulation::size_type                                                size_type;
typedef CGAL::Creator_uniform_2<NT,Point>                                       Creator;

int main(int argc, char** argv)
{
  int N;
  if(argc < 2) {
    std::cout << "usage: " << argv[0] << " [number_of_points_to_insert]" << std::endl;
    std::cout << "Defaulting to 100k points..." << std::endl;
    N = 100000;
  } else {
    N = atoi(argv[1]);
  }

  int N1 = N/2;
  int N2 = N - N1;

  // Generate N random points uniformly distributed with respect to the Euclidean
  // metric in the disk with radius 0.85, which contains the fundamental domain.
  // Some of the points will be outside the octagon, so they will not be inserted.
  std::vector<Point> pts;
  CGAL::Random_points_in_disc_2<Point,Creator> g(0.85);
  std::copy_n(g, N1, std::back_inserter(pts));

  // The triangulation is automatically initialized with the dummy points.
  Triangulation tr;

  // Batch-insert new points in the triangulation. Note that by default, after
  // the insertion of a set of points, an attempt is made to remove all dummy points
  // from the triangulation. Note the third boolean parameter in the call of the
  // `insert()` function which suppresses the automatic removal of dummy points.
  std::cout << "Batch-inserting " << N1 << " random points in the triangulation... "; std::cout.flush();
  size_type N_batch_inserted = tr.insert(pts.begin(), pts.end(), false);
  std::cout << "DONE! " << std::endl;

  // Insert new points in the triangulation one by one. When points are inserted
  // one by one, dummy points are not automatically removed.
  std::cout << "Single-inserting " << N2 << " random points in the triangulation... "; std::cout.flush();
  int N_single_inserted = 0;
  for(int i=0; i<N2; ++i)
  {
    Vertex_handle vh = tr.insert(*(++g));
    if(vh != Vertex_handle())
      N_single_inserted++;
  }
  std::cout << "DONE! " << std::endl;

  // Total number of inserted points.
  std::size_t N_inserted = N_batch_inserted + N_single_inserted;

  // Finally, we try to manually remove all dummy points from the triangulation.
  std::cout << "Cleaning dummy points from the triangulation... "; std::cout.flush();
  tr.try_to_remove_dummy_vertices();
  std::cout << "DONE! " << std::endl;

  // Make sure that the triangulation is valid.
  assert(tr.is_valid());

  std::size_t NV = tr.number_of_vertices();
  std::size_t NF = tr.number_of_faces();
  std::size_t NE = tr.number_of_edges();

  // This function `tr.number_of_dummy_points()` returns the number of dummy points that
  // are currently in the triangulation.
  int NDP = tr.number_of_dummy_points();

  std::cout << std::endl;
  std::cout << "-------------- STATS --------------" << std::endl;
  std::cout << "Random points generated:           " << N << std::endl;
  std::cout << "Vertices in the triangulation:     " << NV << std::endl;
  std::cout << "Dummy points in the triangulation: " << NDP << std::endl;
  std::cout << "Random points inserted:            " << N_inserted << std::endl;
  std::cout << "Random points outside/duplicates:  " << (N - N_inserted) << std::endl;
  std::cout << std::endl;
  std::cout << "---------- COMBINATORICS ----------" << std::endl;
  std::cout << "Number of vertices NV:             " << NV << std::endl;
  std::cout << "Number of faces NF:                " << NF << std::endl;
  std::cout << "Number of edges NE:                " << NE << std::endl;

  // The number of vertices in the triangulation must equal the number of points
  // inserted, plus the dummy points in the triangulation.
  assert(N_inserted + NDP == NV);
  std::cout << "Number of vertices is correct!     " << std::endl;

  // Note that the Euler relation is already verified by the function `is_valid()`.
  assert((2 + NE - NV - NF) / 2 == 2);
  std::cout << "Euler relation verified!           " << std::endl << std::endl;

  return EXIT_SUCCESS;
}
