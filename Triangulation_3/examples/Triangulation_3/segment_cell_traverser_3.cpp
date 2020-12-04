
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>

#include <assert.h>
#include <iostream>
#include <fstream>
#include <string>

#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/Random.h>


// Define the kernel.
typedef CGAL::Exact_predicates_inexact_constructions_kernel     Kernel;
typedef Kernel::Point_3                                         Point_3;

// Define the structure.
typedef CGAL::Delaunay_triangulation_3< Kernel > DT;
typedef DT::Cell_handle                          Cell_handle;
typedef DT::Segment_cell_iterator                Segment_cell_iterator;

int main(int argc, char* argv[])
{
  const char* fname = (argc>1) ? argv[1] : "data/blobby.xyz";

  std::vector<Point_3> points;
  std::ifstream stream(fname);
  if (!stream ||
      !CGAL::read_xyz_points(stream, std::back_inserter(points)))
  {
    std::cerr << "Error: cannot read file " << fname << std::endl;
    return EXIT_FAILURE;
  }

  //bbox
  double xmin = points[0].x();
  double xmax = points[0].x();
  double ymin = points[0].y();
  double ymax = points[0].y();
  double zmin = points[0].z();
  double zmax = points[0].z();

  for(const Point_3& p : points)
  {
    xmin = (std::min)(xmin, p.x());
    ymin = (std::min)(ymin, p.y());
    zmin = (std::min)(zmin, p.z());
    xmax = (std::max)(xmax, p.x());
    ymax = (std::max)(ymax, p.y());
    zmax = (std::max)(zmax, p.z());
  }

  // Construct the Delaunay triangulation.
  DT dt( points.begin(), points.end() );
  assert( dt.is_valid() );

  CGAL::Random rng;
  std::cout << "Random seed is " << CGAL::get_default_random().get_seed() << std::endl;

  ////////////////////////////////////////////////////////////
  // Construct a traverser and use begin/end iterators
  ////////////////////////////////////////////////////////////
  Point_3 p1(rng.get_double(xmin, xmax),
             rng.get_double(ymin, ymax),
             rng.get_double(zmin, zmax));
  Point_3 p2(rng.get_double(xmin, xmax),
             rng.get_double(ymin, ymax),
             rng.get_double(zmin, zmax));

  Segment_cell_iterator ct = dt.segment_traverser_cells_begin(p1, p2);
  Segment_cell_iterator ctend = dt.segment_traverser_cells_end();

  // Count the number of finite cells traversed.
  unsigned int inf = 0, fin = 0;
  for( ; ct != ctend; ++ct )
  {
    if( dt.is_infinite(ct) )
      ++inf;
    else
      ++fin;
  }

  std::cout << "While traversing from " << p1 << " to " << p2 << std::endl;
  std::cout << inf << " infinite and "
            << fin << " finite cells were visited." << std::endl;
  std::cout << std::endl << std::endl;

  ////////////////////////////////////////////////////////////
  // Construct a traverser and use range-iterator
  ////////////////////////////////////////////////////////////

  // Count the number of finite cells traversed.
  inf = 0, fin = 0;
  for (const Cell_handle ch : dt.segment_traverser_cell_handles(p1, p2))
  {
    if (dt.is_infinite(ch))
      ++inf;
    else
      ++fin;
  }

  std::cout << "While traversing from " << p1 << " to " << p2 << std::endl;
  std::cout << inf << " infinite and "
            << fin << " finite cells were visited." << std::endl;
  std::cout << std::endl << std::endl;

  return 0;
}
