#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>

#include <assert.h>
#include <iostream>
#include <fstream>
#include <string>

#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/Random.h>

#define CGAL_TRIANGULATION_3_VERBOSE_TRAVERSER_EXAMPLE

// Define the kernel.
typedef CGAL::Exact_predicates_inexact_constructions_kernel     Kernel;
typedef Kernel::Point_3                                         Point_3;

// Define the structure.
typedef CGAL::Delaunay_triangulation_3< Kernel >  DT;
typedef DT::Segment_simplex_iterator              Segment_simplex_iterator;

int main(int argc, char* argv[])
{
  const char* fname = (argc>1) ? argv[1] : "data/blobby.xyz";
  int nb_seg = (argc > 2) ? atoi(argv[2]) : 100;

  // Reads a .xyz point set file in points.
  // As the point is the second element of the tuple (that is with index 1)
  // we use a property map that accesses the 1st element of the tuple.

  std::vector<Point_3> points;
  std::ifstream stream(fname);
  if (!stream ||
      !CGAL::read_xyz_points(stream, std::back_inserter(points)))
  {
    std::cerr << "Error: cannot read file " << fname << std::endl;
    return EXIT_FAILURE;
  }

    // Construct the Delaunay triangulation.
    DT dt( points.begin(), points.end() );
    assert( dt.is_valid() );

    CGAL::Random rng;
    std::cout << "Random seed is " << CGAL::get_default_random().get_seed() << std::endl;

    unsigned int nb_cells = 0, nb_facets = 0, nb_edges = 0, nb_vertex = 0;

    for (int i = 0; i < nb_seg; ++i)
    {
      // Construct a traverser.
      Point_3 p1(rng.get_double(-0.48, 0.31),
                 rng.get_double(-0.22, 0.22),
                 rng.get_double(-0.19, 0.19));
      Point_3 p2(rng.get_double(-0.48, 0.31),
                 rng.get_double(-0.22, 0.22),
                 rng.get_double(-0.19, 0.19));

#ifdef CGAL_TRIANGULATION_3_VERBOSE_TRAVERSER_EXAMPLE
      std::cout << "Traverser " << (i + 1)
                << "\n\t(" << p1
                << ")\n\t(" << p2 << ")" << std::endl;
#endif

      Segment_simplex_iterator st = dt.segment_traverser_simplices_begin(p1, p2);
      Segment_simplex_iterator stend = dt.segment_traverser_simplices_end();

      for (; st != stend; ++st)
      {
        if (st->dimension() == 3)          ++nb_cells;
        else if (st->dimension() == 2)     ++nb_facets;
        else if (st->dimension() == 1)     ++nb_edges;
        else if (st->dimension() == 0)     ++nb_vertex;
      }

#ifdef CGAL_TRIANGULATION_3_VERBOSE_TRAVERSER_EXAMPLE
      std::cout << "While traversing from " << st.source()
                << " to " << st.target() << std::endl;
      std::cout << "\tfacets   : " << nb_facets << std::endl;
      std::cout << "\tedges    : " << nb_edges << std::endl;
      std::cout << "\tvertices : " << nb_vertex << std::endl;
      std::cout << std::endl << std::endl;
#endif
    }

    std::cout << "Traversing simplices of triangulation with "
      << nb_seg << " segments :" << std::endl;
    std::cout << "\tnb cells     : " << nb_cells << std::endl;
    std::cout << "\tnb facets    : " << nb_facets << std::endl;
    std::cout << "\tnb edges     : " << nb_edges << std::endl;
    std::cout << "\tnb vertices  : " << nb_vertex << std::endl;

    return 0;
}
