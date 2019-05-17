#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_segment_traverser_3.h>

#include <assert.h>
#include <iostream>
#include <fstream>
#include <string>

#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/Random.h>
#include <CGAL/Timer.h>

//#define CGAL_TRIANGULATION_3_VERBOSE_TRAVERSER_EXAMPLE

// Define the kernel.
typedef CGAL::Exact_predicates_inexact_constructions_kernel     Kernel;
typedef Kernel::Point_3                                         Point_3;

// Define the structure.
typedef CGAL::Delaunay_triangulation_3< Kernel >                DT;

typedef DT::Cell_handle                                         Cell_handle;

typedef CGAL::Triangulation_segment_simplex_iterator_3<DT>      Simplex_traverser;

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

    CGAL::default_random = CGAL::Random(0);
    CGAL::Random rng(0);
    CGAL::Timer time;
    time.start();

    unsigned int nb_cells = 0, nb_facets = 0, nb_edges = 0, nb_vertex = 0;
    unsigned int nb_collinear = 0;

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
      Simplex_traverser st(dt, p1, p2);

      // Count the number of finite cells traversed.
      unsigned int inf = 0, fin = 0;
      for (; st != st.end(); ++st)
      {
        if( dt.is_infinite(st) )
            ++inf;
        else
        {
          ++fin;
          if (st.is_facet())       ++nb_facets;
          else if (st.is_edge())   ++nb_edges;
          else if (st.is_vertex()) ++nb_vertex;
          else if (st.is_cell())   ++nb_cells;

          if (st.is_collinear())   ++nb_collinear;

#ifdef CGAL_TRIANGULATION_3_VERBOSE_TRAVERSER_EXAMPLE
          if (st.is_facet())
            std::cout << "facet " << std::endl;
          else if (st.is_edge())
            std::cout << "edge " << std::endl;
          else if (st.is_vertex())
            std::cout << "vertex " << std::endl;
          else
          {
            CGAL_assertion(st.is_cell());
            std::cout << "cell " << std::endl;
          }
#endif
        }
      }

#ifdef CGAL_TRIANGULATION_3_VERBOSE_TRAVERSER_EXAMPLE
      std::cout << "While traversing from " << st.source()
                << " to " << st.target() << std::endl;
      std::cout << "\tinfinite cells : " << inf << std::endl;
      std::cout << "\tfinite cells   : " << fin << std::endl;
      std::cout << "\tfacets   : " << nb_facets << std::endl;
      std::cout << "\tedges    : " << nb_edges << std::endl;
      std::cout << "\tvertices : " << nb_vertex << std::endl;
      std::cout << std::endl << std::endl;
#endif
    }

    time.stop();
    std::cout << "Traversing simplices of triangulation with "
      << nb_seg << " segments took " << time.time() << " seconds."
      << std::endl;
    std::cout << "\tnb cells     : " << nb_cells << std::endl;
    std::cout << "\tnb facets    : " << nb_facets << std::endl;
    std::cout << "\tnb edges     : " << nb_edges << std::endl;
    std::cout << "\tnb vertices  : " << nb_vertex << std::endl;
    std::cout << "\tnb collinear : " << nb_collinear << std::endl;

    return 0;
}