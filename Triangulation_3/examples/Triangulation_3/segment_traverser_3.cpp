#define CGAL_TRIANGULATION_3_TRAVERSER_CHECK_INTERSECTION

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_segment_traverser_3.h>

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
typedef CGAL::Delaunay_triangulation_3< Kernel >                DT;

typedef DT::Cell_handle                                         Cell_handle;

typedef CGAL::Triangulation_segment_cell_iterator_3< DT >       Cell_traverser;
typedef CGAL::Triangulation_segment_simplex_iterator_3<DT>      Simplex_traverser;

int main(int argc, char* argv[])
{
  const char* fname = (argc>1) ? argv[1] : "data/blobby.xyz";

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

  //bbox
  //min (-0.481293,-0.220929,-0.194076), max (0.311532,0.225525,0.198025)

    // Construct the Delaunay triangulation.
    DT dt( points.begin(), points.end() );
    assert( dt.is_valid() );

    CGAL::default_random = CGAL::Random(0);
    CGAL::Random rng(0);

    unsigned int nb_facets = 0, nb_edges = 0, nb_vertex = 0;

    unsigned int nb_seg = 100;
    for (unsigned int i = 0; i < nb_seg; ++i)
    {
      // Construct a traverser.
      Point_3 p1(rng.get_double(-0.48, 0.31),
                 rng.get_double(-0.22, 0.22),
                 rng.get_double(-0.19, 0.19));
      Point_3 p2(rng.get_double(-0.48, 0.31),
                 rng.get_double(-0.22, 0.22),
                 rng.get_double(-0.19, 0.19));

      std::cout << "Traverser " << (i + 1)
        << "\n\t(" << p1
        << ")\n\t(" << p2 << ")" << std::endl;
      Cell_traverser ct(dt, p1, p2);

      // Count the number of finite cells traversed.
      unsigned int inf = 0, fin = 0;
      for( ; ct != ct.end(); ++ct )
      {
        if( dt.is_infinite(ct) )
            ++inf;
        else
        {
          ++fin;

          DT::Locate_type lt;
          int li, lj;
          ct.entry(lt, li, lj);

          switch (lt)
          {
          case DT::Locate_type::FACET:
            ++nb_facets;
            break;
          case DT::Locate_type::EDGE:
            ++nb_edges;
            break;
          case DT::Locate_type::VERTEX:
            ++nb_vertex;
            break;
          default:
            /*when source is in a cell*/
            CGAL_assertion(lt == DT::Locate_type::CELL);
          }
        }
      }

      Simplex_traverser st(dt, p1, p2);
      for (; st != st.end(); ++st)
      {
        if (st.is_vertex())
          std::cout << "vertex " << st.get_vertex()->point() << std::endl;
      }

      std::cout << "While traversing from " << ct.source()
                << " to " << ct.target() << std::endl;
      std::cout << inf << " infinite and "
                << fin << " finite cells were visited." << std::endl;
      std::cout << std::endl << std::endl;
    }

     return 0;
}