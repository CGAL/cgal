
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_segment_traverser_3.h>

#include <assert.h>
#include <iostream>
#include <fstream>
#include <string>
#include <boost/foreach.hpp>

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

typedef CGAL::Triangulation_segment_cell_iterator_3< DT >       Cell_traverser;

int main(int argc, char* argv[])
{
  const std::vector<Point_3> points = { { -2,  0,  0 },
                                        {  2,  0,  0 }, 
                                        {  0,  1, -1 },
                                        {  0, -1, -1 },
                                        {  0,  0,  1 },
                                        {  -10, -10, -10  },
                                        {  -10, 10, -10   },
                                        {  10, 10, -10    },
                                        {  10, -10, -10   },
                                        {  -10, -10, 10   },
                                        {  -10, 10, 10    },
                                        {  10, 10, 10     },
                                        {  10, -10, 10    },
  };
  std::vector<DT::Vertex_handle> vertices;
  vertices.reserve(points.size());
  DT dt;
  for(auto p: points) vertices.push_back(dt.insert(p));
  DT::Cell_handle c;
  assert( dt.is_valid() );
  assert(dt.is_cell(vertices[0], vertices[2], vertices[3], vertices[4], c));
  assert(dt.is_cell(vertices[1], vertices[2], vertices[3], vertices[4], c));

  std::cerr << dt.number_of_finite_cells() << '\n';
  Cell_traverser ct(dt, points[0], points[1]);
  unsigned int nb_cells = 0, nb_facets = 0, nb_edges = 0, nb_vertex = 0;

// Count the number of finite cells traversed.
      unsigned int inf = 0, fin = 0;
      for( ; ct != ct.end(); ++ct )
      {
        std::cerr << "Cell ( ";
        for(int i = 0; i < 4; ++i)
          std::cerr << ct->vertex(i)->point() << "  ";
        std::cerr << " )\n";
        

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
            ++nb_cells;
           CGAL_assertion(lt == DT::Locate_type::CELL);
          }
        }
      }

#ifdef CGAL_TRIANGULATION_3_VERBOSE_TRAVERSER_EXAMPLE
      std::cout << "While traversing from " << ct.source()
                << " to " << ct.target() << std::endl;
      std::cout << inf << " infinite and "
                << fin << " finite cells were visited." << std::endl;
      std::cout << std::endl << std::endl;
#endif

    std::cout << "Triangulation has " << dt.number_of_vertices() << " vertices."
      << std::endl;
    std::cout << "\tnb cells     : " << nb_cells << std::endl;
    std::cout << "\tnb facets    : " << nb_facets << std::endl;
    std::cout << "\tnb edges     : " << nb_edges << std::endl;
    std::cout << "\tnb vertices  : " << nb_vertex << std::endl;

     return 0;
}
