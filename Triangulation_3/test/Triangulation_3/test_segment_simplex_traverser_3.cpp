#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_segment_traverser_3.h>

#include <assert.h>
#include <iostream>
#include <fstream>
#include <string>

#include <CGAL/Random.h>

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
  const std::vector<Point_3> points = { { -2,  0,  0 },
                                        {  2,  0,  0 }, 
                                        {  0,  1,  -1 },
                                        {  0, -1,  -1 },
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
  Simplex_traverser st(dt, vertices[0]->point(), vertices[1]->point());
  
    unsigned int nb_cells = 0, nb_facets = 0, nb_edges = 0, nb_vertex = 0;
    unsigned int nb_collinear = 0;


      // Count the number of finite cells traversed.
      unsigned int inf = 0, fin = 0;
      for (; st != st.end(); ++st)
      {
        if(Cell_handle(st) != Cell_handle() && dt.is_infinite(st) )
          ++inf;
        else
        {
          ++fin;
          if (st.is_facet())       ++nb_facets;
          else if (st.is_edge())   ++nb_edges;
          else if (st.is_vertex()) ++nb_vertex;
          else if (st.is_cell())   ++nb_cells;

          if (st.is_collinear())  {
            ++nb_collinear;
            std::cout << "collinear\n";
          }

          if (st.is_facet()) {
            std::cout << "facet " << std::endl;
            DT::Facet f = *st;
            std::cout << "  ( "
                      << f.first->vertex((f.second+1)&3)->point()
                      << "  "
                      << f.first->vertex((f.second+2)&3)->point()
                      << "  "
                      << f.first->vertex((f.second+3)&3)->point()
                      << " )\n";
          }
          else if (st.is_edge()) {
            std::cout << "edge " << std::endl;
            DT::Edge e = *st;
            std::cout << "  ( "
                      << e.first->vertex(e.second)->point()
                      << "  "
                      << e.first->vertex(e.third)->point()
                      << " )\n";
          }
          else if (st.is_vertex()) {
            std::cout << "vertex " << std::endl;
            DT::Vertex_handle v = *st;
            std::cout << "  ( " << v->point() << " )\n";
          }
          else
          {
            CGAL_assertion(st.is_cell());
            std::cout << "cell \n  ( ";
            DT::Cell_handle ch = *st;
            for(int i = 0; i < 4; ++i)
              std::cout << ch->vertex(i)->point() << "  ";
            std::cout << " )\n";
          }
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

    std::cout << "\tnb cells     : " << nb_cells << std::endl;
    std::cout << "\tnb facets    : " << nb_facets << std::endl;
    std::cout << "\tnb edges     : " << nb_edges << std::endl;
    std::cout << "\tnb vertices  : " << nb_vertex << std::endl;
    std::cout << "\tnb collinear : " << nb_collinear << std::endl;

    assert(nb_cells == 2);
    assert(nb_facets == 1);
    assert(nb_edges == 0);
    assert(nb_vertex == 2);

    return 0;
}
