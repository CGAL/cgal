
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>

#include <assert.h>
#include <iostream>
#include <fstream>


// Define the kernel.
typedef CGAL::Exact_predicates_inexact_constructions_kernel     Kernel;
typedef Kernel::Point_3                                         Point_3;

// Define the structure.
typedef CGAL::Delaunay_triangulation_3< Kernel >  DT;
typedef DT::Cell_handle                           Cell_handle;
typedef DT::Segment_cell_iterator                 Segment_cell_iterator;

int main()
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

  Cell_handle c;
  assert( dt.is_valid() );
  assert( dt.is_cell(vertices[0], vertices[2], vertices[3], vertices[4], c));
  assert( dt.is_cell(vertices[1], vertices[2], vertices[3], vertices[4], c));

  std::cerr << dt.number_of_finite_cells() << '\n';
  Segment_cell_iterator ct = dt.segment_traverser_cells_begin(points[0], points[1]);
  Segment_cell_iterator ctend = dt.segment_traverser_cells_end();

  // Count the number of finite cells traversed.
  unsigned int inf = 0, fin = 0;
  while(ct != ctend)
  {
    std::cerr << "Cell ( ";
    for(int i = 0; i < 4; ++i)
      std::cerr << ct->vertex(i)->point() << "  ";
    std::cerr << " )\n";

    if( dt.is_infinite(ct) )
      ++inf;
    else
      ++fin;

    ++ct;
  }

  std::cout << "While traversing from " << points[0]
            << " to " << points[1] << std::endl;
  std::cout << inf << " infinite and "
            << fin << " finite cells were visited." << std::endl;

  inf = 0;
  fin = 0;
  for (Cell_handle ch : dt.segment_traverser_cell_handles(vertices[2], vertices[3]))
  {
    std::cerr << "Cell ( ";
    for (int i = 0; i < 4; ++i)
      std::cerr << ch->vertex(i)->point() << "  ";
    std::cerr << " )\n";

    if (dt.is_infinite(ch))
      ++inf;
    else
      ++fin;
  }

  std::cout << "While traversing from " << vertices[2]->point()
            << " to " << vertices[3]->point() << std::endl;
  std::cout << inf << " infinite and "
            << fin << " finite cells were visited." << std::endl;

  return 0;
}
