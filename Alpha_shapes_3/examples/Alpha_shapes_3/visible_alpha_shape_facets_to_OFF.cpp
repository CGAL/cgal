#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Alpha_shape_3.h>
#include <CGAL/Alpha_shape_cell_base_3.h>
#include <CGAL/Alpha_shape_vertex_base_3.h>
#include <CGAL/Delaunay_triangulation_3.h>

#include <fstream>
#include <vector>

#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Gt;

typedef CGAL::Alpha_shape_vertex_base_3<Gt>          Vb;
typedef CGAL::Alpha_shape_cell_base_3<Gt>            Fb;
typedef CGAL::Triangulation_data_structure_3<Vb,Fb>  Tds;
typedef CGAL::Delaunay_triangulation_3<Gt,Tds>       Triangulation_3;
typedef CGAL::Alpha_shape_3<Triangulation_3>         Alpha_shape_3;

typedef Gt::Point_3                                  Point;
typedef Alpha_shape_3::Alpha_iterator                Alpha_iterator;

int main()
{
  std::vector<Point> points;

//read input
  std::ifstream is("data/bunny_5000");
  int n;
  is >> n;
  double x, y, z;
  for (int i=0;i<n;++i)
  {
    is >> x >> y >> z;
    points.push_back( Point(x, y, z) );
  }

  std::cerr << points.size() << " points read.\n";
// compute alpha shape
  Alpha_shape_3 as(points.begin(), points.end());
  Alpha_shape_3::NT alpha_solid = as.find_alpha_solid();
  as.set_alpha(alpha_solid);

  std::cerr << "alpha_solid = " << alpha_solid << "\n";
  std::cerr << as.number_of_solid_components() << " number of solid components\n";

// collect alpha-shape facets accessible from the infinity
  // marks the cells that are in the same component as the infinite vertex by flooding
  boost::unordered_set< Alpha_shape_3::Cell_handle > marked_cells;
  std::vector< Alpha_shape_3::Cell_handle > queue;
  queue.push_back( as.infinite_cell() );

  while(!queue.empty())
  {
    Alpha_shape_3::Cell_handle back = queue.back();
    queue.pop_back();

    if ( !marked_cells.insert(back).second ) continue; //already visited

    for (int i=0; i<4; ++i)
    {
      if (as.classify(Alpha_shape_3::Facet(back, i))==Alpha_shape_3::EXTERIOR &&
          marked_cells.count(back->neighbor(i))==0)
        queue.push_back( back->neighbor(i) );
    }
  }

  // filter regular facets to restrict them to those adjacent to a marked cell
  std::vector< Alpha_shape_3::Facet > regular_facets;
  as.get_alpha_shape_facets(std::back_inserter( regular_facets ), Alpha_shape_3::REGULAR );

  std::vector<Alpha_shape_3::Facet> filtered_regular_facets;
  for(Alpha_shape_3::Facet f : regular_facets)
  {
    if ( marked_cells.count(f.first)==1 )
      filtered_regular_facets.push_back(f);
    else
    {
      f = as.mirror_facet(f);
      if ( marked_cells.count(f.first)==1 )
        filtered_regular_facets.push_back(f);
    }
  }

// dump into OFF format
  // assign an id per vertex
  boost::unordered_map< Alpha_shape_3::Vertex_handle, std::size_t> vids;
  points.clear();

  for(Alpha_shape_3::Facet f : filtered_regular_facets)
  {
    for (int i=1;i<4; ++i)
    {
      Alpha_shape_3::Vertex_handle vh = f.first->vertex((f.second+i)%4);
      if (vids.insert( std::make_pair(vh, points.size()) ).second)
        points.push_back( vh->point() );
    }
  }

  // writing
  std::ofstream output("out.off");
  output << "OFF\n " << points.size() << " " << filtered_regular_facets.size() << " 0\n";
  std::copy(points.begin(), points.end(), std::ostream_iterator<Point>(output, "\n"));
  for(const Alpha_shape_3::Facet& f : filtered_regular_facets)
  {
    output << 3;

    for (int i=0;i<3; ++i)
    {
      Alpha_shape_3::Vertex_handle vh = f.first->vertex( as.vertex_triple_index(f.second, i) );
      output << " " << vids[vh];
    }
    output << "\n";
  }

  return 0;
}
