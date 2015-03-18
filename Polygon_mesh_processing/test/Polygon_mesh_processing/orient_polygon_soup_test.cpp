#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>

#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Timer.h>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point_3;
typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef CGAL::Surface_mesh<Point_3> Surface_mesh;

std::istream& read_soup(
  std::istream& stream, 
  std::vector<Point_3>& points,
  std::vector< std::vector<std::size_t> >& polygons) 
{
  CGAL::File_scanner_OFF scanner(stream);
  points.resize(scanner.size_of_vertices());
  polygons.resize(scanner.size_of_facets());
  for (std::size_t i = 0; i < scanner.size_of_vertices(); ++i) 
  {
    double x, y, z, w;
    scanner.scan_vertex( x, y, z, w);
    points[i] = Point_3(x, y, z, w);
    scanner.skip_to_next_vertex( i);
  }
  if(!stream) { return stream; }

  for (std::size_t i = 0; i < scanner.size_of_facets(); ++i) 
  {
    std::size_t no;
    scanner.scan_facet( no, i);
    polygons[i].resize(no);

    for(std::size_t j = 0; j < no; ++j) {
      std::size_t id;
      scanner.scan_facet_vertex_index(id, i);
      if(id < scanner.size_of_vertices())
      {
        polygons[i][j] = id;
      }
      else { return stream; }
    }
  }
  return stream;
}

int main(int,char** ) {
  std::vector<Point_3> points;
  std::vector< std::vector<std::size_t> > polygons;

  std::ifstream input("data/elephant-shuffled.off");
  if ( !input || !read_soup(input, points, polygons)){
    std::cerr << "Error: can not read file.\n";
    return 1;
  }

  bool oriented
    = CGAL::Polygon_mesh_processing::orient_polygon_soup(points, polygons);
  std::cerr << (oriented ? "Oriented." : "Not orientabled.") << std::endl;
  
  if(oriented) {
    Surface_mesh mesh;
    CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(
      points, polygons, mesh);

    Polyhedron poly;
    CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(
      points, polygons, poly);

    std::ofstream out("oriented.off");
    out << poly;
    out.close();
  }
}