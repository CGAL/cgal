#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup_extension.h>


#include <CGAL/algorithm.h>
#include <CGAL/Timer.h>

#include <fstream>
#include <algorithm>
#include <cstdlib>

int main(int argc, char** argv)
{
  typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
  typedef K::Point_3 Point_3;
  typedef CGAL::Surface_mesh<Point_3> Mesh;

  std::vector<Point_3> points;
  std::vector< std::vector<std::size_t> > polygons;
  Mesh ref1;
  const char* input_filename = argc<2 ? "data/blobby-shuffled.off" : argv[1];
  std::ifstream input(input_filename);
  if ( !input){
    std::cerr << "Error: can not read input file.\n";
    return 1;
  }
  typedef  K::Point_3 Point_3;
  CGAL::File_scanner_OFF scanner(input);
  points.resize(scanner.size_of_vertices());
  polygons.resize(scanner.size_of_facets());
  //read points
  for (std::size_t i = 0; i < scanner.size_of_vertices(); ++i)
  {
    double x, y, z, w;
    scanner.scan_vertex( x, y, z, w);
    points[i] = Point_3(x, y, z, w);
    scanner.skip_to_next_vertex( i);
  }
  //read triangles
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
       else
       {
         std::cerr << "Error: input file not valid.\n";
         return 1;
       }
     }
   }
  input.close();

  if(points.size() == 0 || polygons.size()==0)
  {
    std::cerr << "Error: input file not valid.\n";
    return 1;
  }

  const char* reference_filename = argc<2 ? "data/blobby.off" : argv[2];
  input.open(reference_filename);

  if ( !input || !(input >> ref1)){
    std::cerr << "Error: can not read reference file.\n";
    return 1;
  }
  input.close();

  std::cout << "Is the soup a polygon mesh ? : "
            << CGAL::Polygon_mesh_processing::is_polygon_soup_a_polygon_mesh(polygons)
            << std::endl;

  CGAL::Polygon_mesh_processing::orient_triangle_soup_with_reference_triangle_mesh<CGAL::Sequential_tag>(ref1, points, polygons);

  std::cout << "And now ? : "
            << CGAL::Polygon_mesh_processing::is_polygon_soup_a_polygon_mesh(polygons)
            << std::endl;

  CGAL::Polygon_mesh_processing::duplicate_non_manifold_edges_in_polygon_soup(points, polygons);

  std::cout << "And now ? : "
            << CGAL::Polygon_mesh_processing::is_polygon_soup_a_polygon_mesh(polygons)
            << std::endl;

  Mesh poly;
  CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(
        points, polygons, poly);

  typedef boost::property_map<Mesh, CGAL::dynamic_face_property_t<std::size_t> >::type Fccmap;
  Fccmap fccmap = get(CGAL::dynamic_face_property_t<std::size_t>(), poly);

  std::cout<<CGAL::Polygon_mesh_processing::
         connected_components(poly, fccmap)<<" CCs before merge."<<std::endl;
  CGAL::Polygon_mesh_processing::
      merge_reversible_connected_components(poly);


  std::cout<<CGAL::Polygon_mesh_processing::
         connected_components(poly, fccmap)<<" remaining CCs."<<std::endl;
  return 0;
}
