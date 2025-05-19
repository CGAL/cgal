#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Scale_space_surface_reconstruction_3.h>

#include <CGAL/IO/read_points.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>

#include <fstream>
#include <iostream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel     Kernel;

typedef CGAL::Scale_space_surface_reconstruction_3< Kernel >                Reconstruction;
typedef CGAL::Scale_space_reconstruction_3::Weighted_PCA_smoother< Kernel > Smoother;
typedef CGAL::Scale_space_reconstruction_3::Alpha_shape_mesher< Kernel >    Mesher;

typedef Reconstruction::Point                                   Point;
typedef Reconstruction::Facet_const_iterator                    Facet_iterator;

typedef CGAL::Surface_mesh<Point> Surface_mesh;
typedef boost::graph_traits<Surface_mesh>::vertex_descriptor vertex_descriptor;

int main(int argc, char* argv[]) {
  // Read the data
  std::string fname = argc==1?CGAL::data_file_path("points_3/kitten.off"):argv[1];

  std::cerr << "Reading " << std::flush;
  std::vector<Point> points;
  if(!CGAL::IO::read_points(fname, std::back_inserter(points)))
  {
    std::cerr << "Error: cannot read file" << std::endl;
    return EXIT_FAILURE;
  }

  std::cerr << "done: " << points.size() << " points." << std::endl;

  // Construct the mesh in a scale space.
  Reconstruction reconstruct(points.begin(), points.end() );
  Smoother smoother(10, 200 );
  reconstruct.increase_scale(4, smoother);

  Mesher mesher(smoother.squared_radius(),
                false, // Do not separate shells
                true // Force manifold output
               );
  reconstruct.reconstruct_surface(mesher);

  Reconstruction::Point_range smoothed(reconstruct.points());
  Reconstruction::Facet_range polygons(reconstruct.facets());

  CGAL::Polygon_mesh_processing::orient_polygon_soup(smoothed, polygons);

  Surface_mesh mesh;
  CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(smoothed, polygons, mesh);

  // Also store the input points as vertex property
  Surface_mesh::Property_map<vertex_descriptor, Point> original;
  bool created;
  std::tie(original, created) = mesh.add_property_map<vertex_descriptor,Point>("v:original");
  assert(created);

  int i = 0;
  for(auto v : vertices(mesh)){
    put(original, v, points[i++]);
  }


  std::cerr << "Done." << std::endl;

  return EXIT_SUCCESS;
}
