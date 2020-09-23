#include <CGAL/Simple_cartesian.h>

#include <CGAL/Polyhedron_3.h>

#include <CGAL/Surface_mesh_parameterization/IO/File_off.h>
#include <CGAL/Surface_mesh_parameterization/Square_border_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/Discrete_conformal_map_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/parameterize.h>

#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Unique_hash_map.h>

#include <boost/array.hpp>
#include <boost/unordered_set.hpp>

#include <cstdlib>
#include <iostream>
#include <fstream>

typedef CGAL::Simple_cartesian<double>                           Kernel;
typedef Kernel::Point_2                                          Point_2;
typedef Kernel::Point_3                                          Point_3;
typedef CGAL::Polyhedron_3<Kernel>                               PolyMesh;

typedef boost::graph_traits<PolyMesh>::halfedge_descriptor       halfedge_descriptor;
typedef boost::graph_traits<PolyMesh>::vertex_descriptor         vertex_descriptor;
typedef boost::graph_traits<PolyMesh>::face_descriptor           face_descriptor;

typedef boost::graph_traits<PolyMesh>::vertex_iterator           vertex_iterator;

typedef boost::array<vertex_descriptor, 4>                       Vd_array;

typedef CGAL::Unique_hash_map<vertex_descriptor, Point_2>        UV_uhm;
typedef boost::associative_property_map<UV_uhm>                  UV_pmap;

namespace SMP = CGAL::Surface_mesh_parameterization;

bool read_vertices(const PolyMesh& mesh,
                   const char* filename,
                   Vd_array& fixed_vertices)
{
  std::string str = filename;
  if( (str.length()) < 14 || (str.substr(str.length() - 14) != ".selection.txt") ) {
    std::cerr << "Error: vertices must be given by a *.selection.txt file" << std::endl;
    return false;
  }

  std::ifstream in(filename);
  std::string line;
  if(!std::getline(in, line)) {
    std::cerr << "Error: could not read input file: " << filename << std::endl;
    return false;
  }

  // The selection file is a list of integers, so we must build a correspondence
  // between vertices and the integers.
  std::vector<vertex_descriptor> vds;
  vds.reserve(num_vertices(mesh));
  vertex_iterator vi = vertices(mesh).begin(), vi_end = vertices(mesh).end();
  CGAL_For_all(vi, vi_end) {
    vds.push_back(*vi);
  }

  // Get the first line and read the fixed vertex indices
  std::size_t counter = 0;
  std::istringstream point_line(line);
  std::size_t s;
  boost::unordered_set<std::size_t> indices;
  while(point_line >> s) {
    if(s >= vds.size())
    {
      std::cerr << "Error: Vertex index too large" << std::endl;
      return false;
    }

    vertex_descriptor vd = vds[s];
    if(!is_border(vd, mesh)) { // must be on the border
      std::cerr << "Error: vertex is not on the border of the mesh" << std::endl;
      return false;
    }

    if(counter >= 4) { // too many border vertices
      std::cerr << "Error: Too many vertices are fixed" << std::endl;
      return false;
    }

    fixed_vertices[counter++] = vd;
    indices.insert(s);
  }

  if(indices.size() < 4) {
    std::cerr << "Error: at least four unique vertices must be provided" << std::endl;
    return false;
  }

  return true;
}

int main(int argc, char** argv)
{
  std::ifstream in((argc>1) ? argv[1] : "data/nefertiti.off");
  if(!in){
    std::cerr << "Error: problem loading the input data" << std::endl;
    return 1;
  }

  PolyMesh sm;
  in >> sm;

  halfedge_descriptor bhd = CGAL::Polygon_mesh_processing::longest_border(sm).first;

  // The 2D points of the uv parametrisation will be written into this map
  UV_uhm uv_uhm;
  UV_pmap uv_map(uv_uhm);

  const char* filename = (argc > 2) ? argv[2] : "data/square_corners.selection.txt";
  Vd_array vda;
  if(!read_vertices(sm, filename, vda)) {
    std::cerr << "Error: problem loading the square corners" << std::endl;
    return EXIT_FAILURE;
  }

  typedef SMP::Square_border_uniform_parameterizer_3<PolyMesh> Border_parameterizer;
  typedef SMP::Discrete_conformal_map_parameterizer_3<PolyMesh, Border_parameterizer> Parameterizer;

  // Border parameterizers (pick one)
  Border_parameterizer border_param(vda[0], vda[1], vda[2], vda[3]);
//  Border_parameterizer border_param; // the border parameterizer will compute the corner vertices

  SMP::Error_code err = SMP::parameterize(sm, Parameterizer(border_param), bhd, uv_map);

  if(err != SMP::OK) {
    std::cerr << "Error: " << SMP::get_error_message(err) << std::endl;
    return 1;
  }

  std::ofstream out("result.off");
  SMP::IO::output_uvmap_to_off(sm, bhd, uv_map, out);

  return EXIT_SUCCESS;
}
