#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Bbox_3.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/generators.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>

#include <iostream>
#include <fstream>
#include <unordered_map>

namespace PMP = ::CGAL::Polygon_mesh_processing;

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using FT = K::FT;
using Point = K::Point_3;
using Vector = K::Vector_3;
using Iso_cuboid = K::Iso_cuboid_3;

using Mesh = CGAL::Surface_mesh<Point>;
using vertex_descriptor = boost::graph_traits<Mesh>::vertex_descriptor;
using face_descriptor = boost::graph_traits<Mesh>::face_descriptor;

bool invert_and_add_bbox(Mesh& sm)
{
  std::cout << "Inverting and adding a Bbox..." << std::endl;

  PMP::triangulate_faces(sm);

  // check the sanity of the input
  bool has_SI = PMP::does_self_intersect(sm);
  if(has_SI) {
    std::cerr << "Error: input has self intersections" << std::endl;
    return false;
  }

  auto vol_id_map = sm.add_property_map<face_descriptor, std::size_t>().first;
  std::size_t vccn = PMP::volume_connected_components(sm, vol_id_map,
                                                      CGAL::parameters::do_orientation_tests(false));
  std::size_t ccn = PMP::internal::number_of_connected_components(sm);
  if(vccn != ccn) {
    std::cerr << "Error: input has nested connected components" << std::endl;
    return false;
  }

  PMP::orient_to_bound_a_volume(sm);
  PMP::reverse_face_orientations(sm);

  const CGAL::Bbox_3 bb = PMP::bbox(sm, CGAL::parameters::bbox_scaling(100));
  Mesh bbox_mesh;
  CGAL::make_hexahedron(Iso_cuboid(Point(bb.xmin(), bb.ymin(), bb.zmin()),
                                   Point(bb.xmax(), bb.ymax(), bb.zmax())),
                        bbox_mesh,
                        CGAL::parameters::do_not_triangulate_faces(false));

  // if the face:weight pmap exists, get the smallest value
  // as to assign an even smaller value to the bounding box's faces
  auto fwm = sm.property_map<face_descriptor, double>("f:weight");

  double min_weight = std::numeric_limits<double>::max();
  if(fwm)
    for(face_descriptor f : faces(sm))
      min_weight = (std::min)(min_weight, get(*fwm, f));

  std::unordered_map<face_descriptor, face_descriptor> f2f;
  CGAL::copy_face_graph(bbox_mesh, sm,
                        CGAL::parameters::face_to_face_output_iterator(
                          std::inserter(f2f, f2f.end())));

  if(fwm)
    for(const auto& e : f2f)
      put(*fwm, e.second, 0.01 * min_weight);

  if(CGAL::is_empty(sm) || !CGAL::is_closed(sm)) {
    std::cerr << "Error: empty or open output" << std::endl;
    return false;
  }

  return true;
}

bool remove_bbox_and_invert(Mesh& sm)
{
  std::cout << "Removing Bbox and inverting..." << std::endl;

  CGAL_precondition(!CGAL::is_empty(sm));

  auto vpm = get(CGAL::vertex_point, sm);

  vertex_descriptor extreme_v = *(vertices(sm).begin());
  for(vertex_descriptor v : vertices(sm)) {
    if(get(vpm, v).z() > get(vpm, extreme_v).z())
      extreme_v = v;
  }

  face_descriptor extreme_f = face(halfedge(extreme_v, sm), sm);
  CGAL_assertion(extreme_f != boost::graph_traits<Mesh>::null_face());

  std::cout << vertices(sm).size() << " NV " << faces(sm).size() << " NF (with BBOX)" << std::endl;

  std::vector<face_descriptor> fs { extreme_f };
  PMP::remove_connected_components(sm, fs);

  std::cout << vertices(sm).size() << " NV " << faces(sm).size() << " NF" << std::endl;

  PMP::orient_to_bound_a_volume(sm);

  if(CGAL::is_empty(sm)) {
    std::cerr << "Error: empty output" << std::endl;
    return false;
  } else if(!CGAL::is_closed(sm)) {
    std::cerr << "Error: open output" << std::endl;
    return false;
  } else if(!CGAL::is_triangle_mesh(sm)) {
    std::cerr << "Warning: non-triangle output" << std::endl;
    return false;
  }

  return true;
}

int main(int argc, char** argv)
{
  if(argc < 2)
  {
    std::cerr << "Usage: "
              << argv[0] << "\n"
              << "\tinput_filename\n"
              << "\t[invert]\n"
              << "\t[output_filename]\n"
              << std::endl;
    std::cerr << "Input format: any format readable with CGAL::IO::read_polygon_mesh()" << std::endl;
    std::cerr << "'invert' format: 'true'/'false'" << std::endl;
    std::cerr << "Output format: PLY" << std::endl;
    return EXIT_FAILURE;
  }

  const char* input_filename = argv[1];
  const std::string do_add_str = (argc > 2) ? argv[2] : "add";
  const char* output_filename = (argc > 3) ? argv[3] : "out.ply";
  const bool do_add = (do_add_str == "add");

  std::cout << "in: " << input_filename << std::endl;
  std::cout << "add: " << do_add << std::endl;
  std::cout << "out: " << output_filename << std::endl;

  Mesh sm;

  // @tmp, while https://github.com/CGAL/cgal/pull/7874 is not integrated
  // because we don't want to lose the f:weight property map while reading/writing PLY files
  const std::string ext = CGAL::IO::internal::get_file_extension(input_filename);
  if(ext == "ply") {
    std::ifstream in(input_filename);
    if(!in || !CGAL::IO::read_PLY(in, sm)) {
      std::cerr << "Error: failed to read input (PLY)" << std::endl;
      return EXIT_FAILURE;
    }
  } else {
    CGAL_assertion(!do_add); // expecting to be here for removal only
    // PMP in case the mesh is non-manifold
    if(!PMP::IO::read_polygon_mesh(input_filename, sm, CGAL::parameters::verbose(true))) {
      std::cerr << "Error: failed to read input" << std::endl;
      return EXIT_FAILURE;
    }
  }

  std::cout << "Input mesh: " << num_vertices(sm) << " NV " << num_faces(sm) << " NF" << std::endl;

  if(CGAL::is_empty(sm) || !CGAL::is_closed(sm)) {
    std::cerr << "Error: empty or open input" << std::endl;
    return EXIT_FAILURE;
  }

  const bool res = do_add ? invert_and_add_bbox(sm) : remove_bbox_and_invert(sm);

  std::ofstream out(output_filename);
  if(!out) {
    std::cerr << "Error: failed to create output file" << std::endl;
    return EXIT_FAILURE;
  }

  if(!CGAL::IO::write_PLY(out, sm,
                          CGAL::parameters::use_binary_mode(false)
                                           .stream_precision(17)))
  {
    std::cerr << "Error: failed to write" << std::endl;
    return EXIT_FAILURE;
  }

  if(res) {
    std::cout << "Done!" << std::endl;
    return EXIT_SUCCESS;
  } else {
    std::cerr << "Error during the process" << std::endl;
    return EXIT_FAILURE;
  }
}