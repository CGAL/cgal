#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/create_straight_skeleton_3.h>
#include <CGAL/Straight_skeleton_3.h>
#include <CGAL/Straight_skeleton_3/IO/OBJ.h>

#include <CGAL/Real_timer.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/property_map.h>

#include <iostream>

namespace PMP = CGAL::Polygon_mesh_processing;
namespace SS3 = CGAL::Straight_skeletons_3;

using K = CGAL::Exact_predicates_exact_constructions_kernel;
using FT = K::FT;
using Point_3 = K::Point_3;

using Mesh = CGAL::Surface_mesh<Point_3>;
using face_descriptor = boost::graph_traits<Mesh>::face_descriptor;

using Straight_skeleton_3 = CGAL::Straight_skeleton_3<K>;
using NodeSPtr = std::shared_ptr<Straight_skeleton_3::Node>;
using SpokeSPtr = std::shared_ptr<Straight_skeleton_3::Spoke>;
using SheetSPtr = std::shared_ptr<Straight_skeleton_3::Sheet>;

int main(int argc, char** argv)
{
  std::cout.precision(17);
  std::cerr.precision(17);

  const std::filesystem::path filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/oblong.off");
  const FT time_bound = (argc > 2) ? std::stod(argv[2]) : 1;

  Mesh sm;
  if(!CGAL::IO::read_polygon_mesh(filename, sm) || CGAL::is_empty(sm) || !is_valid_face_graph(sm)) {
    std::cerr << "Error: failed to read input " << filename << std::endl;
    return EXIT_FAILURE;
  }

  // start by triangulating the input as to avoid fake coplanar faces
  if (!CGAL::is_triangle_mesh(sm)) {
    std::cout << "Triangulating input mesh..." << std::endl;
    PMP::triangulate_faces(sm);
  }

  std::cout << "Input mesh: " << vertices(sm).size() << " NV " << faces(sm).size() << " NF" << std::endl;


  // every face has weight 2.0, meaning faces move at twice the default speed
  CGAL::Constant_property_map<face_descriptor, FT> face_weight_map(2);

  CGAL::Real_timer timer;
  timer.start();

  // Main call
  auto straight_skeleton = CGAL::create_straight_skeleton_3(sm, time_bound,
                                                            CGAL::parameters::face_weight_map(face_weight_map));

  timer.stop();
  std::cout << "Elapsed: " << timer.time() << std::endl;

  if (!straight_skeleton) {
    std::cerr << "Error: failed to create straight skeleton." << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "3D Straight Skeleton:\n" << std::endl;

  for (const NodeSPtr& n : straight_skeleton->nodes()) {
    std::cout << "Node " << n->id() << " @ " << n->point() << std::endl;
    std::cout << "  Incident sheets:";
    for (const std::weak_ptr<Straight_skeleton_3::Sheet>& sheet_wptr : n->sheets()) {
      if (SheetSPtr sheet = sheet_wptr.lock()) {
        std::cout << " " << sheet->id();
      }
    }
    std::cout << std::endl;
    std::cout << "  Incident spokes:";
    for (const std::weak_ptr<Straight_skeleton_3::Spoke>& spoke_wptr : n->spokes()) {
      if (SpokeSPtr spoke = spoke_wptr.lock()) {
        std::cout << " " << spoke->id();
      }
    }
    std::cout << std::endl;
  }

  for (const SpokeSPtr& sp : straight_skeleton->spokes()) {
    std::cout << "Spoke " << sp->id() << " from node " << sp->source()->id();
    if (sp->has_target()) {
      std::cout << " to node " << sp->target()->id();
    }
    std::cout << " with incident sheets:";
    for (const std::weak_ptr<Straight_skeleton_3::Sheet>& sheet_wptr : sp->sheets()) {
      if (SheetSPtr sheet = sheet_wptr.lock()) {
        std::cout << " " << sheet->id();
      }
    }
    std::cout << std::endl;
  }

  for (const SheetSPtr& sh : straight_skeleton->sheets()) {
    std::cout << "Sheet " << sh->id() << std::endl;
    std::cout << "  Incident nodes:";
    for (const NodeSPtr& node : sh->nodes()) {
      std::cout << " " << node->id();
    }
    std::cout << std::endl;
  }

  // Walking within a sheet with next()
  const SheetSPtr& sheet = *(straight_skeleton->sheets().begin());
  SpokeSPtr spoke = *sheet->spokes().begin();

  while (spoke->has_target()) {
    const NodeSPtr& node = spoke->target();
    std::cout << "Walking in sheet " << sheet->id() << " from spoke " << spoke->id() << " to node " << node->id() << std::endl;
    if (SpokeSPtr next_spoke = spoke->next(sheet, node)) {
      std::cout << "  Next spoke in sheet: " << next_spoke->id() << std::endl;
      spoke = next_spoke;
    } else {
      std::cout << "  No next spoke in sheet." << std::endl;
      break;
    }
  }

  const FT skel_time_bound = - CGAL::abs(time_bound);
  CGAL::Straight_skeletons_3::IO::write_OBJ("skeleton.obj", straight_skeleton, skel_time_bound);

  std::cout << "Done" << std::endl;
  return EXIT_SUCCESS;
}
