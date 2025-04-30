// Copyright (c) 2024 GeometryFactory (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

#include "algo/3d/OutwardMeshOffset.h"

#include "debug.h"
#include "util/Configuration.h"
#include "algo/Controller.h"
#include "data/3d/Polyhedron.h"
#include "algo/3d/PolyhedronTransformation.h"
#include "algo/3d/SelfIntersection.h"
#include "algo/3d/SimpleStraightSkel.h"
#include "db/3d/PLYFile.h"

#include <CGAL/Bbox_3.h>
#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <list>

namespace algo { namespace _3d {

bool
OutwardMeshOffset::
assign_weights(Mesh& sm,
               const char* weights_filename)
{
  // do not change the pmap's name without changing it in 'PLYFile.cpp'
  auto res = sm.add_property_map<face_descriptor, double>("f:weight");
  if(!res.second) {
    std::cerr << "Error: failed to add property map?" << std::endl;
    return false;
  }

  auto fwm = res.first;

  if(!weights_filename) {
    DEBUG_PRINT("No input weights provided; all weights are set to '1'.");
    for(face_descriptor f : faces(sm))
      put(fwm, f, 1.);

    return true;
  }

  std::ifstream weights_in(weights_filename) ;
  std::string x1, x2, y1, y2, z1, z2;
  CGAL::FT vx1, vx2, vy1, vy2, vz1, vz2;
  vx1 = vx2 = vy1 = vy2 = vz1 = vz2 = 0;

  if(!(weights_in >> x1 >> vx1
                  >> x2 >> vx2
                  >> y1 >> vy1
                  >> y2 >> vy2)) {
    std::cerr << "Error: failed to read weights" << std::endl;
    return false;
  }

  if(weights_in >> z1 >> vz1
                >> z2 >> vz2) {
    DEBUG_PRINT("bottom & top weight info detected");
  } else {
    if(vx2 != vy2) {
      std::cerr << "Warning: unknown z-speeds, and x-speed and y-speed differ..." << std::endl;
      // arbitrary choice
      vz2 = vy2;
    } else {
      // assign the uniform speed to the "up" direction
      vz2 = vx2;
    }
  }

  if(vx1 < 0 || vx2 < 0 || vy1 < 0 || vy2 < 0 || vz1 < 0 || vz2 < 0) {
    std::cerr << "Error: negative weights?" << std::endl;
    return false;
  }

  DEBUG_PRINT("x1_val = " << x1_val);
  DEBUG_PRINT("x2_val = " << x2_val);
  DEBUG_PRINT("y1_val = " << y1_val);
  DEBUG_PRINT("y2_val = " << y2_val);
  DEBUG_PRINT("bot_val = " << bot_val);
  DEBUG_PRINT("top_val = " << top_val);

  CGAL::FT eps_weight = std::numeric_limits<double>::max(); // 'double' on purpose
  eps_weight = (std::min)(eps_weight, x1_val);
  eps_weight = (std::min)(eps_weight, x2_val);
  eps_weight = (std::min)(eps_weight, y1_val);
  eps_weight = (std::min)(eps_weight, y2_val);
  eps_weight = (std::min)(eps_weight, bot_val);
  eps_weight = (std::min)(eps_weight, top_val);

  DEBUG_PRINT("min_weight = " << eps_weight);

  if(eps_weight == 0.) {
    std::cerr << "Error: all weights are zero" << std::endl;
    return false;
  }

  eps_weight = 1e-10 * eps_weight;

  if(vx1 == 0.) { DEBUG_PRINT("vx1 to eps weight"); vx1 = eps_weight; }
  if(vx2 == 0.) { DEBUG_PRINT("vx2 to eps weight"); vx2 = eps_weight; }
  if(vy1 == 0.) { DEBUG_PRINT("vy1 to eps weight"); vy1 = eps_weight; }
  if(vy2 == 0.) { DEBUG_PRINT("vy2 to eps weight"); vy2 = eps_weight; }
  if(vz1 == 0.) { DEBUG_PRINT("vz1 to eps weight"); vz1 = eps_weight; }
  if(vz2 == 0.) { DEBUG_PRINT("vz2 to eps weight"); vz2 = eps_weight; }

  const Vector3 east  {  1,  0,  0 }; // x2
  const Vector3 south {  0, -1,  0 }; // y1
  const Vector3 west  { -1,  0,  0 }; // x1
  const Vector3 north {  0,  1,  0 }; // y2
  const Vector3 up    {  0,  0,  1 }; // z2
  const Vector3 down  {  0,  0, -1 }; // z1

  for(face_descriptor f : faces(sm))
  {
    CGAL::FT weight = 1.;

    // internal stuff, we don't need to normalize and introduce inexactness
    Vector3 v = CGAL::NULL_VECTOR;
    CGAL::Polygon_mesh_processing::internal::sum_normals<Point3>(sm, f, get(CGAL::vertex_point, sm), v, CGAL::K());
    CGAL::FT sq_n = v.squared_length();

    if(v.x() == 0 && v.y() == 0) {
      if(v.z() > 0)
        weight = vz2;
      else
        weight = vz1;
    } else {
      if(v.x() >= 0) {
        const CGAL::FT sq_cos = CGAL::square(CGAL::scalar_product(v, west)) / sq_n;
        if(v.y() >= 0) {
          // north east quadrant
          weight = vy2 * (1 - sq_cos) + vx2 * sq_cos;
        } else {
          // south east quadrant
          weight = vy1 * (1 - sq_cos) + vx2 * sq_cos;
        }
      } else { // x < 0
        const CGAL::FT sq_cos = CGAL::square(CGAL::scalar_product(v, east)) / sq_n;
        if(v.y() >= 0) {
          // north west quadrant
          weight = vy2 * (1 - sq_cos) + vx1 * sq_cos;
        } else {
          // south west quadrant
          weight = vy1 * (1 - sq_cos) + vx1 * sq_cos;
        }
      }
    }

    DEBUG_PRINT("face: " << f << " weight: " << weight);

    // @todo currently 'double', but only because the run_and_compare.sh pipeline
    // with multiple .cpp needs to save to a file and the f:weight pmap will not
    // be detected if its value type is, e.g., EPECK::FT
    // Could be CGAL::FT if there were no intermediary saving
    put(fwm, f, CGAL::to_double(weight));
  }

  DEBUG_PRINT("E-W-S-N weights: " << vx1 << " " << vx2 << " " << vy1 << " " << vy2);

  return true;
}

bool
OutwardMeshOffset::
invert_and_add_bbox(Mesh& sm)
{
  namespace PMP = CGAL::Polygon_mesh_processing;

  DEBUG_PRINT("Inverting and adding a Bbox...");

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
  CGAL::make_hexahedron(Iso_cuboid3(Point3(bb.xmin(), bb.ymin(), bb.zmin()),
                                    Point3(bb.xmax(), bb.ymax(), bb.zmax())),
                        bbox_mesh,
                        CGAL::parameters::do_not_triangulate_faces(false));

  // if the face:weight pmap exists, get the smallest value
  // as to assign an even smaller value to the bounding box's faces
  auto fwm = sm.property_map<face_descriptor, double>("f:weight");

  double min_weight = std::numeric_limits<double>::max(); // 'double' on purpose
  if(fwm) {
    for(face_descriptor f : faces(sm))
      min_weight = (std::min)(min_weight, get(*fwm, f));
  }

  std::unordered_map<face_descriptor, face_descriptor> f2f;
  CGAL::copy_face_graph(bbox_mesh, sm,
                        CGAL::parameters::face_to_face_output_iterator(
                          std::inserter(f2f, f2f.end())));

  if(fwm) {
    for(const auto& e : f2f)
      put(*fwm, e.second, 1e-10 * min_weight);
  }

  if(CGAL::is_empty(sm) || !CGAL::is_closed(sm)) {
    std::cerr << "Error: empty or open output" << std::endl;
    return false;
  }

  return true;
}

PolyhedronSPtr
OutwardMeshOffset::
perturb(PolyhedronSPtr polyhedron)
{
    polyhedron = PolyhedronTransformation::perturb(polyhedron);
    PolyhedronTransformation::normalizeFacetPlanes(polyhedron);
    polyhedron = PolyhedronTransformation::shiftFacets(polyhedron, 0.0);
    return polyhedron;
}

bool
OutwardMeshOffset::
remove_bbox_and_invert(Mesh& sm)
{
  namespace PMP = CGAL::Polygon_mesh_processing;

  DEBUG_PRINT("Removing Bbox and inverting...");

  CGAL_precondition(!CGAL::is_empty(sm));

  auto vpm = get(CGAL::vertex_point, sm);

  vertex_descriptor extreme_v = *(vertices(sm).begin());
  for(vertex_descriptor v : vertices(sm)) {
    if(get(vpm, v).z() > get(vpm, extreme_v).z())
      extreme_v = v;
  }

  face_descriptor extreme_f = face(halfedge(extreme_v, sm), sm);
  CGAL_assertion(extreme_f != boost::graph_traits<Mesh>::null_face());

  // std::cout << vertices(sm).size() << " NV " << faces(sm).size() << " NF (with BBOX)" << std::endl;

  std::vector<face_descriptor> fs { extreme_f };
  PMP::remove_connected_components(sm, fs);

  // std::cout << vertices(sm).size() << " NV " << faces(sm).size() << " NF" << std::endl;

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

bool
OutwardMeshOffset::
run(const char* mesh_filename,
    const char* weights_filename,
    const std::list<CGAL::FT>& save_offsets,
    const std::filesystem::path save_path)
{
    util::ConfigurationSPtr config = util::Configuration::getInstance();
    std::string str_conf_file = config->findDefaultFilename();
    if (!config->load(str_conf_file)) {
        std::cerr << "Error: Config file '" << str_conf_file << "' not found." << std::endl;
    }

    Mesh sm;
    if(!CGAL::IO::read_polygon_mesh(mesh_filename, sm)) {
      std::cerr << "Error: failed to read input " << mesh_filename << std::endl;
      return false;
    }

    // convert input files into a single PLY with information
    if(!assign_weights(sm, weights_filename)) {
        std::cerr << "Error: failed to load weights " << weights_filename << std::endl;
        return false;
    }

    // the algorithm shrinks and we want an outward offset: invert the mesh
    // and add a far bounding box as to transform the input into a hole in a volume
    // whose shrinking will be equivalent to an outward offsetting
    if(!invert_and_add_bbox(sm)) {
        std::cerr << "Error: failed to add outer bounding box" << std::endl;
        return false;
    }

    // convert from CGAL::Surface_mesh to the Skeleton's Mesh Data Structure
    PolyhedronSPtr polyhedron = db::_3d::PLYFile::load(sm);

    // apply the perturbation
    polyhedron = perturb(polyhedron);

    // Failure here is likely from a bad perturbation (after all, there is a epsilon
    // probability that we create something that is degenerate).
    // @todo try again with another perturbation, or smarter (iterative) perturbation
    if (!polyhedron || !polyhedron->isConsistent()) {
        std::cerr << "Error: failed to build polyhedron (bad perturbation?)" << std::endl;
        return false;
    }

#if 0 // this could be kept, but I don't want to pay the heavy runtime cost for now
    if (SelfIntersection::hasSelfIntersectingSurface(polyhedron)) {
        std::cerr << "Error: input has self intersections (post perturbation)" << std::endl;
        return false;
    }
#endif

    // run the skeleton code
    algo::ControllerSPtr controller = { };
    SimpleStraightSkelSPtr algoskel3d =
      SimpleStraightSkel::create(polyhedron, controller, save_offsets, save_path);
    algoskel3d->run();

    for (const CGAL::FT& save_value : save_offsets) {
        DEBUG_PRINT("Post process @ " << save_value)

        // convert from the Skeleton's Mesh Data Structure to CGAL::Surface_mesh
        // @todo avoid the file I/O and either work on the SS3 data structure, or
        // add BGL <--> SS3 DS conversions
        std::stringstream tmp_ss;
        tmp_ss << save_path.string() << "/offset_" << save_value << "_triangulated.obj";
        clear(sm);
        if(!CGAL::IO::read_polygon_mesh(tmp_ss.str(), sm)) {
            std::cerr << "Error: failed to read temporary file " << tmp_ss.str() << std::endl;
            return false;
        }

        // remove the far bounding box and invert the result
        if (!remove_bbox_and_invert(sm)) {
            std::cerr << "Error: failed to remove outer bounding box" << std::endl;
            return false;
        }

        // write the final result
        std::stringstream out_ss;
        out_ss << save_path.string() << "/result_" << save_value << ".obj";

        // @todo this writes 'double', not CGAL::FT (aka EPECK::FT)
        if (!CGAL::IO::write_polygon_mesh(out_ss.str(), sm, CGAL::parameters::stream_precision(17))) {
            std::cerr << "Error: failed to write result " << out_ss.str() << std::endl;
            return false;
        }
    }

    DEBUG_PRINT("Offset meshes generated");
    return true;
}

} }
