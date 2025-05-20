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
#include "db/3d/AbstractFile.h"
#include "db/3d/OBJFile.h"
#include "db/3d/PLYFile.h"

#include <CGAL/Bbox_3.h>
#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/random_perturbation.h>
#include <CGAL/Polygon_mesh_processing/remesh_planar_patches.h>
#include <CGAL/Polygon_mesh_processing/region_growing.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <list>
#include <random>
#include <vector>

namespace algo { namespace _3d {

// @tmp
namespace utils {

template<typename PolygonMesh, typename Values>
void save_colored_mesh(PolygonMesh& pmesh,
                       const Values& values,
                       const std::string fullpath)
{
  using Color = CGAL::IO::Color;

  std::cout << "Save " << fullpath << std::endl;

  srand(static_cast<unsigned int>(time(NULL)));

  using face_descriptor = typename boost::graph_traits<PolygonMesh>::face_descriptor;

  using value_type = typename CGAL::cpp20::remove_cvref<decltype(values[face_descriptor()])>::type;

  using Face_property_color = CGAL::dynamic_face_property_t<Color>;
  using Face_color_map = typename boost::property_map<PolygonMesh, Face_property_color>::type;
  Face_color_map face_color = get(Face_property_color(), pmesh);

  // get a unique vector of values
  std::vector<value_type> unique_values;
  for (auto f : faces(pmesh)) {
    unique_values.push_back(values[f]);
  }

  std::sort(unique_values.begin(), unique_values.end());
  unique_values.erase(std::unique(unique_values.begin(), unique_values.end()), unique_values.end());

  std::map<std::size_t, CGAL::Color> colors;
  for (const auto& value : unique_values) {
    colors[value] = Color(static_cast<unsigned char>(rand() % 256),
                          static_cast<unsigned char>(rand() % 256),
                          static_cast<unsigned char>(rand() % 256));
    std::cout << " value " << value << " has color " << colors[value] << std::endl;
  }

  for (auto f : faces(pmesh)) {
    std::cout << "face " << f << " with value " << values[f] << " gets color " << colors[values[f]] << std::endl;
    put(face_color, f, colors[values[f]]);
  }

  std::ofstream out(fullpath);
  CGAL::IO::write_PLY(out, pmesh, CGAL::parameters::face_color_map(face_color));
}

} // namespace utils

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
  std::string x1_str, x2_str, y1_str, y2_str, bot_str, top_str;
  CGAL::FT x1_val, x2_val, y1_val, y2_val, bot_val, top_val;
  x1_val = x2_val = y1_val = y2_val = bot_val = top_val = 0;

  if(!(weights_in >> x1_str >> x1_val
                  >> x2_str >> x2_val
                  >> y1_str >> y1_val
                  >> y2_str >> y2_val)) {
    std::cerr << "Error: failed to read weights" << std::endl;
    return false;
  }

  CGAL_assertion(x1_str == "x1:" && x2_str == "x2:" && y1_str == "y1:" && y2_str == "y2:");

  if(weights_in >> bot_str >> bot_val
                >> top_str >> top_val) {
    DEBUG_PRINT("bottom & top weight info detected");
    CGAL_assertion(bot_str == "bottom:" && top_str == "top:");
  } else {
    if(x2_val != y2_val) {
      std::cerr << "Warning: unknown z-speeds, and x-speed and y-speed differ..." << std::endl;
      // arbitrary choice
      top_val = y2_val;
    } else {
      // assign the uniform speed to the top
      top_val = x2_val;
    }
  }

  if(x1_val < 0 || x2_val < 0 || y1_val < 0 || y2_val < 0 || bot_val < 0 || top_val < 0) {
    std::cerr << "Error: negative weights are not allowed" << std::endl;
    return false;
  }

  DEBUG_PRINT("x1_val = " << x1_val);
  DEBUG_PRINT("x2_val = " << x2_val);
  DEBUG_PRINT("y1_val = " << y1_val);
  DEBUG_PRINT("y2_val = " << y2_val);
  DEBUG_PRINT("bot_val = " << bot_val);
  DEBUG_PRINT("top_val = " << top_val);

  CGAL::FT eps_weight = std::numeric_limits<double>::max(); // 'double' on purpose
  if(x1_val > 0) eps_weight = (std::min)(eps_weight, x1_val);
  if(x2_val > 0) eps_weight = (std::min)(eps_weight, x2_val);
  if(y1_val > 0) eps_weight = (std::min)(eps_weight, y1_val);
  if(y2_val > 0) eps_weight = (std::min)(eps_weight, y2_val);
  if(bot_val > 0) eps_weight = (std::min)(eps_weight, bot_val);
  if(top_val > 0) eps_weight = (std::min)(eps_weight, top_val);

  DEBUG_PRINT("min_weight = " << eps_weight);

  if(eps_weight == std::numeric_limits<double>::max()) {
    std::cerr << "Error: all weights are zero" << std::endl;
    return false;
  }

  // @todo handle true zero
  eps_weight = 1e-10 * eps_weight;

  if(x1_val == 0) { DEBUG_PRINT("x1_val to eps weight " << eps_weight); x1_val = eps_weight; }
  if(x2_val == 0) { DEBUG_PRINT("x2_val to eps weight " << eps_weight); x2_val = eps_weight; }
  if(y1_val == 0) { DEBUG_PRINT("y1_val to eps weight " << eps_weight); y1_val = eps_weight; }
  if(y2_val == 0) { DEBUG_PRINT("y2_val to eps weight " << eps_weight); y2_val = eps_weight; }
  if(bot_val == 0) { DEBUG_PRINT("bot_val to eps weight " << eps_weight); bot_val = eps_weight; }
  if(top_val == 0) { DEBUG_PRINT("top_val to eps weight " << eps_weight); top_val = eps_weight; }

  for(face_descriptor f : faces(sm))
  {
    // internal stuff, we don't need to normalize and introduce inexactness
    Vector3 v = CGAL::NULL_VECTOR;
    CGAL::Polygon_mesh_processing::internal::sum_normals<Point3>(sm, f, get(CGAL::vertex_point, sm), v, CGAL::K());
    CGAL::FT sq_n = v.squared_length();

    DEBUG_PRINT("face: " << f << " normal: " << v);

#if 1
    CGAL::FT sq_cos_x = CGAL::square(v.x()) / sq_n;
    CGAL::FT sq_cos_y = CGAL::square(v.y()) / sq_n;
    CGAL::FT sq_cos_z = CGAL::square(v.z()) / sq_n;

    // The weight is a weighted sum of all speed contributions, where the weights are
    // the squared cosines of the angles between the normal and the axes
    CGAL::FT weight_x = (v.x() >= 0) ? x1_val : x2_val;
    CGAL::FT weight_y = (v.y() >= 0) ? y1_val : y2_val;
    CGAL::FT weight_z = (v.z() >= 0) ? top_val : bot_val;
    CGAL::FT weight = weight_x*sq_cos_x + weight_y*sq_cos_y + weight_z*sq_cos_z;
#else
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
#endif

    DEBUG_PRINT("face: " << f << " weight: " << weight);

    // @todo currently 'double', but only because the run_and_compare.sh pipeline
    // with multiple .cpp needs to save to a file and the f:weight pmap will not
    // be detected if its value type is, e.g., EPECK::FT
    // Could be CGAL::FT if there were no intermediary saving
    put(fwm, f, CGAL::to_double(weight));
    CGAL_postcondition(get(fwm, f) != 0);
  }

  DEBUG_PRINT("E-W-S-N weights: " << x1_val << " " << x2_val << " " << y1_val << " " << y2_val);

  utils::save_colored_mesh(sm, fwm, "results/weighted.ply");

  return true;
}

bool
OutwardMeshOffset::
invert_and_add_bbox(Mesh& sm)
{
  namespace PMP = CGAL::Polygon_mesh_processing;

  DEBUG_PRINT("Inverting and adding a Bbox...");

  auto fwm = sm.property_map<face_descriptor, double>("f:weight");
  CGAL_assertion(bool(fwm));

  struct Weight_setter_visitor
    : public CGAL::Polygon_mesh_processing::Triangulate_faces::Default_visitor<Mesh>
  {
    Mesh::Property_map<Mesh::Face_index, double> property;
    double weight = 1.0;
    void before_subface_creations(face_descriptor f_old) { weight = get(property, f_old); }
    void after_subface_created(face_descriptor f_new) { put(property, f_new, weight); }
  };

  Weight_setter_visitor visitor;
  visitor.property = *fwm;

  PMP::triangulate_faces(sm, CGAL::parameters::visitor(visitor));

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

#define CGAL_SS3_DO_NOT_USE_ENCLOSING_BBOX
#ifndef CGAL_SS3_DO_NOT_USE_ENCLOSING_BBOX
  const CGAL::Bbox_3 bb = PMP::bbox(sm, CGAL::parameters::bbox_scaling(100));
  Mesh bbox_mesh;
  CGAL::make_hexahedron(Iso_cuboid3(Point3(bb.xmin(), bb.ymin(), bb.zmin()),
                                    Point3(bb.xmax(), bb.ymax(), bb.zmax())),
                        bbox_mesh,
                        CGAL::parameters::do_not_triangulate_faces(false));

  // if the face:weight pmap exists, get the smallest value
  // as to assign an even smaller value to the bounding box's faces
  double min_weight = std::numeric_limits<double>::max(); // 'double' on purpose
  for(face_descriptor f : faces(sm)) {
    min_weight = (std::min)(min_weight, get(*fwm, f));
  }
  DEBUG_PRINT("min weight: " << min_weight)

  std::unordered_map<face_descriptor, face_descriptor> f2f;
  CGAL::copy_face_graph(bbox_mesh, sm,
                        CGAL::parameters::face_to_face_output_iterator(
                          std::inserter(f2f, f2f.end())));


  if(fwm) {
    for(const auto& e : f2f)
      put(*fwm, e.second, 1e-10 * min_weight);
  }
#endif

  if(CGAL::is_empty(sm) || !CGAL::is_closed(sm)) {
    std::cerr << "Error: empty or open output" << std::endl;
    return false;
  }

  return true;
}

bool
OutwardMeshOffset::
remove_bbox_and_invert(Mesh& sm)
{
  namespace PMP = CGAL::Polygon_mesh_processing;

  DEBUG_PRINT("Removing Bbox and inverting...");

  std::cout << vertices(sm).size() << " NV " << faces(sm).size() << " NF (with BBOX)" << std::endl;

  CGAL_precondition(!CGAL::is_empty(sm));
  CGAL_precondition(CGAL::is_valid_face_graph(sm) && sm.is_valid());

  if(CGAL::is_empty(sm)) {
    std::cerr << "Error: empty output" << std::endl;
    return false;
  } else if(!CGAL::is_closed(sm)) {
    std::cerr << "Error: open output" << std::endl;
    return false;
  }

  // @fixme this shouldn't be needed, but that's because files are currently written
  // using the point position as their ID, even if vertices have disappeared.
  // Thus, placeholders are required, and these are isolated vertices.
  // It's a gimmick to make it easier to debug (otherwise I need to use point positions
  // and not IDs...).
  CGAL::Polygon_mesh_processing::remove_isolated_vertices(sm);

#ifndef CGAL_SS3_DO_NOT_USE_ENCLOSING_BBOX
  auto vpm = get(CGAL::vertex_point, sm);

  vertex_descriptor extreme_v = *(vertices(sm).begin());
  for(vertex_descriptor v : vertices(sm)) {
    if(get(vpm, v).z() > get(vpm, extreme_v).z())
      extreme_v = v;
  }

  face_descriptor extreme_f = face(halfedge(extreme_v, sm), sm);
  CGAL_assertion(extreme_f != boost::graph_traits<Mesh>::null_face());


  std::vector<face_descriptor> fs { extreme_f };
  PMP::remove_connected_components(sm, fs);

  std::cout << vertices(sm).size() << " NV " << faces(sm).size() << " NF" << std::endl;

  if(CGAL::is_empty(sm)) {
    std::cerr << "Error: empty output" << std::endl;
    return false;
  }
#endif

  PMP::orient_to_bound_a_volume(sm);

  return true;
}

PolyhedronSPtr
OutwardMeshOffset::
convert(Mesh& sm)
{
    namespace PMP = CGAL::Polygon_mesh_processing;

    DEBUG_PRINT("Converting mesh...");

    bool merge_faces = false;

    util::ConfigurationSPtr config = util::Configuration::getInstance();
    std::string section("main");
    if (config->isLoaded() &&
        config->contains(section, "merge_coplanar_faces") &&
        config->getBool(section, "merge_coplanar_faces")) {
        merge_faces = true;
    }

    if (!merge_faces) {
        return db::_3d::PLYFile::load(sm);
    }

    // Use shape detection to analyze the mesh
    std::vector<std::size_t> region_ids(num_faces(sm));
    boost::vector_property_map<Plane3> plane_map; // supporting planes of the regions detected

    // detect planar regions in the mesh
    // @todo growing should:
    // - use the .ini value of 'epsilon_coplanarity'
    // - stop if it merges faces with different weights
    // - give an error for adjacent coplanar faces have different weights
    std::size_t nb_regions =
        PMP::region_growing_of_planes_on_faces(sm,
                                              CGAL::make_random_access_property_map(region_ids),
                                              CGAL::parameters::cosine_of_maximum_angle(0.98)
                                                                .region_primitive_map(plane_map)
                                                                .maximum_distance(0.001));

    utils::save_colored_mesh(sm, region_ids, "results/regions_2.ply");

    // detect corner vertices on the boundary of planar regions
    std::vector<std::size_t> corner_id_map(num_vertices(sm), -1); // corner status of vertices
    std::vector<bool> ecm(num_edges(sm), false); // mark edges at the boundary of regions

    PMP::detect_corners_of_regions(sm,
                                   CGAL::make_random_access_property_map(region_ids),
                                   nb_regions,
                                   CGAL::make_random_access_property_map(corner_id_map),
                                   CGAL::parameters::cosine_of_maximum_angle(0.98).
                                                     maximum_distance(0.001).
                                                     edge_is_constrained_map(CGAL::make_random_access_property_map(ecm)));

    std::map<edge_descriptor, EdgeWPtr> e2e;
    PolyhedronSPtr polyhedron = db::_3d::PLYFile::load(sm, {}, e2e);

    std::cout << "On load, " << polyhedron->facets().size() << " facets" << std::endl;

    // If everything is degree 3 in the region growing, merge the facets
    bool all_degree_3 = true;
    vertex_iterator vit = vertices(sm).begin(), vend = vertices(sm).end();
    std::size_t cid = 0;
    for (; vit!=vend; ++vit, ++cid) {
        std::set<std::size_t> incident_regions;
        for (face_descriptor f : CGAL::faces_around_target(halfedge(*vit, sm), sm)) {
            // could break early but it's useful to know how many regions were detected
            incident_regions.insert(region_ids[f]);
        }

        if (incident_regions.size() > 3) {
            DEBUG_PRINT("Region corner with degree " << incident_regions.size() << " at " << sm.point(*vit));
            all_degree_3 = false;
            break;
        }
    }

    if (all_degree_3) {
        // merge the facets incident to an unconstrained edge (i.e., the edge is interior to a region)
        for (edge_descriptor e: edges(sm)) {
          if (ecm[e]) {
              continue;
          }

          EdgeSPtr edge = e2e[e].lock();
          if (!edge) {
              continue;
          }
          DEBUG_PRINT("Merging facets " << edge->getFacetL()->getID() << " and " << edge->getFacetR()->getID());
          DEBUG_PRINT("SM check: " << sm.point(source(e, sm)) << " " << sm.point(target(e, sm)));
          DEBUG_PRINT("SS check: " << *(edge->getVertexSrc()->getPoint()) << " " << *(edge->getVertexDst()->getPoint()));

          db::_3d::AbstractFile::mergeFacets(edge, polyhedron);
        }
    }

    db::_3d::AbstractFile::removeVerticesDegLt3(polyhedron);

    std::cout << "Converted, " << polyhedron->facets().size() << " facets" << std::endl;
    db::_3d::OBJFile::save("results/converted.obj", polyhedron, false /*do not triangulate*/);

    return polyhedron;
}

PolyhedronSPtr
OutwardMeshOffset::
preprocess(Mesh& sm)
{
    DEBUG_PRINT("Preprocessing mesh...");

    CGAL_precondition(CGAL::is_triangle_mesh(sm));

    PolyhedronSPtr polyhedron = convert(sm);
    CGAL_assertion(polyhedron && polyhedron->isConsistent());

    polyhedron = PolyhedronTransformation::perturb(polyhedron);
    CGAL_assertion(polyhedron && polyhedron->isConsistent());

    if (!polyhedron || !polyhedron->isConsistent()) {
        // Failure here is likely from a bad perturbation (after all, there is a probabily
        // epsilon that we create something that is degenerate).
        // @todo try again with another perturbation, or smarter (iterative) perturbation
        std::cerr << "Error: invalid polyhedron (bad perturbation?)" << std::endl;
        return {};
    }

    return polyhedron;
}

bool
OutwardMeshOffset::
run(const char* mesh_filename,
    const char* weights_filename,
    const std::list<CGAL::FT>& save_offsets,
    const std::filesystem::path save_path)
{
    namespace PMP = ::CGAL::Polygon_mesh_processing;

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

    auto fwm = sm.property_map<face_descriptor, double>("f:weight");
    CGAL_USE(fwm);
    CGAL_assertion(bool(fwm));
    CGAL_assertion_code(for (face_descriptor f : faces(sm)))
    CGAL_assertion(get(*fwm, f) != 0);

    // the algorithm shrinks and we want an outward offset: invert the mesh
    // and add a far bounding box as to transform the input into a hole in a volume
    // whose shrinking will be equivalent to an outward offsetting
    if(!invert_and_add_bbox(sm)) {
        std::cerr << "Error: failed to add outer bounding box" << std::endl;
        return false;
    }

    CGAL_assertion(!PMP::does_self_intersect(sm));

    PolyhedronSPtr polyhedron = preprocess(sm);

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
