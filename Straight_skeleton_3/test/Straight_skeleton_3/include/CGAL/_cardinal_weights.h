// Copyright (c) 2024-2025 GeometryFactory (France)
// This file is part of CGAL (www.cgal.org)
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_STRAIGHT_SKELETON_3_TEST_CARDINAL_WEIGHTS_H
#define CGAL_STRAIGHT_SKELETON_3_TEST_CARDINAL_WEIGHTS_H

#include "CGAL/_color_input.h"

#include <CGAL/Straight_skeleton_3/internal/debug.h>

#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/IO/Color.h>
#include <CGAL/number_utils.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>

#include <boost/graph/graph_traits.hpp>

#include <fstream>
#include <iostream>
#include <string>

namespace CGAL {
namespace Straight_skeletons_3 {
namespace utils {

// This is an helper function that computes the weight associated to a face
// using its normal and the provided weights in the x, y and z directions.
// Weights are written in the property map `fwm`.
//
// \tparam PolygonMesh must be a model of `FaceListGraph`
// \tparam FaceWeightMap must be a model of `ReadWritePropertyMap`
//
// \param weights_filename the name of the file containing the weights
//                         in the format:
//                         x1: <value> x2: <value> y1: <value> y2: <value>
//                         [bottom: <value> top: <value>]
// \param pmesh the polygon mesh to assign weights to
// \param fwm the face weight map to write the weights to
//
// \return `true` if all weights could be read, are valid, and were successfully assigned; `false` otherwise.
template <typename PolygonMesh, typename FaceWeightMap>
bool assign_cardinal_weights(const char* weights_filename,
                             const PolygonMesh& pmesh,
                             FaceWeightMap fwm)
{
  using face_descriptor = typename boost::graph_traits<PolygonMesh>::face_descriptor;

  using GT = typename CGAL::GetGeomTraits<PolygonMesh>::type;
  using FT = typename GT::FT;
  using Point_3 = typename GT::Point_3;
  using Vector_3 = typename GT::Vector_3;

  std::ifstream weights_in(weights_filename);

  if (!weights_filename || !weights_in) {
    CGAL_SS3_TRACE_V(1, "Warning: no input weights provided; all weights are set to '1'.");
    for (face_descriptor f : faces(pmesh))
      put(fwm, f, 1.);

    return true;
  }

  std::string x1_str, x2_str, y1_str, y2_str, bot_str, top_str;
  FT x1_val, x2_val, y1_val, y2_val, bot_val, top_val;
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
    CGAL_SS3_TRACE_V(8, "bottom & top weight info detected");
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

  CGAL_SS3_TRACE_V(8, "x1_val = " << x1_val);
  CGAL_SS3_TRACE_V(8, "x2_val = " << x2_val);
  CGAL_SS3_TRACE_V(8, "y1_val = " << y1_val);
  CGAL_SS3_TRACE_V(8, "y2_val = " << y2_val);
  CGAL_SS3_TRACE_V(8, "bot_val = " << bot_val);
  CGAL_SS3_TRACE_V(8, "top_val = " << top_val);

  if(is_zero(x1_val) && is_zero(x2_val) &&
     is_zero(y1_val) && is_zero(y2_val) &&
     is_zero(bot_val) && is_zero(top_val)) {
    std::cerr << "Error: all weights are zero" << std::endl;
    return false;
  }

  FT eps_weight = std::numeric_limits<double>::max(); // 'double' on purpose
  if(x1_val > 0) eps_weight = (std::min)(eps_weight, x1_val);
  if(x2_val > 0) eps_weight = (std::min)(eps_weight, x2_val);
  if(y1_val > 0) eps_weight = (std::min)(eps_weight, y1_val);
  if(y2_val > 0) eps_weight = (std::min)(eps_weight, y2_val);
  if(bot_val > 0) eps_weight = (std::min)(eps_weight, bot_val);
  if(top_val > 0) eps_weight = (std::min)(eps_weight, top_val);

  CGAL_SS3_TRACE_V(8, "min_weight = " << eps_weight);

  // @todo handle true zero
  eps_weight = 1e-10 * eps_weight;

  // std::cout << "HACK" << std::endl;
  // eps_weight = 0;

  if(x1_val == 0) { CGAL_SS3_TRACE_V(16, "x1_val to eps weight " << eps_weight); x1_val = eps_weight; }
  if(x2_val == 0) { CGAL_SS3_TRACE_V(16, "x2_val to eps weight " << eps_weight); x2_val = eps_weight; }
  if(y1_val == 0) { CGAL_SS3_TRACE_V(16, "y1_val to eps weight " << eps_weight); y1_val = eps_weight; }
  if(y2_val == 0) { CGAL_SS3_TRACE_V(16, "y2_val to eps weight " << eps_weight); y2_val = eps_weight; }
  if(bot_val == 0) { CGAL_SS3_TRACE_V(16, "bot_val to eps weight " << eps_weight); bot_val = eps_weight; }
  if(top_val == 0) { CGAL_SS3_TRACE_V(16, "top_val to eps weight " << eps_weight); top_val = eps_weight; }

  for (face_descriptor f : faces(pmesh))
  {
    // internal stuff, we don't need to normalize and introduce inexactness
    Vector_3 v = CGAL::NULL_VECTOR;
    CGAL::Polygon_mesh_processing::internal::sum_normals<Point_3>(
        pmesh, f, get(CGAL::vertex_point, pmesh), v, GT());
    FT sq_n = v.squared_length();

    CGAL_SS3_TRACE_V(16, "facet: " << f << " normal: " << v);

#if 1
    FT sq_cos_x = CGAL::square(v.x()) / sq_n;
    FT sq_cos_y = CGAL::square(v.y()) / sq_n;
    FT sq_cos_z = CGAL::square(v.z()) / sq_n;

    // The weight is a weighted sum of all speed contributions, where the weights are
    // the squared cosines of the angles between the normal and the axes
    FT weight_x = (v.x() >= 0) ? x1_val : x2_val;
    FT weight_y = (v.y() >= 0) ? y1_val : y2_val;
    FT weight_z = (v.z() >= 0) ? top_val : bot_val;
    FT weight = weight_x*sq_cos_x + weight_y*sq_cos_y + weight_z*sq_cos_z;
#else
    if(v.x() == 0 && v.y() == 0) {
      if(v.z() > 0)
        weight = vz2;
      else
        weight = vz1;
    } else {
      if(v.x() >= 0) {
        const FT sq_cos = CGAL::square(CGAL::scalar_product(v, west)) / sq_n;
        if(v.y() >= 0) {
          // north east quadrant
          weight = vy2 * (1 - sq_cos) + vx2 * sq_cos;
        } else {
          // south east quadrant
          weight = vy1 * (1 - sq_cos) + vx2 * sq_cos;
        }
      } else { // x < 0
        const FT sq_cos = CGAL::square(CGAL::scalar_product(v, east)) / sq_n;
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

    CGAL_SS3_TRACE_V(16, "facet: " << f << " weight: " << weight);

    // @todo currently 'double', but only because the run_and_compare.sh pipeline
    // with multiple .cpp needs to save to a file and the f:weight pmap will not
    // be detected if its value type is, e.g., EPECK::FT
    // Could be FT if there were no intermediary saving
    put(fwm, f, CGAL::to_double(weight));
    CGAL_postcondition(get(fwm, f) > 0);
  }

  CGAL_SS3_TRACE_V(8, "E-W-S-N weights: " << x1_val << " " << x2_val << " " << y1_val << " " << y2_val);

  utils::save_colored_mesh(pmesh, fwm, "results/weighted.ply");

  return true;
}

} // namespace utils
} // namespace Straight_skeletons_3
} // namespace CGAL

#endif // CGAL_STRAIGHT_SKELETON_3_TEST_CARDINAL_WEIGHTS_H
