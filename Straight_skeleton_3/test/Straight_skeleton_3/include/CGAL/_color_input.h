// Copyright (c) 2024-2025 GeometryFactory (France)
// This file is part of CGAL (www.cgal.org)
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_STRAIGHT_SKELETON_3_TEST_COLOR_INPUT_H
#define CGAL_STRAIGHT_SKELETON_3_TEST_COLOR_INPUT_H

#include <CGAL/boost/graph/IO/PLY.h>
#include <CGAL/IO/Color.h>
#include <CGAL/Named_function_parameters.h>
#include <CGAL/type_traits.h>

#include <boost/graph/graph_traits.hpp>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <string>

namespace CGAL {
namespace Straight_skeletons_3 {
namespace utils {

// Simple helper function to draw a mesh whose faces are colored according to the weights (speeds).
template<typename PolygonMesh, typename Values>
void save_colored_mesh(const PolygonMesh& pmesh,
                       const Values& values,
                       const std::string fullpath)
{
  using Color = CGAL::IO::Color;

  using face_descriptor = typename boost::graph_traits<PolygonMesh>::face_descriptor;

  using value_type = typename CGAL::cpp20::remove_cvref<decltype(values[face_descriptor()])>::type;

  std::cout << "Saving colored mesh to " << fullpath << std::endl;

  // get a unique vector of values
  std::vector<value_type> unique_values;
  for (auto f : faces(pmesh)) {
    unique_values.push_back(values[f]);
  }

  std::sort(unique_values.begin(), unique_values.end());
  unique_values.erase(std::unique(unique_values.begin(), unique_values.end()), unique_values.end());

  srand(static_cast<unsigned int>(time(NULL)));

  std::map<std::size_t, CGAL::Color> colors;
  for (const auto& value : unique_values) {
    colors[value] = Color(static_cast<unsigned char>(rand() % 256),
                          static_cast<unsigned char>(rand() % 256),
                          static_cast<unsigned char>(rand() % 256));
    std::cout << " value " << value << " has color " << colors[value] << std::endl;
  }

  auto& nc_pmesh = const_cast<PolygonMesh&>(pmesh);
  auto face_color = nc_pmesh.template add_property_map<face_descriptor, Color>("f:color").first;

  for (auto f : faces(pmesh)) {
    std::cout << "facet " << f << " with value " << values[f]
              << " gets color " << colors[values[f]] << std::endl;
    put(face_color, f, colors[values[f]]);
  }

  std::ofstream out(fullpath);
  CGAL::IO::write_PLY(out, pmesh, CGAL::parameters::face_color_map(face_color));
}

} // namespace utils
} // namespace Straight_skeletons_3
} // namespace CGAL

#endif // CGAL_STRAIGHT_SKELETON_3_TEST_COLOR_INPUT_H
