// Copyright (c) 2024-2025 GeometryFactory (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

/**
 * file   util/Configuration.h
 * author Gernot Walzl
 * date   2012-10-30
 */

#ifndef CGAL_STRAIGHT_SKELETON_3_INTERNAL_CONFIGURATION_H
#define CGAL_STRAIGHT_SKELETON_3_INTERNAL_CONFIGURATION_H

#include <CGAL/license/Straight_skeleton_3.h>

#include <CGAL/Straight_skeleton_3/internal/debug.h>
#include <CGAL/Straight_skeleton_3/IO/StringFuncs.h>

#include <CGAL/assertions.h>
#include <CGAL/tss.h>

#include <fstream>
#include <map>
#include <memory>
#include <sstream>
#include <string>

namespace CGAL {
namespace Straight_skeletons_3 {

class Configuration;

using ConfigurationSPtr = std::shared_ptr<Configuration>;

/*!
 * \ingroup PkgStraightSkeleton3Classes
 *
 * \class CGAL::Straight_skeletons_3::Configuration
 *
 * \cgalAdvancedBegin
 * This class is for advanced users.
 * \cgalAdvancedEnd
 *
 * \brief The class allows advanced users to fine tune the 3D straight skeleton construction algorithm.
 *
 * The class `Configuration` contains a set of parameters that can be used to tweak the behavior
 * of the 3D straight skeleton construction algorithm. These parameters are either set to their
 * default values or can be loaded and overridden at runtime from an `INI` configuration file,
 * provided it follows the format detailed below.
 *
 * \cgalHeading{%Configuration File Entries}
 *
 * The configuration file is organized in sections, each containing key-value pairs.
 * The following entries are recognized:
 *
 * [Preprocessing]
 * - **`truncate_precision`:** Positive floating-point value controlling numerical truncation precision.
 *                             Coordinates of the input polyhedron are rounded to the nearest
 *                             multiple of this value. %Default is `1e-7`.
 * - **`translate_and_scale_polyhedron`:** Boolean (`TRUE`/`FALSE`) to enable translation and scaling
 *                                         of the input polyhedron to fit within a unit cube.
 *                                         %Default is `FALSE`.
 * - **`merge_coplanar_faces`:** Boolean to enable or disable merging of (almost) coplanar faces.
 *                               %Default is `TRUE`.
 * - **`coplanarity_epsilon`:** Non-negative floating-point value for coplanarity checks.
 *                              This value controls the threshold for considering two facet normals
 *                              as coplanar: if the squared Euclidean distance between their normalized
 *                              normals is less than `coplanarity_epsilon` squared, the facets are
 *                              considered coplanar and may be merged.
 *                              %Default is `1e-7`.
 * - **`perturbation_epsilon`:** Positive floating-point value used to create the pertubation to ensure
 *                               general position. The value is used to nudge the coordinates of the planes.
 *                               %Default is `1e-10`.
 * - **`check_degenerate_configuration`:** Boolean to enable or disable checking if the result of the
 *                                         perturbation mechanism is indeed a valid polyhedron (and try again
 *                                         in the extremely unlikely event that the perturbation
 *                                         created a degenerate configuration).
 *                                         %Default is `FALSE`.
 *
 * [Algorithm]
 * - **`vertex_splitter`:** String to choose the high-degree vertex splitting strategy. Available options:
 *                          `Combi_vertex_splitter`, `Convex_vertex_splitter`.
 *                           %Default is `Combi_vertex_splitter`.
 * - **`selected_combinatorial_split`:** Positive index of the valid combinatorial split combination
 *                                       to use, when `Combi_vertex_splitter` is selected.
 *                                       %Default is `0`.
 * - **`convex_split_optimization`:** String (`min` or `max`) used to specify optimization strategy, when
 *                                    `Convex_vertex_splitter` is selected. The measure being optimized
 *                                    is the number of convex edges in the resulting split configuration.
 *                                    %Default is `max`.
 * - **`edge_event`:** String to choose how a degenerate edge event is handled. Available options: `convex`, `reflex`, `split`.
 *                     %Default is `convex`.
 * - **`const_offset`** Positive floating-point value used to insert periodic events at fixed time intervals.
 *                      %Default is `0` (disabled).
 * - **`stop_after_last_save_event`:** Boolean to enable or disable stopping immediately as soon as the
 *                                     last save event has been reached.
 *                                     %Default is `TRUE`.
 *
 * \sa `CGAL::create_straight_skeleton_3()`
 * \sa `CGAL::create_straight_skeleton_and_offset_polyhedra_3()`
 */
class Configuration
{
#ifndef DOXYGEN_RUNNING
public:
  static ConfigurationSPtr get_instance()
  {
    CGAL_STATIC_THREAD_LOCAL_VARIABLE_0(ConfigurationSPtr, instance);
    if (!instance) {
      instance = std::make_shared<Configuration>();
    }
    return instance;
  }

  void parse(std::istream& input)
  {
    properties_.clear();
    std::string section;
    std::string line;
    while (std::getline(input, line)) {
      line = IO::StringFuncs::trim(line);
      if (line.empty()) {
        continue;
      }
      if (IO::StringFuncs::startsWith(line, "[") && IO::StringFuncs::endsWith(line, "]")) {
        section = line.substr(1, line.length()-2);
        continue;
      }
      if (!section.empty()) {
        std::size_t pos = line.find("=");
        if (pos != std::string::npos) {
          std::string key = IO::StringFuncs::trim(line.substr(0, pos));
          std::string value = IO::StringFuncs::trim(line.substr(pos+1, line.length()-pos-1));
          std::string mapkey = section + "." + key;
          properties_[mapkey] = value;
        }
      }
    }
    // namespace pod = boost::program_options::detail;
    // std::set<std::string> options;
    // options.insert("*");
    // for (pod::config_file_iterator i(input, options), e; i != e; ++i) {
    //   properties_[i->string_key] = i->value[0];
    // }
  }

  bool load(const std::string& filename)
  {
    CGAL_SS3_IO_TRACE(filename);
    bool result = false;
    std::ifstream input(filename.c_str());
    if (input.is_open()) {
      parse(input);
      result = true;
      input.close();
    } else {
      CGAL_SS3_IO_TRACE("Error: Configuration file not found.");
    }
    return result;
  }

  void load_default_values()
  {
    // Set default values for all parameters
    properties_.clear();
    properties_["Preprocessing.truncate_precision"] = "1e-7";
    properties_["Preprocessing.translate_and_scale_polyhedron"] = "FALSE";
    properties_["Preprocessing.merge_coplanar_faces"] = "TRUE";
    properties_["Preprocessing.coplanarity_epsilon"] = "1e-7";
    properties_["Preprocessing.perturbation_epsilon"] = "1e-10";
    properties_["Preprocessing.check_degenerate_configuration"] = "FALSE";

    properties_["Algorithm.vertex_splitter"] = "Combi_vertex_splitter";
    properties_["Algorithm.selected_combinatorial_split"] = "0";
    properties_["Algorithm.convex_split_optimization"] = "max";
    properties_["Algorithm.edge_event"] = "convex";
    properties_["Algorithm.const_offset"] = "0";
    properties_["Algorithm.stop_after_last_save_event"] = "TRUE";
  }

  bool is_loaded() const
  {
    bool result = (properties_.size() > 0);
    return result;
  }

  bool contains(const std::string& section, const std::string& key)
  {
    std::string mapkey = section + "." + key;
    bool result = (properties_.find(mapkey) != properties_.end());
    return result;
  }

  std::string get_string(const std::string& section, const std::string& key)
  {
    CGAL_precondition(is_loaded());
    std::string result;
    std::string mapkey = section + "." + key;
    if (properties_.find(mapkey) == properties_.end()) {
      // map does not contain this key
      CGAL_SS3_IO_TRACE("Error: key=" << mapkey << " not found.");
    } else {
      result = properties_[mapkey];
    }
    return result;
  }

  int get_int(const std::string& section, const std::string& key)
  {
    CGAL_precondition(is_loaded());
    int result = 0;
    std::string value = get_string(section, key);
    if (value.length() != 0) {
      result = atoi(value.c_str());
    }
    return result;
  }

  double get_double(const std::string& section, const std::string& key)
  {
    CGAL_precondition(is_loaded());
    double result = 0.0;
    std::string value = get_string(section, key);
    if (value.length() != 0) {
      result = atof(value.c_str());
    }
    return result;
  }

  template <typename FT>
  FT get_FT(const std::string& section, const std::string& key)
  {
    CGAL_precondition(is_loaded());
    FT result = 0;
    std::string value = get_string(section, key);
    if (value.length() != 0) {
      std::istringstream iss(value.c_str());
      iss >> result;
    }
    return result;
  }

  bool get_Boolean(const std::string& section, const std::string& key)
  {
    CGAL_precondition(is_loaded());
    bool result = false;
    std::string value = get_string(section, key);
    if (value.length() != 0) {
      if (value.compare("1") == 0 ||
          value.compare("t") == 0 ||
          value.compare("T") == 0 ||
          value.compare("true") == 0 ||
          value.compare("TRUE") == 0 ||
          value.compare("True") == 0) {
        result = true;
      }
    }
    return result;
  }

  std::map<std::string, std::string> properties_;
#endif /* DOXYGEN_RUNNING */
};



} // namespace Straight_skeletons_3
} // namespace CGAL

#endif /* CGAL_STRAIGHT_SKELETON_3_INTERNAL_CONFIGURATION_H */
