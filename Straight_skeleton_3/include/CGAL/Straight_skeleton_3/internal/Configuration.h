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
#include <CGAL/unordered_flat_map.h>
#include <CGAL/tss.h>

#include <charconv>
#include <filesystem>
#include <fstream>
#include <memory>
#include <sstream>
#include <string>

namespace CGAL {
namespace Straight_skeletons_3 {

class Configuration;

using ConfigurationSPtr = std::shared_ptr<Configuration>;

/* <-- not documented at the moment -->
 *
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
 * Entries that are absent from the configuration file keep their default values.
 *
 * \cgalHeading{%Configuration File Entries}
 *
 * The configuration file is organized in sections, each containing key-value pairs.
 * The following entries are recognized:
 *
 * [Preprocessing]
 * - **`translate_and_scale_polyhedron`:** Currently unused. %Default is `FALSE`.
 * - **`merge_coplanar_faces`:** Boolean to enable or disable merging of (almost) coplanar faces.
 *                               %Default is `TRUE`.
 * - **`coplanarity_epsilon`:** Non-negative floating-point value for coplanarity checks, when
 *                              `merge_coplanar_faces` is enabled. Adjacent faces are merged into
 *                              a region if all their vertices are within a distance of
 *                              `coplanarity_epsilon` times the length of the diagonal of the bounding box
 *                              of the input from the supporting plane of the region.
 *                              %Default is `1e-7`.
 * - **`perturbation_epsilon`:** Positive floating-point value used to create the perturbation to ensure
 *                               general position. The value is used to nudge the coefficients of the planes.
 *                               %Default is `1e-10`.
 * - **`check_degenerate_configuration`:** Boolean to enable or disable checking if the result of the
 *                                         perturbation mechanism is indeed a valid polyhedron (and try again
 *                                         in the extremely unlikely event that the perturbation
 *                                         created a degenerate configuration).
 *                                         %Default is `FALSE`.
 * - **`seed`:** Non-negative integer used to seed the random number generator of the perturbation.
 *               %Default is `0`.
 *
 * [Algorithm]
 * - **`vertex_splitter`:** String to choose the high-degree vertex splitting strategy. Available options:
 *                          `Combi_vertex_splitter`, `Convex_vertex_splitter`.
 *                           %Default is `Combi_vertex_splitter`.
 * - **`selected_combinatorial_split`:** Non-negative index of the valid combinatorial split combination
 *                                       to use, when `Combi_vertex_splitter` is selected.
 *                                       %Default is `0`.
 * - **`convex_split_optimization`:** String (`min` or `max`) used to specify optimization strategy, when
 *                                    `Convex_vertex_splitter` is selected. The measure being optimized
 *                                    is the number of convex edges in the resulting split configuration.
 *                                    %Default is `max`.
 * - **`edge_event`:** String to choose how a degenerate edge event is handled. Available options: `convex`, `reflex`, `flip`.
 *                     %Default is `convex`.
 * - **`const_offset`:** Negative floating-point value used to insert periodic events at fixed time intervals
 *                      (time is negative during the propagation).
 *                      %Default is `0` (disabled).
 * - **`save_at_stable_times`:** Non-negative time delta used to avoid saving too close to an event.
 *                               A value of `0` disables this behavior.
 *                               %Default is `0`.
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

  // entries of the stream override existing values, other entries are left unchanged
  void parse(std::istream& input)
  {
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

  bool load(const std::filesystem::path& filename)
  {
    bool result = false;
    std::ifstream input(filename);
    if (input.is_open()) {
      parse(input);
      result = true;
      input.close();
    } else {
      CGAL_SS3_IO_TRACE_V(1, "Warning: Configuration file not found.");
    }
    return result;
  }

  void load_default_values()
  {
    // Set default values for all parameters
    properties_.clear();
    properties_["Preprocessing.translate_and_scale_polyhedron"] = "FALSE";
    properties_["Preprocessing.merge_coplanar_faces"] = "TRUE";
    properties_["Preprocessing.coplanarity_epsilon"] = "1e-7";
    properties_["Preprocessing.perturbation_epsilon"] = "1e-10";
    properties_["Preprocessing.check_degenerate_configuration"] = "FALSE";
    properties_["Preprocessing.seed"] = "0";

    properties_["Algorithm.vertex_splitter"] = "Combi_vertex_splitter";
    properties_["Algorithm.selected_combinatorial_split"] = "0";
    properties_["Algorithm.convex_split_optimization"] = "max";
    properties_["Algorithm.edge_event"] = "convex";
    properties_["Algorithm.const_offset"] = "0";
    properties_["Algorithm.save_at_stable_times"] = "0";
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
      CGAL_SS3_IO_TRACE_V(1, "Warning: key=" << mapkey << " not found.");
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

  unsigned long long get_ull(const std::string& section, const std::string& key)
  {
    CGAL_precondition(is_loaded());
    unsigned long long result = 0;
    std::string value = get_string(section, key);
    if (value.length() != 0) {
      const auto [ptr, ec] = std::from_chars(value.data(), value.data() + value.size(), result);
      if (ec != std::errc() || ptr != value.data() + value.size()) {
        CGAL_warning(false);
        return 0;
      }
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

  CGAL::unordered_flat_map<std::string, std::string> properties_;
#endif /* DOXYGEN_RUNNING */
};



} // namespace Straight_skeletons_3
} // namespace CGAL

#endif /* CGAL_STRAIGHT_SKELETON_3_INTERNAL_CONFIGURATION_H */
