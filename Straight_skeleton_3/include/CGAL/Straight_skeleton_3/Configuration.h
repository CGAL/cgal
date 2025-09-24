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

#include <CGAL/Straight_skeleton_3/internal/debug.h>
#include <CGAL/Straight_skeleton_3/IO/StringFuncs.h>

#include <fstream>
#include <map>
#include <memory>
#include <sstream>
#include <string>

namespace CGAL {
namespace Straight_skeletons_3 {

class Configuration;

using ConfigurationSPtr = std::shared_ptr<Configuration>;

/**
 * \ingroup PkgStraightSkeleton3Ref
 *
 * \class CGAL::Straight_skeletons_3::Configuration
 *
 * \brief Handles the reading and storage of configuration options for the 3D straight skeleton package.
 *
 * These parameters are fine-tuning options for the 3D straight skeleton algorithm for advanced users.
 *
 * This class loads and stores configuration parameters from an INI configuration file. The default
 * provided configuration file is called `StraightSkel.ini` and is located in the `examples` and `test`
 * directories. These configuration files can be directly modified by users, or another file
 * can be loaded explicitly by passing its path using the named parameter `config_file_path`,
 * provided it follows the format detailed below.
 *
 * \cgalHeading{%Configuration File Entries}
 *
 * The configuration file is organized in sections, each containing key-value pairs.
 * The following entries are recognized:
 *
 * [Preprocessing]
 * - `translate_and_scale_polyhedron`: Boolean (`TRUE`/`FALSE`) to enable translation and scaling
 *                                     of the input polyhedron to fit within a unit cube.
 * - `truncate_precision`: Positive floating-point value controlling numerical truncation precision.
 *                         Coordinates of the input polyhedron are rounded to the nearest
 *                         multiple of this value.
 * - `merge_coplanar_faces`: Boolean to enable or disable merging of (almost) coplanar faces.
 * - `coplanarity_epsilon`: Non-negative floating-point value for coplanarity checks.
 *                          This value controls the threshold for considering two facet normals
 *                          as coplanar: if the squared Euclidean distance between their normalized
 *                          normals is less than `coplanarity_epsilon` squared, the facets are
 *                          considered coplanar and may be merged.
 * - `perturbation_epsilon`: Positive floating-point value used to create the pertubation to ensure
 *                           general position. The value is used to nudge the coordinates of the planes.
 * - `check_degenerate_configuration`: Boolean to enable or disable checking if the result of the *                                     perturbation mechanism is indeed a valid polyhedron (and try again
 *                                     in the extremely unlikely event that the perturbation
 *                                     created a degenerate configuration).
 *
 * [Algorithm]
 * - `vertex_splitter`: String to choose the high-degree vertex splitting strategy. Available options:
 *                      `CombiVertexSplitter`, `ConvexVertexSplitter`.
 * - `selected_combinatorial_split`: Positive index of the valid combinatorial split combination
 *                                   to use, when `CombiVertexSplitter` is selected.
 * - `convex_split_optimization`: String (`min` or `max`) used to specify optimization strategy, when
 *                                `ConvexVertexSplitter` is selected. The measure being optimized
 *                                is the number of convex edges in the resulting split configuration.
 * - `stop_after_last_save_event`: Boolean to enable or disable stopping immediately as soon as the
 *                                 last save event has been reached.
 *
 * \sa `CGAL::Straight_skeletons_3::face_offset()`
 */
class Configuration
{
#ifndef DOXYGEN_RUNNING
public:
  static ConfigurationSPtr getInstance()
  {
    if (!instance_) {
      instance_ = ConfigurationSPtr(new Configuration());
    }
    return instance_;
  }

  // Seek a file called 'StraightSkel.ini', either in the working directory
  std::string findDefaultFilename()
  {
    std::string name("StraightSkel");
    std::string result = name + ".ini";
    std::string home(getenv("HOME"));
    std::string sysconfdir("/etc");
    std::string filenames[2];
    filenames[0] = home+"/."+name+"/"+name+".ini";
    filenames[1] = sysconfdir+"/"+name+"/"+name+".ini";
    for (unsigned int i = 0; i < 2; ++i) {
      std::ifstream input(filenames[i].c_str());
      if (input.is_open()) {
        input.close();
        result = filenames[i];
        break;
      }
    }
    return result;
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
      CGAL_SS3_IO_TRACE("Error: Config file not found.");
    }
    return result;
  }

  bool isLoaded() const
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

  std::string getString(const std::string& section, const std::string& key)
  {
    std::string result;
    std::string mapkey = section + "." + key;
    if (properties_.find(mapkey) == properties_.end()) {
      // map does not contain this key
      CGAL_SS3_IO_TRACE("Warning: key=" << mapkey << " not found.");
    } else {
      result = properties_[mapkey];
    }
    return result;
  }

  int getInt(const std::string& section, const std::string& key)
  {
    int result = 0;
    std::string value = getString(section, key);
    if (value.length() != 0) {
      result = atoi(value.c_str());
    }
    return result;
  }

  double getDouble(const std::string& section, const std::string& key)
  {
    double result = 0.0;
    std::string value = getString(section, key);
    if (value.length() != 0) {
      result = atof(value.c_str());
    }
    return result;
  }

  template <typename FT>
  FT getFT(const std::string& section, const std::string& key)
  {
    FT result = 0.0;
    std::string value = getString(section, key);
    if (value.length() != 0) {
      std::istringstream iss(value.c_str());
      iss >> result;
    }
    return result;
  }

  bool getBool(const std::string& section, const std::string& key)
  {
    bool result = false;
    std::string value = getString(section, key);
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

  static ConfigurationSPtr instance_;

  std::map<std::string, std::string> properties_;
#endif /* DOXYGEN_RUNNING */
};

#ifndef DOXYGEN_RUNNING
ConfigurationSPtr Configuration::instance_ = ConfigurationSPtr();
#endif /* DOXYGEN_RUNNING */

} // namespace Straight_skeletons_3
} // namespace CGAL

#endif /* CGAL_STRAIGHT_SKELETON_3_INTERNAL_CONFIGURATION_H */
