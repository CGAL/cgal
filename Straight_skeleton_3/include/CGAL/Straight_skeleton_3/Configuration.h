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
 * @file   util/Configuration.h
 * @author Gernot Walzl
 * @date   2012-10-30
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

class Configuration
{
public:
  static ConfigurationSPtr getInstance()
  {
    if (!instance_) {
      instance_ = ConfigurationSPtr(new Configuration());
    }
    return instance_;
  }

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
};

ConfigurationSPtr Configuration::instance_ = ConfigurationSPtr();

} // namespace Straight_skeletons_3
} // namespace CGAL

#endif /* CGAL_STRAIGHT_SKELETON_3_INTERNAL_CONFIGURATION_H */
