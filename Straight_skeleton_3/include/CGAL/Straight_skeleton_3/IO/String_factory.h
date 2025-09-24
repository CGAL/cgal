// Copyright (c) 2024 GeometryFactory (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

/**
 * file   util/StringFactory.h
 * author Gernot Walzl
 * date   2011-12-19
 */

#ifndef CGAL_STRAIGHT_SKELETON_3_INTERNAL_IO_STRING_FACTORY_H
#define CGAL_STRAIGHT_SKELETON_3_INTERNAL_IO_STRING_FACTORY_H

#include <boost/date_time/local_time/local_time.hpp>

#include <sstream>
#include <string>

namespace CGAL {
namespace Straight_skeletons_3 {
namespace IO {

class StringFactory
{
public:
  static std::string fromBoolean(bool value) {
    return (value ? "true" : "false");
  }

  static std::string fromInteger(int value) {
    std::stringstream sstr;
    sstr << value;
    return sstr.str();
  }

  static std::string fromFloat(float value) {
    std::stringstream sstr;
    sstr.precision(7);
    sstr << value;
    return sstr.str();
  }

  static std::string fromFloatArr(int length, float value[]) {
    std::string result;
    result += "<";
    for (int i = 0; i < length; ++i) {
      if (i > 0) {
        result += ", ";
      }
      result += fromFloat(value[i]);
    }
    result += ">";
    return result;
  }

  static std::string fromDouble(double value) {
    std::stringstream sstr;
    sstr.precision(17);
    sstr << value;
    return sstr.str();
  }

  static std::string fromDoubleArr(int length, double value[]) {
    std::string result;
    result += "<";
    for (int i = 0; i < length; ++i) {
      if (i > 0) {
        result += ", ";
      }
      result += fromDouble(value[i]);
    }
    result += ">";
    return result;
  }

  static std::string fromPointer(const void* value) {
    std::stringstream sstr;
    sstr << value;
    return sstr.str();
  }

  static std::string replaceAll(const std::string& str,
                                const std::string& search,
                                const std::string& replace) {
    std::string result = str;
    size_t pos = result.find(search);
    while (pos != std::string::npos) {
      result.replace(pos, search.size(), replace);
      pos = result.find(search, pos+replace.size());
    }
    return result;
  }

  static std::string now(const std::string& format) {
    boost::local_time::local_time_facet* facet =
            new boost::local_time::local_time_facet(format.c_str());
    std::stringstream date_stream;
    date_stream.imbue(std::locale(date_stream.getloc(), facet));
    date_stream << boost::local_time::local_microsec_clock::local_time(
            boost::local_time::time_zone_ptr());
    return date_stream.str();

    // time_t rawtime;
    // time(&rawtime);
    // struct tm * timeinfo = localtime(&rawtime);
    // char result[256];
    // strftime(result, sizeof(result), format.c_str(), timeinfo);
    // return string(result);
  }
};

} // namespace IO
} // namespace Straight_skeletons_3
} // namespace CGAL

#endif /* CGAL_STRAIGHT_SKELETON_3_INTERNAL_IO_STRING_FACTORY_H */
