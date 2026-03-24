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

#include <CGAL/license/Straight_skeleton_3.h>

#include <boost/date_time/local_time/local_time.hpp>

#include <sstream>
#include <string>

namespace CGAL {
namespace Straight_skeletons_3 {
namespace IO {

class String_factory
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

  static std::string fromDouble(double value) {
    std::stringstream sstr;
    sstr.precision(17);
    sstr << value;
    return sstr.str();
  }

  static std::string fromPointer(const void* value) {
    std::stringstream sstr;
    sstr << value;
    return sstr.str();
  }
};

} // namespace IO
} // namespace Straight_skeletons_3
} // namespace CGAL

#endif /* CGAL_STRAIGHT_SKELETON_3_INTERNAL_IO_STRING_FACTORY_H */
