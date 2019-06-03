// Copyright (c) 2019
// GeometryFactory (France)
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
// 
//
// Author(s)     : Simon Giraudot
 
#ifndef CGAL_KERNEL_HASH_FUNCTIONS_H
#define CGAL_KERNEL_HASH_FUNCTIONS_H

namespace CGAL
{

template <typename K>
inline std::enable_if_t<std::is_floating_point<typename K::FT>::value, std::size_t>
hash_value (const Point_3<K>& point)
{
  std::size_t result = boost::hash_value(point.x());
  boost::hash_combine(result, boost::hash_value(point.y()));
  boost::hash_combine(result, boost::hash_value(point.z()));
  return result;
}

} //namespace CGAL

// overloads of std::hash used for using std::unordered_[set/map] on CGAL Kernel objects
namespace std
{

template <typename K> struct hash<CGAL::Point_3<K> >
{
  std::size_t operator() (const CGAL::Point_3<K>& point) const
  {
    return CGAL::hash_value<K> (point);
  }
};

}

#endif  // CGAL_KERNEL_HASH_FUNCTIONS_H
