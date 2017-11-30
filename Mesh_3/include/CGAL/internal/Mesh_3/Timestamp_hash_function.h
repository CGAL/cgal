// Copyright (c) 2014-2017 GeometryFactory Sarl (France)
// All rights reserved.
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
//
// Author(s)     : Laurent Rineau,
//                 Mael Rouxel-Labbé

#ifndef CGAL_TIMESTAMP_HASH_FUNCTION_H
#define CGAL_TIMESTAMP_HASH_FUNCTION_H

#include <CGAL/Has_timestamp.h>

#include <boost/functional/hash.hpp>
#include <boost/utility/enable_if.hpp>

#include <cstddef>
#include <iterator>

namespace CGAL {
namespace Mesh_3 {
namespace internal {

// Hash function for unordered sets/maps with keys of type Vertex_handle
template <typename Vertex_handle>
struct Timestamp_hash_function
{
  template <typename VH>
  std::size_t operator()(VH vh,
                         typename boost::enable_if_c<
                           CGAL::internal::Has_timestamp<
                             typename std::iterator_traits<VH>::value_type>::value>::type* = NULL) const
  {
    return vh->time_stamp();
  }

  template <typename VH>
  std::size_t operator()(VH vh,
                         typename boost::disable_if_c<
                           CGAL::internal::Has_timestamp<
                             typename std::iterator_traits<VH>::value_type>::value>::type* = NULL) const
  {
    return boost::hash_value(&*vh);
  }
};

} // namespace internal
} // namespace Mesh_3
} // namespace CGAL

#endif // CGAL_TIMESTAMP_HASH_FUNCTION_H
