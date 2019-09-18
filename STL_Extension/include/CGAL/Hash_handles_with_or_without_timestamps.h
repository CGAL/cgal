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
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
//
// Author(s)     : Laurent Rineau,
//                 Mael Rouxel-Labbé

#ifndef CGAL_HASH_HANDLES_WITH_OR_WITHOUT_TIMESTAMPS_H
#define CGAL_HASH_HANDLES_WITH_OR_WITHOUT_TIMESTAMPS_H

#include <CGAL/Has_timestamp.h>

#include <boost/functional/hash.hpp>

#include <cstddef>
#include <iterator>

namespace CGAL {

struct Hash_handles_with_or_without_timestamps
{
  template <typename Handle>
  std::size_t operator()(const Handle h) const
  {
    typedef typename std::iterator_traits<Handle>::value_type Type;

    return hash(h, Boolean_tag<CGAL::internal::Has_timestamp<Type>::value>());
  }

  template<typename Handle>
  std::size_t hash(const Handle h, Tag_false) const
  {
    return boost::hash_value(&*h);
  }

  template<typename Handle>
  std::size_t hash(const Handle h, Tag_true) const
  {
    return h->time_stamp();
  }
};

} // namespace CGAL

#endif // CGAL_HASH_HANDLES_WITH_OR_WITHOUT_TIMESTAMPS_H
