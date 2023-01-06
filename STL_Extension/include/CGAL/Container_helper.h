// Copyright (c) 2014  GeometryFactory (France).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Andreas Fabri

#ifndef CGAL_CONTAINER_HELPER_H
#define CGAL_CONTAINER_HELPER_H

#include <CGAL/array.h>
#include <CGAL/Has_member.h>

#include <CGAL/assertions.h>
#include <boost/mpl/logical.hpp>

#include <cstddef>
#include<type_traits>

namespace CGAL {
namespace internal {

CGAL_GENERATE_MEMBER_DETECTOR(size);
CGAL_GENERATE_MEMBER_DETECTOR(resize);

// Typical container
template <class Container>
void resize(Container& c, std::size_t size,
            std::enable_if_t<has_resize<Container>::value>* = nullptr)
{
  c.resize(size);
}

// Container without a resize() function, but with a size() function (e.g. an array)
template <class Container>
void resize(Container& CGAL_assertion_code(array), std::size_t CGAL_assertion_code(size),
            std::enable_if_t<
              boost::mpl::and_<
                boost::mpl::not_<has_resize<Container> >,
                                 has_size<Container> >::value >* = nullptr)
{
  CGAL_assertion(array.size() == size);
}

// A class with neither resize() nor size(), can't enforce size (it better be correct!)
template <class Container>
void resize(Container&, std::size_t,
            std::enable_if_t<
              !boost::mpl::or_<has_resize<Container>,
                               has_size<Container> >::value >* = nullptr)
{
}

// Plenty of times we ask for a model of RandomAccessContainer and use reserve, but some models
// of the concept do not have a reserve (e.g. deque)
CGAL_GENERATE_MEMBER_DETECTOR(reserve);

// Container with 'reserve'
template <class Container>
void reserve(Container& c, std::size_t size,
             std::enable_if_t<has_reserve<Container>::value>* = nullptr)
{
  c.reserve(size);
}

// Container with no 'reserve'
template <class Container>
void reserve(Container&, std::size_t,
             std::enable_if_t<!has_reserve<Container>::value>* = nullptr)
{
}

} // namespace internal
} // namespace CGAL

#endif // CGAL_CONTAINER_HELPER_H
