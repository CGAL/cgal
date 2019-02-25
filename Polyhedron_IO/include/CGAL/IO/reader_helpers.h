// Copyright (c) 2015 GeometryFactory
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
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_IO_READER_HELPERS_H
#define CGAL_IO_READER_HELPERS_H

#include <CGAL/array.h>
#include <CGAL/assertions.h>
#include <CGAL/Has_member.h>

#include <boost/mpl/logical.hpp>
#include <boost/utility/enable_if.hpp>

namespace CGAL{

namespace IO {

namespace internal {

CGAL_GENERATE_MEMBER_DETECTOR(size);
CGAL_GENERATE_MEMBER_DETECTOR(resize);

// Typical container
template <class Container>
void resize(Container& c, std::size_t size,
            typename boost::enable_if_c<has_resize<Container>::value>::type* = NULL)
{
  c.resize(size);
}

// Container without a resize() function, but with a size() function (e.g. an array)
template <class Container>
void resize(Container& CGAL_assertion_code(array), std::size_t CGAL_assertion_code(size),
            typename boost::enable_if<
              boost::mpl::and_<
                boost::mpl::not_<has_resize<Container> >,
                                 has_size<Container> > >::type* = NULL)
{
  CGAL_assertion(array.size() == size);
}

// A class with neither resize() nor size(), can't enforce size (it better be correct!)
template <class Container>
void resize(Container&, std::size_t,
            typename boost::disable_if<
              boost::mpl::or_<has_resize<Container>,
                              has_size<Container> > >::type* = NULL)
{
}

// Ideally this should be a std::is_constructible(double, double, double) but boost::is_constructible
// is not safe to use without CXX11
template <typename Kernel>
void fill_point(const double x, const double y, const double z, CGAL::Point_3<Kernel>& pt)
{
  pt = CGAL::Point_3<Kernel>(x, y, z);
}

template <typename Point_3>
void fill_point(const double x, const double y, const double z, Point_3& pt)
{
  // just in case something weirder than arrays or CGAL points are used as points...
  resize(pt, 3);

  pt[0] = x; pt[1] = y; pt[2] = z;
}

} // end namespace internal

} // end namespace IO

} // namespace CGAL

#endif // CGAL_IO_READER_HELPERS_H
