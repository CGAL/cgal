// Copyright (c) 2015 GeometryFactory
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_IO_READER_HELPERS_H
#define CGAL_IO_READER_HELPERS_H

#include <CGAL/array.h>
#include <CGAL/assertions.h>
#include <CGAL/Has_member.h>
#include <CGAL/Point_3.h>

#include <boost/mpl/logical.hpp>
#include <boost/utility/enable_if.hpp>
#include <CGAL/Container_helper.h>
#include <CGAL/Point_3.h>

namespace CGAL{
namespace IO {
namespace internal {

// Ideally this should be a std::is_constructible(double, double, double) but boost::is_constructible
// is not safe to use without CXX11
template <typename Kernel>
void fill_point(const double x, const double y, const double z, const double w, CGAL::Point_3<Kernel>& pt)
{
  pt = CGAL::Point_3<Kernel>(x, y, z, w);
}

template <typename Point_3>
void fill_point(const double x, const double y, const double z, const double w, Point_3& pt)
{
  // just in case something weirder than arrays or CGAL points are used as points...
  CGAL::internal::resize(pt, 3);

  pt[0] = x/w; pt[1] = y/w; pt[2] = z/w;
}

} // end namespace internal
} // end namespace IO
} // namespace CGAL

#endif // CGAL_IO_READER_HELPERS_H
