// Copyright (c) 2019
// GeometryFactory (France).  All rights reserved.
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
// Author(s)     : Maxime Gimeno

#ifndef CGAL_INTERSECTIONS_BBOX_2_BBOX_2_H
#define CGAL_INTERSECTIONS_BBOX_2_BBOX_2_H

#include <CGAL/Bbox_2.h>

namespace CGAL {

bool
inline
do_intersect(const CGAL::Bbox_2& c,
             const CGAL::Bbox_2& bbox)
{
  return CGAL::do_overlap(c, bbox);
}

typename boost::optional< typename boost::variant<Bbox_2> >
inline
intersection(const CGAL::Bbox_2& a,
             const CGAL::Bbox_2& b)
{
  typedef typename boost::variant<Bbox_2> variant_type;
  typedef typename boost::optional<variant_type> result_type;

  if(!do_intersect(a, b))
    return result_type();

  double xmin = (std::max)(a.xmin(), b.xmin());
  double xmax = (std::min)(a.xmax(), b.xmax());

  double ymin = (std::max)(a.ymin(), b.ymin());
  double ymax = (std::min)(a.ymax(), b.ymax());

  return result_type(std::forward<Bbox_2>(Bbox_2(xmin, ymin, xmax, ymax)));
}

} // namespace CGAL

#endif // CGAL_INTERSECTIONS_BBOX_2_BBOX_2_H
