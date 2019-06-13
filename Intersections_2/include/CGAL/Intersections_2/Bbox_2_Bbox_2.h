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
#include <CGAL/Iso_rectangle_2.h>
#include <CGAL/kernel_assertions.h>
#include <CGAL/number_utils.h>
#include <CGAL/Intersections_2/Iso_rectangle_2_Iso_rectangle_2.h>

namespace CGAL {

bool
inline
do_intersect(const CGAL::Bbox_2& c,
             const CGAL::Bbox_2& bbox)
{
  return CGAL::do_overlap(c, bbox);
}

typename boost::optional< typename
boost::variant< Bbox_2> >
inline
intersection(const CGAL::Bbox_2& a,
             const CGAL::Bbox_2& b) {

  typedef typename
  boost::variant< Bbox_2> variant_type;
  typedef typename boost::optional< variant_type > Result_type;
  if(!do_intersect(a,b))
  {
    return Result_type();
  }

  double xmin, xmax, ymin, ymax;
  xmin = (std::max)(a.xmin(), b.xmin());
  xmax = (std::min)(a.xmax(), b.xmax());

  ymin = (std::max)(a.ymin(), b.ymin());
  ymax = (std::min)(a.ymax(), b.ymax());

  return Result_type(std::forward<Bbox_2>(Bbox_2(xmin, ymin, xmax, ymax)));
}

}
#endif // BBOX_2_BBOX_2_H
