// Copyright (c) 2003-2008  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
// Author(s)     : Monique Teillaud, Sylvain Pion

// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (ECG - Effective Computational Geometry for Curves and Surfaces) 
// and a STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)

#ifndef CGAL_IO_QT_WIDGET_CIRCULAR_ARC_ENDPOINT_2_H
#define CGAL_IO_QT_WIDGET_CIRCULAR_ARC_ENDPOINT_2_H

#include <CGAL/IO/Qt_widget.h>
#include <CGAL/Circular_arc_point_2.h>

namespace CGAL {

template < typename CK >
CGAL::Qt_widget &
operator<<(CGAL::Qt_widget & widget, const CGAL::Circular_arc_point_2<CK> &p)
{
  typedef typename CK::Point_2   Point_2;
  return widget << Point_2(to_double(p.x()), to_double(p.y()));
}

} // namespace CGAL

#endif // CGAL_IO_QT_WIDGET_CIRCULAR_ARC_ENDPOINT_2_H
