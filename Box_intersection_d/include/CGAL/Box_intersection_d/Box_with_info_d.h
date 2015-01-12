// Copyright (c) 2004  Max-Planck-Institute Saarbruecken (Germany).
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
//
// Author(s)     : Lutz Kettner  <kettner@mpi-sb.mpg.de>
//                 Andreas Meyer <ameyer@mpi-sb.mpg.de>

#ifndef CGAL_BOX_INTERSECTION_D_BOX_WITH_INFO_D_H
#define CGAL_BOX_INTERSECTION_D_BOX_WITH_INFO_D_H

#include <CGAL/basic.h>
#include <CGAL/Box_intersection_d/Box_d.h>

namespace CGAL {

namespace Box_intersection_d {

template<class NT_, int N, class Info_> 
class Box_with_info_d : public Box_d< NT_, N, ID_FROM_BOX_ADDRESS> {
protected:
    Info_ m_info;
public:
    typedef Box_d< NT_, N, ID_FROM_BOX_ADDRESS> Base;
    typedef NT_                      NT;
    typedef Info_                     Info;
    Box_with_info_d() {}
    Box_with_info_d( Info h) : m_info(h) {}
    Box_with_info_d( bool complete, Info h): Base(complete), m_info(h) {}
    Box_with_info_d(NT l[N], NT h[N], Info n) : Base( l, h), m_info(n) {}
    Box_with_info_d( const Bbox_2& b, Info h) : Base( b), m_info(h) {}
    Box_with_info_d( const Bbox_3& b, Info h) : Base( b), m_info(h) {}
    Info info() const { return m_info; }
};

} // end namespace Box_intersection_d


} //namespace CGAL

#endif
