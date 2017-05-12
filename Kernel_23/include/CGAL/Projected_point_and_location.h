// Copyright (c) 1999,2002,2005
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany)
// and Tel-Aviv University (Israel).  All rights reserved.
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
//
//
// Author(s)     : Andreas Fabri


#ifndef CGAL_PROJECTED_POINT_AND_LOCATION_H
#define CGAL_PROJECTED_POINT_AND_LOCATION_H

namespace CGAL {
  
template <typename Point>
struct Projected_point_and_location {
  Point projected_point;
  unsigned int dimension;
  unsigned int index;
  operator Point() const { return projected_point; }
};
  

} // namespace

#endif // CGAL_PROJECTED_POINT_AND_LOCATION_H
