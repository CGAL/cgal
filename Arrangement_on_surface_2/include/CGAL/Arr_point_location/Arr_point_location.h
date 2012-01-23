// Copyright (c) 2012 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: $
// $Id: $
// 
//
// Author(s) : Efi Fogel   <efif@post.tau.ac.il>

#ifndef CGAL_ARR_POINT_LOCATION_H
#define CGAL_ARR_POINT_LOCATION_H

// The macro CGAL_POINT_LOCATION_VERSION controls which version of the
// point location is used. Currently two values are supported:
// 1. Point location with CGAL::Object
// 2. Point location with boost::optional<boost::variant<...> >
// The default value is 2.

#if !defined(CGAL_POINT_LOCATION_VERSION)
#define CGAL_INTERSECTION_VERSION 2
#endif

namespace CGAL {

} //namespace CGAL

#endif
