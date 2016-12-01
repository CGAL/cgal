// Copyright (c) 2013-2015  The University of Western Sydney, Australia.
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
// Authors: Weisheng Si, Quincy Tse

/*! \file Cone_spanners_2_enum.h
 *
 * This header defines enumerators for the cone spanners functors.
 */

#ifndef CONE_SPANNERS_2_ENUM_H
#define CONE_SPANNERS_2_ENUM_H

namespace CGAL {
  /*! \ingroup PkgConeBasedSpanners

   \brief An enum of the types of cone spanners.
   */
  enum Half { EVEN_CONES = 0, ODD_CONES = 1, ALL_CONES = 2 };

}  // namespace CGAL

#endif
