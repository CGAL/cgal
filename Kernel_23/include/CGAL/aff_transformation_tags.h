// Copyright (c) 1999  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
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
// SPDX-License-Identifier: LGPL-3.0+
// 
//
// Author(s)     : Andreas Fabri
 

#ifndef CGAL_AFF_TRANSFORMATION_TAGS_H
#define CGAL_AFF_TRANSFORMATION_TAGS_H

#include <CGAL/config.h>

namespace CGAL {

class Translation {};
class Rotation {};
class Scaling {};
class Reflection {};
class Identity_transformation {};

#ifndef CGAL_HEADER_ONLY

CGAL_EXPORT extern const Translation              TRANSLATION;
CGAL_EXPORT extern const Rotation                 ROTATION;
CGAL_EXPORT extern const Scaling                  SCALING;
CGAL_EXPORT extern const Reflection               REFLECTION;
CGAL_EXPORT extern const Identity_transformation  IDENTITY;

#endif

} //namespace CGAL
#ifdef CGAL_HEADER_ONLY
#include <CGAL/aff_transformation_tags_impl.h>
#endif // CGAL_HEADER_ONLY

#endif // CGAL_AFF_TRANSFORMATION_TAGS_H
