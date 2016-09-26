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
// 
//
// Author(s)     : Andreas Fabri
//                 Stefan Schirra

#ifndef CGAL_ORIGIN_H
#define CGAL_ORIGIN_H

#include <CGAL/config.h>

namespace CGAL {

class Origin
{};

class Null_vector
{};

#ifndef CGAL_HEADER_ONLY

CGAL_EXPORT extern const Origin ORIGIN;
CGAL_EXPORT extern const Null_vector NULL_VECTOR;

#endif

} //namespace CGAL

#ifdef CGAL_HEADER_ONLY
#include <CGAL/Origin_impl.h>
#endif // CGAL_HEADER_ONLY

#endif // CGAL_ORIGIN_H
