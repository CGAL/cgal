// Copyright (c) 2006 Tel-Aviv University (Israel).
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
// $URL$
// $Id$
// 
//
// Author(s): Ron Wein          <wein@post.tau.ac.il>

#ifndef CGAL_ARR_ENUM_H
#define CGAL_ARR_ENUM_H

#include <CGAL/enum.h>

CGAL_BEGIN_NAMESPACE

enum Curve_end
{
  MIN_END,
  MAX_END
};

enum Boundary_type
{
  MINUS_INFINITY = -1,
  NO_BOUNDARY,
  PLUS_INFINITY
};

CGAL_END_NAMESPACE

#endif

