// Copyright (c) 1997-2001  
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
// Author(s)     : Sven Schoenherr <sven@inf.ethz.ch>

#ifndef CGAL_OPTIMISATION_DEBUG_H
#define CGAL_OPTIMISATION_DEBUG_H

// macro definitions
// =================

// debug
// -----
#if (    defined( CGAL_OPTIMISATION_NO_DEBUG) \
      || defined( CGAL_NO_DEGUG) || defined( NDEBUG))
#  define  CGAL_optimisation_debug  if ( 0)
#else
#  define  CGAL_optimisation_debug  if ( 1)
#endif // optimisation debug

#endif // CGAL_OPTIMISATION_DEBUG_H

// ===== EOF ==================================================================
