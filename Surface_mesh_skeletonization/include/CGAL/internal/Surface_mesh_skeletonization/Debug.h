// Copyright (c) 2013  GeometryFactory (France). All rights reserved.
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
// Author(s)     : Xiang Gao <gaox@ethz.ch>
//

#ifndef CGAL_MCFSKEL_DEBUG_H
#define CGAL_MCFSKEL_DEBUG_H

/// @cond CGAL_DOCUMENT_INTERNAL

/**
 * @file Debug.h
 * @brief This file contains some macro used to hide/show the log for debugging
 * purpose.
 *
 */

// enable debugging output statement
// this is for locating bugs
//#define CGAL_MCFSKEL_DEBUG

// enable info output statement
// this is for tuning parameters
//#define CGAL_MCFSKEL_INFO

#ifdef CGAL_MCFSKEL_DEBUG
#define MCFSKEL_DEBUG(x) x
#else
#define MCFSKEL_DEBUG(x)
#endif

#ifdef CGAL_MCFSKEL_INFO
#define MCFSKEL_INFO(x) x
#else
#define MCFSKEL_INFO(x)
#endif

/// @endcond

#endif // CGAL_MCFSKEL_DEBUG_H
