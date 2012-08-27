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
// Author(s)     : Matthias Baesken



#ifndef CGAL_LEDA_BASIC_H
#define CGAL_LEDA_BASIC_H

#include <CGAL/config.h>

#ifdef CGAL_USE_LEDA
// The following is needed for LEDA 4.4 due to min/max problems...
#  define LEDA_NO_MIN_MAX_TEMPL

#if CGAL_LEDA_VERSION < 500
#include <LEDA/basic.h>
#else
#include <LEDA/system/basic.h>
#endif

#ifdef LEDA_NAMESPACE
#  define CGAL_LEDA_SCOPE  leda
#else
#  define CGAL_LEDA_SCOPE 
#endif


#endif


#endif
