/****************************************************************************
 * Core Library Version 1.7, August 2004
 * Copyright (c) 1995-2004 Exact Computation Project
 * All rights reserved.
 *
 * This file is part of CGAL (www.cgal.org).
 * You can redistribute it and/or modify it under the terms of the GNU
 * Lesser General Public License as published by the Free Software Foundation,
 * either version 3 of the License, or (at your option) any later version.
 *
 * Licensees holding a valid commercial license may use this file in
 * accordance with the commercial license agreement provided with the
 * software.
 *
 * This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
 * WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 *
 *
 * File: BigFloat.cpp
 * Synopsis:
 *       BigFloat numbers with error bounds 
 *
 *       EXACTNESS PROPERTY:
 *       ==================
 *       For BigFloats that are exact (i.e., error=0),
 *       addition/subtraction and multiplication return the
 *       exact result (i.e., error=0).  We also introduce the operation
 *       div2(), which simply divides a BigFloat by 2,
 *       but this again preserves exactness.  Such exactness
 *       properties are used in our Newton iteration/Sturm Sequences.
 *
 * Written by 
 *       Chee Yap <yap@cs.nyu.edu>
 *       Chen Li <chenli@cs.nyu.edu>
 *       Zilin Du <zilin@cs.nyu.edu>
 *
 * WWW URL: http://cs.nyu.edu/exact/
 * Email: exact@cs.nyu.edu
 *
 * $URL$
 * $Id$
 ***************************************************************************/

#ifndef CGAL_HEADER_ONLY

#include <CGAL/CORE/BigFloat.h>
#include <CGAL/CORE/BigFloat_impl.h>

#endif
