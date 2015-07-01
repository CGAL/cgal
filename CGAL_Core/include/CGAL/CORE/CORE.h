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
 * File: CORE.h
 * Synopsis:
 *      The main inclusion file for the Core Library system.
 *      All "Core programs" must include this file.
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

#ifndef _CORE_CORE_H_
#define _CORE_CORE_H_

#include <CGAL/CORE/CoreDefs.h>
#include <CGAL/CORE/Timer.h>

// User can still access machine types:
typedef double machine_double;
typedef long machine_long;

#ifndef CORE_LEVEL
#   define CORE_LEVEL  DEFAULT_CORE_LEVEL
#endif

#if CORE_LEVEL  == 1
#   define Real double
#   define Expr double
#elif CORE_LEVEL  == 2
#   include <CGAL/CORE/Real.h>
#   undef long
#   undef double
#   define long Real
#   define double Real
#   define Expr Real
#elif CORE_LEVEL  == 3
#   include <CGAL/CORE/Expr.h>
#   undef long
#   undef double
#   define long Expr
#   define double Expr
#   define Real Expr
#elif CORE_LEVEL == 4
#   include <CGAL/CORE/Expr.h>
#endif

#endif // _CORE_CORE_H_
