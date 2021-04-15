/****************************************************************************
 * Core Library Version 1.7, August 2004
 * Copyright (c) 1995-2004 Exact Computation Project
 * All rights reserved.
 *
 * This file is part of CGAL (www.cgal.org).
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
 * SPDX-License-Identifier: LGPL-3.0-or-later
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
