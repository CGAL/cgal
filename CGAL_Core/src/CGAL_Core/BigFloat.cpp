/****************************************************************************
 * Core Library Version 1.7, August 2004
 * Copyright (c) 1995-2004 Exact Computation Project
 * All rights reserved.
 *
 * This file is part of CGAL (www.cgal.org).
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
 * SPDX-License-Identifier: LGPL-3.0-or-later
 ***************************************************************************/

#ifndef CGAL_HEADER_ONLY

#include <CGAL/CORE/BigFloat.h>
#include <CGAL/CORE/BigFloat_impl.h>

#endif
