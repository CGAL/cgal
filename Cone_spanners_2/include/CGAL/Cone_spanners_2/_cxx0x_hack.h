// Copyright (c) 2013  The University of Western Sydney, Australia.
// All rights reserved.
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
//
// Author(s)     : Quincy Tse, Weisheng Si

/** @file _cxx0x_hack.h
 *
 * This header declares constants to help make code cross-compile between
 * C++11 and prior versions.
 */

#ifndef _CXX0X_HACK_H
#define _CXX0X_HACK_H

#if defined(__GXX_EXPERIMENTAL_CXX0X__) || __cplusplus >= 201103L
#define GXX11
// XXX Hacks for g++-4.5
#define nullptr NULL
#ifndef _GLIBCXX_NOEXCEPT
#define _GLIBCXX_NOEXCEPT throw ()
#endif
#else
#define nullptr NULL
#include <CGAL/result_of.h>
#define _GLIBCXX_NOEXCEPT throw ()
#endif // CXX11

#endif // _CXX0X_HACK_H
