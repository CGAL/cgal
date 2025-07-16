// Copyright (c) 2024 GeometryFactory (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

/**
 * @file   debug.h
 * @author Gernot Walzl
 * @date   2011-12-22
 *
 * based on http://www.digitalpeer.com/id/debug
 */

#ifndef DEBUG_H
#define DEBUG_H

#include "config.h"
#include "util/StackTrace.h"

#include <cassert>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#ifdef DEBUG

#define DBGOUT std::cout

/* just a helper for code location */
#define LOC DBGOUT << "[DEBUG] " << __FILE__ << ":" << __LINE__ << ": ";

#ifndef DEBUG_PRINT_VERBOSITY
#  define DEBUG_PRINT_VERBOSITY 1024
#endif

/* macro using var args */
#define DEBUG_CODE(code) code
#define DEBUG_PRINT(m) /*LOC*/ DBGOUT << m << std::endl;
#define DEBUG_PRINT_V(l,m) if (l <= DEBUG_PRINT_VERBOSITY) /*LOC*/ DBGOUT << m << std::endl;
#define DEBUG_PRINT_IF(c,l,m) if ( (c) && l <= DEBUG_PRINT_VERBOSITY) /*LOC*/ DBGOUT << m << std::endl;

/* macros that warn if a smart pointer is invalid */
#define DEBUG_SPTR(sptr) if (!sptr) { LOC DBGOUT << "shared pointer is invalid: " << (#sptr) << std::endl; util::StackTrace::print(DBGOUT); }
#define DEBUG_WPTR(wptr) if (wptr.expired()) { LOC DBGOUT << "weak pointer is expired: " << (#wptr) << std::endl; util::StackTrace::print(DBGOUT); }

#else

/* when debug isn't defined all the macro calls do absolutely nothing */
#define DEBUG_CODE(code)
#define DEBUG_PRINT(m)
#define DEBUG_PRINT_V(l,m)
#define DEBUG_PRINT_IF(c,l,m)
#define DEBUG_SPTR(sptr)
#define DEBUG_WPTR(wptr)

#endif // DEBUG

#endif /* DEBUG_H */
