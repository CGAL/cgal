/****************************************************************************
 * Core Library Version 1.7, August 2004
 * Copyright (c) 1995-2004 Exact Computation Project
 * All rights reserved.
 *
 * This file is part of CGAL (www.cgal.org).
 *
 * File: Expr.h
 * Synopsis: a class of Expression in Level 3
 *
 * Written by
 *       Koji Ouchi <ouchi@simulation.nyu.edu>
 *       Chee Yap <yap@cs.nyu.edu>
 *       Igor Pechtchanski <pechtcha@cs.nyu.edu>
 *       Vijay Karamcheti <vijayk@cs.nyu.edu>
 *       Chen Li <chenli@cs.nyu.edu>
 *       Zilin Du <zilin@cs.nyu.edu>
 *       Sylvain Pion <pion@cs.nyu.edu>
 *       Vikram Sharma<sharma@cs.nyu.edu>
 *
 * WWW URL: http://cs.nyu.edu/exact/
 * Email: exact@cs.nyu.edu
 *
 * $URL$
 * $Id$
 * SPDX-License-Identifier: LGPL-3.0-or-later
 ***************************************************************************/

// Author(s)     : Andreas Fabri

#ifndef CGAL_CORE_EXPORT_H
#define CGAL_CORE_EXPORT_H

#include <CGAL/config.h>
#include <CGAL/export/helpers.h>

// If CGAL_EXPORTS is defined, one are building the CGAL library, and we do
// not want artificial dll-imports of Core symbols (because of
// auto-linking).
#if ( ! defined(CGAL_EXPORTS) ) && defined(CGAL_BUILD_SHARED_LIBS) \
  && ( ! defined(CGAL_HEADER_ONLY) )

#  if defined(CGAL_Core_EXPORTS) // defined by CMake or in cpp files of the dll

#    define CGAL_CORE_EXPORT CGAL_DLL_EXPORT
#    define CGAL_CORE_EXPIMP_TEMPLATE

#  else // not CGAL_Core_EXPORTS

#    define CGAL_CORE_EXPORT CGAL_DLL_IMPORT
#    define CGAL_CORE_EXPIMP_TEMPLATE extern

#  endif // not CGAL_CORE_EXPORTS

#else // not CGAL_BUILD_SHARED_LIBS

#  define CGAL_CORE_EXPORT
#  define CGAL_CORE_EXPIMP_TEMPLATE

#endif // not CGAL_BUILD_SHARED_LIBS

#endif //  CGAL_CORE_EXPORT_H


