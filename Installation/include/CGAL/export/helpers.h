// Copyright (c) 2011 GeometryFactory (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Laurent Rineau

#ifndef CGAL_EXPORT_HELPERS_H
#define CGAL_EXPORT_HELPERS_H

#if defined(CGAL_HEADER_ONLY) && ! defined(CGAL_USE_Qt5_RESOURCES)
#  define CGAL_DLL_IMPORT
#  define CGAL_DLL_EXPORT
#  define CGAL_DLL_LOCAL

#else // !CGAL_HEADER_ONLY
#  if defined(_WIN32) || defined(__CYGWIN__)
#    define CGAL_DLL_IMPORT __declspec(dllimport)
#    define CGAL_DLL_EXPORT __declspec(dllexport)
#    define CGAL_DLL_LOCAL
#  else
    #if __GNUC__ >= 4
      #define CGAL_DLL_IMPORT __attribute__ ((visibility ("default")))
      #define CGAL_DLL_EXPORT __attribute__ ((visibility ("default")))
      #define CGAL_DLL_LOCAL  __attribute__ ((visibility ("hidden")))
    #else
      #define CGAL_DLL_IMPORT
      #define CGAL_DLL_EXPORT
      #define CGAL_DLL_LOCAL
    #endif
#  endif

#endif // CGAL_HEADER_ONLY

#endif // CGAL_EXPORT_HELPERS_H
