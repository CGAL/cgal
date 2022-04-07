// Copyright (c) 2016 GeometryFactory (France)
//  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial

#ifndef CGAL_TSS_H
#define CGAL_TSS_H

#include <CGAL/config.h>

#if defined( CGAL_HAS_THREADS )

#  define CGAL_STATIC_THREAD_LOCAL_VARIABLE_0(TYPE, VAR)       \
  static thread_local TYPE VAR

#  define CGAL_STATIC_THREAD_LOCAL_VARIABLE(TYPE, VAR, ARG1)       \
  static thread_local TYPE VAR(ARG1)

#  define CGAL_STATIC_THREAD_LOCAL_VARIABLE_2(TYPE, VAR, ARG1, ARG2)       \
  static thread_local TYPE VAR(ARG1,ARG2)

#  define CGAL_STATIC_THREAD_LOCAL_VARIABLE_3(TYPE, VAR, ARG1, ARG2, ARG3) \
  static thread_local TYPE VAR(ARG1,ARG2,ARG3)

#  define CGAL_STATIC_THREAD_LOCAL_VARIABLE_4(TYPE, VAR, ARG1, ARG2, ARG3, ARG4) \
  static thread_local TYPE VAR(ARG1,ARG2,ARG3,ARG4)

#else // not CGAL_HAS_THREADS

#  define CGAL_STATIC_THREAD_LOCAL_VARIABLE_0(TYPE, VAR) static TYPE VAR

#  define CGAL_STATIC_THREAD_LOCAL_VARIABLE(TYPE, VAR,ARG1) static TYPE VAR(ARG1)

#  define CGAL_STATIC_THREAD_LOCAL_VARIABLE_2(TYPE, VAR,ARG1,ARG2) static TYPE VAR(ARG1,ARG2)

#  define CGAL_STATIC_THREAD_LOCAL_VARIABLE_3(TYPE, VAR,ARG1,ARG2,ARG3) static TYPE VAR(ARG1,ARG2,ARG3)

#  define CGAL_STATIC_THREAD_LOCAL_VARIABLE_4(TYPE, VAR,ARG1,ARG2,ARG3,ARG4) static TYPE VAR(ARG1,ARG2,ARG3,ARG4)
#endif // not CGAL_HAS_THREADS

#endif // CGAL_TSS_H
