/****************************************************************************
 * Core Library Version 1.7, August 2004
 * Copyright (c) 1995-2004 Exact Computation Project
 * All rights reserved.
 *
 * This file is part of CGAL (www.cgal.org).
 *
 * $URL$
 * $Id$
 * SPDX-License-Identifier: LGPL-3.0-or-later
 ***************************************************************************/

#ifndef _CORE_CONFIG_H_
#define _CORE_CONFIG_H_

// disable debug
//#define CORE_DISABLE_DEBUG

// disable inline functions
//#define CORE_DISABLE_INLINE

// disable memory pool
//#define CORE_DISABLE_MEMPOOL

// debug reference counting
//#define CORE_RC_DEBUG 1

#include <CGAL/auto_link/CORE.h>

#include <CGAL/export/CORE.h>

#ifdef CGAL_TEST_SUITE
// disabled for the testsuite to avoid `w`
#define CGAL_CORE_warning_msg(X ,Y)
// if (!(X)) CGAL_error_msg(Y)
#else
#define CGAL_CORE_warning_msg(X ,Y) CGAL_warning_msg(X ,Y)
#endif


#endif // _CORE_CONFIG_H_
