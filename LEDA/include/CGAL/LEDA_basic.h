// Copyright (c) 1999
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Matthias Baesken



#ifndef CGAL_LEDA_BASIC_H
#define CGAL_LEDA_BASIC_H

#include <CGAL/config.h>

#ifdef CGAL_USE_LEDA
// The following is needed for LEDA 4.4 due to min/max problems...
#  define LEDA_NO_MIN_MAX_TEMPL

#ifdef CGAL_CXX20
#  define STREAM_DUMMY // disable stream_dummy() function that is not used and using features removed from c++20
// We cannot undef STREAM_DUMMY as LEDA/internal/PREAMBULE.h is not protected
#endif
#include <LEDA/system/basic.h>



#ifdef LEDA_NAMESPACE
#  define CGAL_LEDA_SCOPE  leda
#else
#  define CGAL_LEDA_SCOPE
#endif


#endif


#endif
