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

#include <LEDA/system/basic.h>

#ifdef LEDA_NAMESPACE
#  define CGAL_LEDA_SCOPE  leda
#else
#  define CGAL_LEDA_SCOPE
#endif


#endif


#endif
