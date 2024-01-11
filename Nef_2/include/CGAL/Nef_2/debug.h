// Copyright (c) 2000  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Michael Seel <seel@mpi-sb.mpg.de>

#ifndef CGAL_NEF_2_DEBUG_H
#define CGAL_NEF_2_DEBUG_H

#include <CGAL/license/Nef_2.h>


#include <iostream>

#ifdef NDEBUG
#undef CGAL_USE_TRACE
#endif

#ifndef CGAL_NEF_DEBUG
#define CGAL_NEF_DEBUG 1
#endif

#ifdef CGAL_USE_TRACE
  static int debugthread=1;
#endif

#undef CGAL_NEF_TRACE
#undef CGAL_NEF_TRACEN
#undef CGAL_NEF_TRACEV
#undef CGAL_NEF_CTRACE
#undef CGAL_NEF_CTRACEN
#undef CGAL_NEF_SETDTHREAD

#ifdef CGAL_USE_TRACE
#define CGAL_NEF_SETDTHREAD(l) debugthread=l
#else
#define CGAL_NEF_SETDTHREAD(l)
#endif

#ifdef CGAL_USE_TRACE
#define CGAL_NEF_TRACE(t) if((debugthread%CGAL_NEF_DEBUG)==0) \
    { std::cerr<<" "<<t; }
#else
#define CGAL_NEF_TRACE(t) (static_cast<void>(0))
#endif

#ifdef CGAL_USE_TRACE
#define CGAL_NEF_TRACEV(t) if((debugthread%CGAL_NEF_DEBUG)==0) \
    { std::cerr<<" "<<#t<<" = "<<(t)<<std::endl; }
#else
#define CGAL_NEF_TRACEV(t) (static_cast<void>(0))
#endif

#ifdef CGAL_USE_TRACE
#define CGAL_NEF_TRACEN(t) if((debugthread%CGAL_NEF_DEBUG)==0) \
    { std::cerr<< " "<<t<<std::endl; }
#else
#define CGAL_NEF_TRACEN(t) (static_cast<void>(0))
#endif

#ifdef CGAL_USE_TRACE
#define CGAL_NEF_CTRACE(b,t) if(b) {std::cerr<<" "<<t;} else {std::cerr<<" 0"}
#else
#define CGAL_NEF_CTRACE(b,t) (static_cast<void>(0))
#endif

#ifdef CGAL_USE_TRACE
#define CGAL_NEF_CTRACEN(b,t) if(b){ std::cerr<<" "<<t<<"\n";} else {std::cerr<<" 0\n"}
#else
#define CGAL_NEF_CTRACEN(b,t) (static_cast<void>(0))
#endif

#endif // CGAL_NEF_2_DEBUG_H
