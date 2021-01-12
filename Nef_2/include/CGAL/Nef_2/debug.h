// Copyright (c) 2000  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Michael Seel <seel@mpi-sb.mpg.de>

#ifndef CGAL_NEF_2_DEBUG_H
#define CGAL_NEF_2_DEBUG_H

#include <iostream>

#ifdef NDEBUG
#undef CGAL_NEF_DEBUG
#define CGAL_NEF_DEBUG 1
#endif

#ifndef CGAL_NEF_DEBUG
#define CGAL_NEF_DEBUG 1
#endif

#ifndef NDEBUG
  static int debugthread=1;
#endif

#undef CGAL_NEF_TRACE
#undef CGAL_NEF_TRACEN
#undef CGAL_NEF_TRACEV
#undef CGAL_NEF_CTRACE
#undef CGAL_NEF_CTRACEN
#undef CGAL_NEF_SETDTHREAD

#ifndef NDEBUG
#define CGAL_NEF_SETDTHREAD(l) debugthread=l
#else
#define CGAL_NEF_SETDTHREAD(l)
#endif

#ifndef NDEBUG
#define CGAL_NEF_TRACE(t) if((debugthread%CGAL_NEF_DEBUG)==0) \
 std::cerr<<" "<<t; \
 std::cerr.flush()
#else
#define CGAL_NEF_TRACE(t)
#endif

#ifndef NDEBUG
#define CGAL_NEF_TRACEV(t) if((debugthread%CGAL_NEF_DEBUG)==0) \
 std::cerr<<" "<<#t<<" = "<<(t)<<std::endl; \
 std::cerr.flush()
#else
#define CGAL_NEF_TRACEV(t)
#endif

#ifndef NDEBUG
#define CGAL_NEF_TRACEN(t) if((debugthread%CGAL_NEF_DEBUG)==0) \
 std::cerr<< " "<<t<<std::endl;
#else
#define CGAL_NEF_TRACEN(t)
#endif

#ifndef NDEBUG
#define CGAL_NEF_CTRACE(b,t) if(b) std::cerr<<" "<<t; else std::cerr<<" 0"
#else
#define CGAL_NEF_CTRACE(b,t)
#endif

#ifndef NDEBUG
#define CGAL_NEF_CTRACEN(b,t) if(b) std::cerr<<" "<<t<<"\n"; else std::cerr<<" 0\n"
#else
#define CGAL_NEF_CTRACEN(b,t)
#endif

#endif // CGAL_NEF_2_DEBUG_H
