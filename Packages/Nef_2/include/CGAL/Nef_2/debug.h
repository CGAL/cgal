// Copyright (c) 1997-2000  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Michael Seel <seel@mpi-sb.mpg.de>

#ifndef CGAL_DEBUG_H
#define CGAL_DEBUG_H

#include <iostream>

#ifdef NDEBUG
#undef CGAL_NEF_DEBUG
#define CGAL_NEF_DEBUG 0
#endif

#if CGAL_NEF_DEBUG>0
  static int debugthread=1;
#endif

#undef CGAL_NEF_TRACE
#undef CGAL_NEF_TRACEN
#undef CGAL_NEF_TRACEV
#undef CGAL_NEF_CTRACE
#undef CGAL_NEF_CTRACEN
#undef CGAL_NEF_SETDTHREAD

#if CGAL_NEF_DEBUG>0
#define CGAL_NEF_SETDTHREAD(l) debugthread=l
#else
#define CGAL_NEF_SETDTHREAD(l)
#endif

#if CGAL_NEF_DEBUG>0
#define CGAL_NEF_TRACE(t) if((debugthread%CGAL_NEF_DEBUG)==0) \
 std::cerr<<" "<<t; \
 std::cerr.flush()
#else
#define CGAL_NEF_TRACE(t) 
#endif

#if CGAL_NEF_DEBUG>0
#define CGAL_NEF_TRACEV(t) if((debugthread%CGAL_NEF_DEBUG)==0) \
 std::cerr<<" "<<#t<<" = "<<(t)<<std::endl; \
 std::cerr.flush()
#else
#define CGAL_NEF_TRACEV(t) 
#endif

#if CGAL_NEF_DEBUG>0
#define CGAL_NEF_TRACEN(t) if((debugthread%CGAL_NEF_DEBUG)==0) \
 std::cerr<<" "<<t<<std::endl; \
 std::cerr.flush()
#else
#define CGAL_NEF_TRACEN(t) 
#endif

#if CGAL_NEF_DEBUG>0
#define CGAL_NEF_CTRACE(b,t) if(b) std::cerr<<" "<<t; else std::cerr<<" 0"
#else
#define CGAL_NEF_CTRACE(b,t) 
#endif

#if CGAL_NEF_DEBUG>0
#define CGAL_NEF_CTRACEN(b,t) if(b) std::cerr<<" "<<t<<"\n"; else std::cerr<<" 0\n"
#else
#define CGAL_NEF_CTRACEN(b,t) 
#endif

#endif // CGAL_DEBUG_H
