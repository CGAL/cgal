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

#undef TRACE
#undef TRACEN
#undef TRACEV
#undef CTRACE
#undef CTRACEN
/*
#undef ASSERT
*/

#if CGAL_NEF_DEBUG>0
namespace CGAL {
static int debugthread=3141592;
}
#endif

#if CGAL_NEF_DEBUG>0
#define SETDTHREAD(l) debugthread=l
#else
#define SETDTHREAD(l)
#endif

#if CGAL_NEF_DEBUG>0
#define TRACE(t)   if((debugthread%CGAL_NEF_DEBUG)==0)\
 std::cerr<<" "<<t;std::cerr.flush()
#else
#define TRACE(t) 
#endif

#if CGAL_NEF_DEBUG>0
#define TRACEV(t)  if((debugthread%CGAL_NEF_DEBUG)==0)\
 std::cerr<<" "<<#t<<" = "<<(t)<<std::endl;std::cerr.flush()
#else
#define TRACEV(t) 
#endif

#if CGAL_NEF_DEBUG>0
#define TRACEN(t)  if((debugthread%CGAL_NEF_DEBUG)==0)\
 std::cerr<<" "<<t<<std::endl;std::cerr.flush()
#else
#define TRACEN(t) 
#endif

#if CGAL_NEF_DEBUG>0
#define CTRACE(b,t)  if(b) std::cerr<<" "<<t; else std::cerr<<" 0"
#else
#define CTRACE(b,t) 
#endif

#if CGAL_NEF_DEBUG>0
#define CTRACEN(b,t) if(b) std::cerr<<" "<<t<<"\n"; else std::cerr<<" 0\n"
#else
#define CTRACEN(b,t) 
#endif

#ifndef _ASSERT
#define  CGAL_NEF_ASSERT(cond,fstr) 
#else
#define CGAL_NEF_ASSERT(cond,fstr)   \
  if (!(cond)) {       \
    std::cerr<<"   ASSERT:   "<< #fstr << std::endl; \
    std::cerr<<"   COND:     "<< #cond << std::endl; \
    std::cerr<<"   POSITION: "<<__FILE__<<" at line "<<__LINE__<<std::endl; \
    abort();           \
  }
#endif


#endif //CGAL_DEBUG_H
