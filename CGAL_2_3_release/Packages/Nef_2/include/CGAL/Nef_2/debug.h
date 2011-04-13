// ============================================================================
//
// Copyright (c) 1997-2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision$
// release_date  : $CGAL_Date$
//
// file          : include/CGAL/Nef_2/debug.h
// package       : Nef_2 
// chapter       : Nef Polyhedra
//
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Michael Seel <seel@mpi-sb.mpg.de>
// maintainer    : Michael Seel <seel@mpi-sb.mpg.de>
// coordinator   : Michael Seel <seel@mpi-sb.mpg.de>
//
// implementation: simple debugging macros
// ============================================================================

#ifndef CGAL_DEBUG_H
#define CGAL_DEBUG_H

#include <iostream>

#undef TRACE
#undef TRACEN
#undef TRACEV
#undef CTRACE
#undef CTRACEN
#undef ASSERT

static int debugthread=3141592;

#if _DEBUG>0
#define SETDTHREAD(l) debugthread=l
#else
#define SETDTHREAD(l)
#endif

#if _DEBUG>0
#define TRACE(t)   if((debugthread%_DEBUG)==0)\
 std::cerr<<" "<<t;std::cerr.flush()
#else
#define TRACE(t) 
#endif

#if _DEBUG>0
#define TRACEV(t)  if((debugthread%_DEBUG)==0)\
 std::cerr<<" "<<#t<<" = "<<(t)<<std::endl;std::cerr.flush()
#else
#define TRACEV(t) 
#endif

#if _DEBUG>0
#define TRACEN(t)  if((debugthread%_DEBUG)==0)\
 std::cerr<<" "<<t<<std::endl;std::cerr.flush()
#else
#define TRACEN(t) 
#endif

#if _DEBUG>0
#define CTRACE(b,t)  if(b) std::cerr<<" "<<t; else std::cerr<<" 0"
#else
#define CTRACE(b,t) 
#endif

#if _DEBUG>0
#define CTRACEN(b,t) if(b) std::cerr<<" "<<t<<"\n"; else std::cerr<<" 0\n"
#else
#define CTRACEN(b,t) 
#endif

#ifndef _ASSERT
#define  ASSERT(cond,fstr) 
#else
#define ASSERT(cond,fstr)   \
  if (!(cond)) {       \
    std::cerr<<"   ASSERT:   "<< #fstr << endl; \
    std::cerr<<"   COND:     "<< #cond << endl; \
    std::cerr<<"   POSITION: "<<__FILE__<<" at line "<<__LINE__<<std::endl; \
    abort();           \
  }
#endif


#endif //CGAL_DEBUG_H


