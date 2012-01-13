// Copyright (c) 1997-2000  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Michael Seel <seel@mpi-sb.mpg.de>
//
#ifndef CGAL_KERNEL_D_DEBUG_H
#define CGAL_KERNEL_D_DEBUG_H

#include <iostream>
#include <string>
#include <sstream>

#undef CGAL_KD_TRACE
#undef CGAL_KD_TRACEN
#undef CGAL_KD_TRACEV
#undef CGAL_KD_CTRACE
#undef CGAL_KD_CTRACEN
#undef CGAL_KD_ASSERT

#if CGAL_KD_DEBUG>0
namespace CGAL {
namespace Kernel_d_internal {
static int debugthread=3141592;
} // Kernel_d_internal 
} // CGAL
#endif

#if CGAL_KD_DEBUG>0
#define CGAL_KD_SETDTHREAD(l) CGAL::Kernel_d_internal::debugthread=l
#else
#define CGAL_KD_SETDTHREAD(l)
#endif

#if CGAL_KD_DEBUG>0
#define CGAL_KD_TRACE(t)   if((CGAL::Kernel_d_internal::debugthread%CGAL_KD_DEBUG)==0)\
 std::cerr<<" "<<t;std::cerr.flush()
#else
#define CGAL_KD_TRACE(t) 
#endif

#if CGAL_KD_DEBUG>0
#define CGAL_KD_TRACEV(t)  if((CGAL::Kernel_d_internal::debugthread%CGAL_KD_DEBUG)==0)\
 std::cerr<<" "<<#t<<" = "<<(t)<<std::endl;std::cerr.flush()
#else
#define CGAL_KD_TRACEV(t) 
#endif

#if CGAL_KD_DEBUG>0
#define CGAL_KD_TRACEN(t)  if((CGAL::Kernel_d_internal::debugthread%CGAL_KD_DEBUG)==0)\
 std::cerr<<" "<<t<<std::endl;std::cerr.flush()
#else
#define CGAL_KD_TRACEN(t) 
#endif

#if CGAL_KD_DEBUG>0
#define CGAL_KD_CTRACE(b,t)  if(b) std::cerr << " " << t; else std::cerr << " 0"
#else
#define CGAL_KD_CTRACE(b,t) 
#endif

#if CGAL_KD_DEBUG>0
#define CGAL_KD_CTRACEN(b,t)  if(b) std::cerr<< " " <<t<<"\n"; else std::cerr<<" 0\n"
#else
#define CGAL_KD_CTRACEN(b,t) 
#endif

#ifndef CGAL_KD__ASSERT
#define  CGAL_KD_ASSERT(cond,fstr) 
#else
#define CGAL_KD_ASSERT(cond,fstr)   \
  if (!(cond)) {       \
    std::cerr << "   ASSERT:   " << #fstr << std::endl; \
    std::cerr << "   COND:     " << #cond << std::endl; \
    std::cerr << "   POSITION: " << __FILE__ << " at line "<< __LINE__ \
              << std::endl; \
    CGAL_error(); \
  }
#endif

#endif //CGAL_KERNEL_D_DEBUG_H
