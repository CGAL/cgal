// Copyright (c) 2016 GeometryFactory Sarl (France)
//  All rights reserved.
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

#ifndef CGAL_ATOMIC_H
#define CGAL_ATOMIC_H

#include <CGAL/config.h>

#ifdef CGAL_HAS_THREADS

#  ifdef CGAL_CAN_USE_CXX11_ATOMIC
#    include <atomic>
#    define CGAL_ATOMIC_NS std
#  else // not CGAL_CAN_USE_CXX11_ATOMIC
#    if BOOST_VERSION >= 105300
#      include <boost/atomic.hpp>
#      define CGAL_ATOMIC_NS boost
#    else // BOOST_VERSION < 105300
#      define CGAL_NO_ATOMIC "Boost.Atomic was introduced in Boost-1.53".
#    endif // BOOST_VERSION < 105300
#  endif // not CGAL_CAN_USE_CXX11_ATOMIC

#  ifndef CGAL_NO_ATOMIC
   namespace CGAL {
     namespace cpp11 {
       using CGAL_ATOMIC_NS ::atomic;

       using CGAL_ATOMIC_NS ::memory_order_relaxed;
       using CGAL_ATOMIC_NS ::memory_order_consume;
       using CGAL_ATOMIC_NS ::memory_order_acquire;
       using CGAL_ATOMIC_NS ::memory_order_release;
       using CGAL_ATOMIC_NS ::memory_order_acq_rel;
       using CGAL_ATOMIC_NS ::memory_order_seq_cst;

       using CGAL_ATOMIC_NS ::atomic_thread_fence;
     }
   }
#  endif // CGAL_ATOMIC_NS
#else
#  define CGAL_NO_ATOMIC "No atomic because CGAL_NO_THREADS is defined."
#endif // CGAL_HAS_THREADS

#endif // CGAL_ATOMIC_H
