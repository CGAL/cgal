// Copyright (c) 2016 GeometryFactory (France)
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
// SPDX-License-Identifier: LGPL-3.0+

#ifndef CGAL_MUTEX_H
#define CGAL_MUTEX_H

#include <CGAL/config.h>

#ifdef CGAL_HAS_THREADS
#  ifdef CGAL_CAN_USE_CXX11_MUTEX
#    include <mutex>
#    define CGAL_MUTEX_NS std::
#  else // not CGAL_CAN_USE_CXX11_MUTEX
#    include <boost/thread/mutex.hpp>
#    if BOOST_VERSION < 105300
       // before Boost.Thread 1.53, `boost::lock_guard` was in ../locks.hpp
#      include <boost/thread/locks.hpp>
#    else
#      include <boost/thread/lock_guard.hpp>
#    endif
#    define CGAL_MUTEX_NS boost::
#  endif // not CGAL_CAN_USE_CXX11_MUTEX

   namespace CGAL {
     namespace cpp11 {
       using CGAL_MUTEX_NS mutex;
       using CGAL_MUTEX_NS lock_guard;
       using CGAL_MUTEX_NS unique_lock;
     }
   }

#  define CGAL_MUTEX CGAL::cpp11::mutex
#  define CGAL_SCOPED_LOCK(M) CGAL::cpp11::unique_lock<CGAL::cpp11::mutex> lock(M)

#endif // CGAL_HAS_THREADS

#endif // CGAL_MUTEX_H
