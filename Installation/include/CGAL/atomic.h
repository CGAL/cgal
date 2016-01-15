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

#ifndef CGAL_ATOMIC_H
#define CGAL_ATOMIC_H

#include <CGAL/config.h>

#ifdef CGAL_HAS_THREADS
#  ifdef CGAL_CAN_USE_CXX11_ATOMIC
#    include <atomic>

     namespace CGAL {
       namespace cpp11 {
         using std::atomic;
       }
     }

#  else // not CGAL_CAN_USE_CXX11_ATOMIC
#    include <boost/thread/atomic.hpp>

     namespace CGAL {
       namespace cpp11 {
         using boost::atomic;
       }
     }

#  endif // not CGAL_CAN_USE_CXX11_ATOMIC
#endif // CGAL_HAS_THREADS

#endif // CGAL_ATOMIC_H

