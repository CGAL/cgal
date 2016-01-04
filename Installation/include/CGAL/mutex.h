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

#ifndef CGAL_MUTEX_H
#define CGAL_MUTEX_H

#include <CGAL/config.h>

#ifdef CGAL_HAS_THREADS
#ifdef CGAL_CAN_USE_CXX11_MUTEX
#include <mutex>
#define CGAL_MUTEX std::mutex
#define CGAL_SCOPED_LOCK(M) std::unique_lock<std::mutex> scoped_lock(M)
#else
#include <boost/thread/mutex.hpp>
#define CGAL_MUTEX boost::mutex
#define CGAL_SCOPED_LOCK(M) boost::mutex::scoped_lock scoped_lock(M)
#endif
#endif
#endif // CGAL_MUTEX_H
