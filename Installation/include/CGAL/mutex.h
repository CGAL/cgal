// Copyright (c) 2016 GeometryFactory (France)
//  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial

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
