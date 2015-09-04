#ifndef CGAL_MUTEX_H
#define CGAL_MUTEX_H

#include <CGAL/config.h>

#ifdef BOOST_MSVC
#include <mutex>
#define CGAL_MUTEX std::mutex
#define CGAL_SCOPED_LOCK(M) std::unique_lock<std::mutex> scoped_lock(M)

#else
#include <boost/thread/mutex.hpp>
#define CGAL_MUTEX boost::mutex
#define CGAL_SCOPED_LOCK(M) boost::mutex::scoped_lock scoped_lock(M)

#endif

#endif // CGAL_MUTEX_H
