#ifndef CGAL_TSS_H
#define CGAL_TSS_H

#include <CGAL/config.h>

#ifdef CGAL_HAS_THREADS
#ifdef BOOST_MSVC
#include <thread>
#define CGAL_THREAD_LOCAL  __declspec( thread )
#else
# include <boost/thread/tss.hpp>
#define CGAL_THREAD_LOCAL thread_local 
#endif
#endif
#endif // CGAL_MUTEX_H
