#ifndef CGAL_TSS_H
#define CGAL_TSS_H

#include <CGAL/config.h>

#if defined( CGAL_HAS_THREADS )
  #ifdef  CGAL_CAN_USE_CXX11_THREAD_LOCAL
    #include <thread>
    #define CGAL_THREAD_LOCAL thread_local
  #else
    #ifdef BOOST_MSVC
      #include <thread>
      #define CGAL_THREAD_LOCAL  __declspec( thread )
    #else
      #include <boost/thread/tss.hpp>
      #define CGAL_THREAD_LOCAL thread_local
    #endif
  #endif
#endif

#endif // CGAL_TSS_H
