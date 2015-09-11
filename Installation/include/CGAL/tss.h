#ifndef CGAL_TSS_H
#define CGAL_TSS_H

#include <CGAL/config.h>

#if defined( CGAL_HAS_THREADS )
  #ifdef  CGAL_CAN_USE_CXX11_THREAD_LOCAL
    #pragma message ( "Use keyword thread_local" )
    #include <thread>
    #define CGAL_THREAD_LOCAL thread_local
  #else
    #ifdef BOOST_MSVC
      #include <thread>
      #pragma message ("Use __declspec( thread )" )
      #define CGAL_THREAD_LOCAL  __declspec( thread )
    #else
      #pragma message ("Use thread_local from boost")
      #define CGAL_USE_BOOST_THREAD
      #include <boost/thread/tss.hpp>
      #define CGAL_THREAD_LOCAL thread_local
    #endif
  #endif

#ifdef CGAL_USE_BOOST_THREAD

#define CGAL_THREAD_LOCAL_DECLARE(TYPE, VAR)  static boost::thread_specific_ptr<TYPE> VAR##_ptr

#define CGAL_THREAD_LOCAL_DECLARE_PTR(TYPE, VAR)  static boost::thread_specific_ptr<TYPE> VAR

#define CGAL_THREAD_LOCAL_IS_UNINITIALIZED(X) X.get() == NULL

#define CGAL_THREAD_LOCAL_INITIALIZE_PTR(VAR, VAL) if(VAR.get() == NULL) {VAR.reset(VAL);}

#define CGAL_THREAD_LOCAL_INITIALIZE(TYPE, VAR, VAL) if(VAR##_ptr.get() == NULL) {VAR##_ptr.reset(new TYPE(VAL));}

#define CGAL_THREAD_LOCAL_SET(VAR, VAL) VAR.reset(VAL)

#define CGAL_THREAD_LOCAL_GET_PTR(VAR) VAR.get()

#define CGAL_THREAD_LOCAL_GET(TYPE, VAR) const TYPE& VAR =  * VAR##_ptr.get()

#else

#define CGAL_THREAD_LOCAL_DECLARE(TYPE, VAR)  static CGAL_THREAD_LOCAL TYPE* VAR##_ptr = NULL

#define CGAL_THREAD_LOCAL_DECLARE_PTR(TYPE, VAR)  static CGAL_THREAD_LOCAL TYPE* VAR = NULL

#define CGAL_THREAD_LOCAL_IS_UNINITIALIZED(X) X == NULL

#define CGAL_THREAD_LOCAL_INITIALIZE_PTR(VAR, VAL) if(VAR == NULL) { VAR = VAL; }

#define CGAL_THREAD_LOCAL_INITIALIZE(TYPE, VAR, VAL) if(VAR##_ptr == NULL) { VAR##_ptr = new TYPE(VAL); }

#define CGAL_THREAD_LOCAL_SET(VAR, VAL) VAR = VAL

#define CGAL_THREAD_LOCAL_GET_PTR(VAR) VAR

#define CGAL_THREAD_LOCAL_GET(TYPE, VAR) const  TYPE& VAR = *VAR##_ptr
#endif

#else 

#define CGAL_THREAD_LOCAL_INITIALIZE(TYPE, VAR, VAL) static const TYPE VAR = VAL;  
#define CGAL_THREAD_LOCAL_DECLARE(TYPE, VAR)
#define CGAL_THREAD_LOCAL_GET(TYPE, VAR)

#endif

#endif // CGAL_TSS_H
