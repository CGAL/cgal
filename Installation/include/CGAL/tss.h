#ifndef CGAL_TSS_H
#define CGAL_TSS_H

#include <CGAL/config.h>

#if defined( CGAL_HAS_THREADS )
#  ifdef  CGAL_CAN_USE_CXX11_THREAD_LOCAL
//#    pragma message ( "Use keyword thread_local" )
#  else
//#    pragma message ("Use thread_local from boost")
#    define CGAL_USE_BOOST_THREAD
#    include <boost/thread/tss.hpp>
#  endif


#  ifdef CGAL_USE_BOOST_THREAD

#    define CGAL_THREAD_LOCAL_VARIABLE(TYPE, VAR,VAL)               \
       boost::thread_specific_ptr<TYPE> VAR##_ptr;                   \
       if(VAR##_ptr.get() == NULL) {VAR##_ptr.reset(new TYPE(VAL));} \
       TYPE& VAR =  * VAR##_ptr.get()

#  else

#    define CGAL_THREAD_LOCAL_VARIABLE(TYPE, VAR,VAL)       \
       thread_local TYPE VAR(VAL)

#  endif

#else 

#  define CGAL_THREAD_LOCAL_VARIABLE(TYPE, VAR,VAL) TYPE VAR = VAL

#endif

#endif // CGAL_TSS_H
