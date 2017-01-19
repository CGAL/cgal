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
#    define CGAL_STATIC_THREAD_LOCAL_USE_BOOST 1

#    define CGAL_STATIC_THREAD_LOCAL_VARIABLE_0(TYPE, VAR)               \
       static boost::thread_specific_ptr<TYPE> VAR##_ptr;                   \
       if(VAR##_ptr.get() == NULL) {VAR##_ptr.reset(new TYPE());} \
       TYPE& VAR =  * VAR##_ptr.get()

#    define CGAL_STATIC_THREAD_LOCAL_VARIABLE(TYPE, VAR, ARG1)               \
       static boost::thread_specific_ptr<TYPE> VAR##_ptr;                   \
       if(VAR##_ptr.get() == NULL) {VAR##_ptr.reset(new TYPE(ARG1));} \
       TYPE& VAR =  * VAR##_ptr.get()

#    define CGAL_STATIC_THREAD_LOCAL_VARIABLE_2(TYPE, VAR, ARG1, ARG2)       \
       static boost::thread_specific_ptr<TYPE> VAR##_ptr;                   \
       if(VAR##_ptr.get() == NULL) {VAR##_ptr.reset(new TYPE(ARG1,ARG2));}  \
       TYPE& VAR =  * VAR##_ptr.get()

#    define CGAL_STATIC_THREAD_LOCAL_VARIABLE_3(TYPE, VAR, ARG1, ARG2, ARG3)    \
       static boost::thread_specific_ptr<TYPE> VAR##_ptr;                       \
       if(VAR##_ptr.get() == NULL) {VAR##_ptr.reset(new TYPE(ARG1,ARG2,ARG3));} \
       TYPE& VAR =  * VAR##_ptr.get()

#    define CGAL_STATIC_THREAD_LOCAL_VARIABLE_4(TYPE, VAR, ARG1, ARG2, ARG3, ARG4)   \
       static boost::thread_specific_ptr<TYPE> VAR##_ptr;                            \
       if(VAR##_ptr.get() == NULL) {VAR##_ptr.reset(new TYPE(ARG1,ARG2,ARG3,ARG4));} \
       TYPE& VAR =  * VAR##_ptr.get()



#  else // not CGAL_USE_BOOST_THREAD, -> use C++11 thread_local

#    define CGAL_STATIC_THREAD_LOCAL_VARIABLE_0(TYPE, VAR)       \
       static thread_local TYPE VAR

#    define CGAL_STATIC_THREAD_LOCAL_VARIABLE(TYPE, VAR, ARG1)       \
       static thread_local TYPE VAR(ARG1)

#    define CGAL_STATIC_THREAD_LOCAL_VARIABLE_2(TYPE, VAR, ARG1, ARG2)       \
  static thread_local TYPE VAR(ARG1,ARG2)

#    define CGAL_STATIC_THREAD_LOCAL_VARIABLE_3(TYPE, VAR, ARG1, ARG2, ARG3) \
  static thread_local TYPE VAR(ARG1,ARG2,ARG3)

#    define CGAL_STATIC_THREAD_LOCAL_VARIABLE_4(TYPE, VAR, ARG1, ARG2, ARG3, ARG4) \
  static thread_local TYPE VAR(ARG1,ARG2,ARG3,ARG4)

#  endif // not CGAL_USE_BOOST_THREAD

#else // not CGAL_HAS_THREADS

#  define CGAL_STATIC_THREAD_LOCAL_VARIABLE_0(TYPE, VAR) static TYPE VAR

#  define CGAL_STATIC_THREAD_LOCAL_VARIABLE(TYPE, VAR,ARG1) static TYPE VAR(ARG1)

#  define CGAL_STATIC_THREAD_LOCAL_VARIABLE_2(TYPE, VAR,ARG1,ARG2) static TYPE VAR(ARG1,ARG2)

#  define CGAL_STATIC_THREAD_LOCAL_VARIABLE_3(TYPE, VAR,ARG1,ARG2,ARG3) static TYPE VAR(ARG1,ARG2,ARG3)

#  define CGAL_STATIC_THREAD_LOCAL_VARIABLE_4(TYPE, VAR,ARG1,ARG2,ARG3,ARG4) static TYPE VAR(ARG1,ARG2,ARG3,ARG4)
#endif // not CGAL_HAS_THREADS

#endif // CGAL_TSS_H
