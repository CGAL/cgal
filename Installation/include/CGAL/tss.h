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

#    define CGAL_STATIC_THREAD_LOCAL_VARIABLE(TYPE, VAR, ARG1)               \
       static boost::thread_specific_ptr<TYPE> VAR##_ptr;                   \
       if(VAR##_ptr.get() == NULL) {VAR##_ptr.reset(new TYPE(ARG1));} \
       TYPE& VAR =  * VAR##_ptr.get()

#  else

#    define CGAL_STATIC_THREAD_LOCAL_VARIABLE(TYPE, VAR, ARG1)       \
       static thread_local TYPE VAR(ARG1)

#  endif

#else 

#  define CGAL_STATIC_THREAD_LOCAL_VARIABLE(TYPE, VAR,ARG1) static TYPE VAR(ARG1)

#endif

#endif // CGAL_TSS_H
