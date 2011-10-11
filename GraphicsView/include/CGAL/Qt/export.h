// Copyright (c) 2011 GeometryFactory (France). All rights reserved.
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
// 
//
// Author(s)     : Andreas Fabri
#ifndef CGAL_QT4_EXPORT_H
#define CGAL_QT4_EXPORT_H

#include <boost/config.hpp>
#include <CGAL/config.h>

#if defined(BOOST_MSVC) && defined(CGAL_BUILD_SHARED_LIBS)

#if defined(CGAL_Qt4_EXPORTS) // add by CMake or in cpp files of the dll
#define	CGAL_QT4_EXPORT __declspec (dllexport)
#define CGAL_QT4_EXPIMP_TEMPLATE
#else
#define CGAL_QT4_EXPORT __declspec (dllimport)
#define CGAL_QT4_EXPIMP_TEMPLATE extern
#endif
 
#else 

#define  CGAL_QT4_EXPORT
#define CGAL_QT4_EXPIMP_TEMPLATE 
#endif

#endif //  CGAL_QT4_EXPORT_H


