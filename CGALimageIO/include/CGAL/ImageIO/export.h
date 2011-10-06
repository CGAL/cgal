// Copyright (c) 2011 GeometryFactory (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
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

#ifndef CGAL_IMAGEIO_EXPORT_H
#define CGAL_IMAGEIO_EXPORT_H

#include <boost/config.hpp>

#if defined(BOOST_MSVC) && defined(CGAL_BUILD_SHARED_LIB)

#if defined(CGAL_ImageIO_EXPORTS) // add by CMake or in cpp files of the dll
#define	CGAL_IMAGEIO_EXPORT __declspec (dllexport)
#define CGAL_IMAGEIO_EXPIMP_TEMPLATE
#else
#define CGAL_IMAGEIO_EXPORT __declspec (dllimport)
#define CGAL_IMAGEIO_EXPIMP_TEMPLATE extern
#endif
 
#else 

#define  CGAL_IMAGEIO_EXPORT
#define CGAL_IMAGEIO_EXPIMP_TEMPLATE 
#endif

#endif //  CGAL_IMAGEIO_EXPORT_H


