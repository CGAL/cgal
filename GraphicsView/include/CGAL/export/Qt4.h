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

#include <CGAL/config.h>
#include <CGAL/export/helpers.h>

#if defined(CGAL_BUILD_SHARED_LIBS)

#  if defined(CGAL_Qt4_EXPORTS) // defined by CMake or in cpp files of the dll

#    define CGAL_QT4_EXPORT CGAL_DLL_EXPORT
#    define CGAL_QT4_EXPIMP_TEMPLATE

#  else // not CGAL_Qt4_EXPORTS

#    define CGAL_QT4_EXPORT CGAL_DLL_IMPORT
#    define CGAL_QT4_EXPIMP_TEMPLATE extern

#  endif // not CGAL_QT4_EXPORTS

#else // not CGAL_BUILD_SHARED_LIBS

#  define CGAL_QT4_EXPORT
#  define CGAL_QT4_EXPIMP_TEMPLATE

#endif // not CGAL_BUILD_SHARED_LIBS

#endif //  CGAL_QT4_EXPORT_H


