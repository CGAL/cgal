// Copyright (c) 2007 GeometryFactory (France). All rights reserved.
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
// Author(s)     : Fernando Cacciola (fernando.cacciola@geometryfactory.com)

//
// This file is for non-windows compilers running in non-windows platforms (such as cl.exe via cygwin)
// The real file (for MSVC) is in /CGAL/config/msvc/CGAL/auto_link.h
//
#if defined(CGAL_LIB_NAME)
#  undef CGAL_LIB_NAME
#endif

#if defined(CGAL_AUTO_LINK_NOMANGLE)
#  undef CGAL_AUTO_LINK_NOMANGLE
#endif
