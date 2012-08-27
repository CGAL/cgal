// Copyright (c) 2007 GeometryFactory (France). All rights reserved.
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
// Author(s)     : Fernando Cacciola (fernando.cacciola@geometryfactry.com)

#ifndef CGAL_AUTO_LINK_IMAGEIO_H
#define CGAL_AUTO_LINK_IMAGEIO_H

#include <CGAL/config.h>

#ifndef CGAL_NO_AUTOLINK_IMAGEIO
#if ( ! defined( CGAL_EXPORTS ) ) && (! defined ( CGAL_ImageIO_EXPORTS ) ) 

// If CGAL_EXPORTS is defined it means that we are building the CGAL
// library as a DLL. The CGAL.dll does not really depend on CGAL_ImageIO,
// whatever the header inclusion graph says.

#define CGAL_LIB_NAME CGAL_ImageIO
#include <CGAL/auto_link/auto_link.h>

#endif // CGAL_EXPORTS
#endif // CGAL_NO_AUTOLINK_IMAGEIO

#endif // CGAL_AUTO_LINK_IMAGEIO_H
