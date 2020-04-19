// Copyright (c) 2007 GeometryFactory (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
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
