// Copyright (c) 2011 GeometryFactory (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Laurent Rineau

#ifndef CGAL_QT_CONFIG_H
#define CGAL_QT_CONFIG_H

#include <CGAL/license/GraphicsView.h>


#include <QtCore/qglobal.h>

#if defined(CGAL_Qt6_DLL)
#  if defined(CGAL_Qt6_EXPORTS)
#    define CGAL_QT_EXPORT Q_DECL_EXPORT
#  else
#    define CGAL_QT_EXPORT Q_DECL_IMPORT
#  endif
#else
// empty definition
#  define CGAL_QT_EXPORT
#endif

#endif // CGAL_QT_CONFIG_H
