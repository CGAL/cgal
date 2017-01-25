// Copyright (c) 2011  GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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
// Author(s)     : Laurent Rineau <Laurent.Rineau@geometryfactory.com>

#ifndef CGAL_QT_RESOURCES_H
#define CGAL_QT_RESOURCES_H

#include <CGAL/license/GraphicsView.h>


#include <CGAL/export/Qt.h>

// cannot use namespaces because of the Q_INIT_RESOURCE macro
CGAL_QT_EXPORT void CGAL_Qt_init_resources();

#define CGAL_QT_INIT_RESOURCES do { CGAL_Qt_init_resources(); } while(0)
// The do{}while(0) trick is used to make that macro value a regular
// statement and not a compound statement.

#ifdef CGAL_HEADER_ONLY
#include <CGAL/Qt/resources_impl.h>
#endif // CGAL_HEADER_ONLY

#endif // CGAL_QT_RESOURCES_H
