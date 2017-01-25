// Copyright (c) 2009-2011  GeometryFactory Sarl (France)
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
// Author(s)     : Laurent RINEAU

#ifndef SCENE_ITEM_CONFIG_H
#define SCENE_ITEM_CONFIG_H

#include <CGAL/license/Three.h>


#include <QtCore/qglobal.h>

#ifdef demo_framework_EXPORTS
#  define scene_item_EXPORTS
#endif

#ifdef scene_item_EXPORTS
#  define SCENE_ITEM_EXPORT Q_DECL_EXPORT
#else
#  define SCENE_ITEM_EXPORT Q_DECL_IMPORT
#endif

#endif // SCENE_ITEM_CONFIG_H
