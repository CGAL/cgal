// Copyright (c) 2005  Tel-Aviv University (Israel).
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
// $URL: svn+ssh://scm.gforge.inria.fr/svn/cgal/branches/features/gsoc2012-Arrangement_on_surface_2-demo-atsui/Arrangement_on_surface_2/demo/Arrangement_on_surface_2/qt_layer.h $
// $Id: qt_layer.h 67117 2012-01-13 18:14:48Z lrineau $
//
//
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>

#ifndef CGAL_QT_LAYER_H
#define CGAL_QT_LAYER_H

/*! class Qt_layer is the main layer in the program.
 *  all the tab widget are attached to it.
 */
#include <CGAL/IO/Qt_widget_layer.h>

#include "cgal_types.h"

class QTabWidget;

class Qt_layer : public CGAL::Qt_widget_layer
{
public:
    Qt_layer( QTabWidget * );
	void draw();

private:
	QTabWidget *myBar;
};

#endif
