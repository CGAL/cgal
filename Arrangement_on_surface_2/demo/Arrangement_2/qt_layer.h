// Copyright (c) 2005  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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
