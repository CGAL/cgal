// Copyright (c) 2001  Tel-Aviv University (Israel).
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
// Author(s)     : Eli Packer <elip@post.tau.ac.il>

#ifndef CGAL_SNAP_ROUNDING_2_TOOLBAR_H
#define CGAL_SNAP_ROUNDING_2_TOOLBAR_H

#include "cgal_types.h"
#include <CGAL/IO/Qt_widget.h>
#include "segment_input_layer_with_snapping.h"

#include <qobject.h>
#include <qtoolbutton.h>
#include <qtoolbar.h>
#include <qbuttongroup.h>
#include <qmainwindow.h>

#include <list>

class Tools_toolbar : public QToolBar
{
  Q_OBJECT
public:
  Tools_toolbar(CGAL::Qt_widget * w, QMainWindow * mw,
    std::list<Segment_2> * l1);
  ~Tools_toolbar(){};

signals:
  void new_object(CGAL::Object);

private:
  QToolButton * but[10];
  QButtonGroup * button_group;
  CGAL::Qt_widget * widget;
  int nr_of_buttons;

  Segment_input_layer<Rep>         segment_layer;
};

#endif
