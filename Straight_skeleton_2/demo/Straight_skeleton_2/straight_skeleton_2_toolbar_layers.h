// Copyright (c) 2005, 2006 Fernando Luis Cacciola Carballal. All rights reserved.
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
// Author(s)     : Fernando Cacciola <fernando_cacciola@ciudad.com.ar>
//

#ifndef CGAL_STRAIGHTSKELETON_2_TOOLBAR_LAYERS_H
#define CGAL_STRAIGHTSKELETON_2_TOOLBAR_LAYERS_H

#include <CGAL/IO/Qt_widget.h>

#include <qobject.h>
#include <qmainwindow.h>
#include <qtoolbutton.h>
#include <qtoolbar.h>
#include <qstatusbar.h>
#include <qbuttongroup.h>


template <class T> class Qt_layer_show_regions;
template <class T> class Qt_layer_show_skeleton;
class Qt_layer_show_progress ;

class Layers_toolbar : public QToolBar
{
  Q_OBJECT
public:
  Layers_toolbar( CGAL::Qt_widget*        w
                , QMainWindow*            mw
                , demo::Regions  const&   in
                , demo::SSkelPtr const&   sskel
                , demo::Regions  const&   out
                );
  ~Layers_toolbar();

  Qt_layer_show_progress& get_progress_layer() { return *progress ; }

private:
  QToolButton         *but[4];
  CGAL::Qt_widget     *widget;
  QMainWindow         *window;
  QButtonGroup        *button_group;
  int                 nr_of_buttons;

  Qt_layer_show_regions <demo::Regions> *showI;
  Qt_layer_show_skeleton<demo::SSkel>   *showSSkel;
  Qt_layer_show_regions <demo::Regions> *showO;
  Qt_layer_show_progress                *progress ;

};//end class

#endif
