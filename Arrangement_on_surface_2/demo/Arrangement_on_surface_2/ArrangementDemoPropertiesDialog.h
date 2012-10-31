// Copyright (c) 2012  Tel-Aviv University (Israel).
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
// $URL: $
// $Id: $
//
// Author(s)     : Alex Tsui <alextsui05@gmail.com>

#ifndef ARRANGEMENT_DEMO_PROPERTIES_DIALOG_H
#define ARRANGEMENT_DEMO_PROPERTIES_DIALOG_H

#include <QDialog>

class ArrangementDemoWindow;

namespace Ui
{
  class ArrangementDemoPropertiesDialog;
}

class ArrangementDemoPropertiesDialog : public QDialog
{
  Q_OBJECT
  public:
  // keep this in order with the ui layout
  enum PropertyKey {
    EDGE_COLOR_KEY,
    EDGE_WIDTH_KEY,
    VERTEX_COLOR_KEY,
    VERTEX_RADIUS_KEY,
    ENVELOPE_EDGE_COLOR_KEY,
    ENVELOPE_EDGE_WIDTH_KEY,
    ENVELOPE_VERTEX_COLOR_KEY,
    ENVELOPE_VERTEX_RADIUS_KEY,
    VERTICAL_RAY_EDGE_COLOR_KEY,
    VERTICAL_RAY_EDGE_WIDTH_KEY,
    DELETE_CURVE_MODE_KEY,
    GRID_SIZE_KEY,
    GRID_COLOR_KEY
  };

  ArrangementDemoPropertiesDialog( ArrangementDemoWindow* parent_ = 0,
                                   Qt::WindowFlags f = 0 );
  QVariant property( int index );

protected:
  void setupUi( );
  void updateUi( );
    
  ArrangementDemoWindow* parent;
  Ui::ArrangementDemoPropertiesDialog* ui;
}; // class ArrangementDemoPropertiesDialog

#endif // ARRANGEMENT_DEMO_PROPERTIES_DIALOG_H
