// Copyright (c) 2012, 2020 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Alex Tsui <alextsui05@gmail.com>
//            Ahmed Essam <theartful.ae@gmail.com>

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
    EDGE_COLOR_KEY,               /*!< color key  */
    EDGE_WIDTH_KEY,               /*!< width key  */
    VERTEX_COLOR_KEY,             /*!< vertex color  */
    VERTEX_RADIUS_KEY,            /*!< vertex radius  */
    ENVELOPE_EDGE_COLOR_KEY,      /*!< envelope color  */
    ENVELOPE_EDGE_WIDTH_KEY,      /*!< envelope size  */
    ENVELOPE_VERTEX_COLOR_KEY,    /*!< envelope vertex color  */
    ENVELOPE_VERTEX_RADIUS_KEY,   /*!< color key  */
    VERTICAL_RAY_EDGE_COLOR_KEY,  /*!< shooting ray color  */
    VERTICAL_RAY_EDGE_WIDTH_KEY,  /*!< shooting ray size  */
    GRID_COLOR_KEY                /*!< color of the grid  */
  };

  ArrangementDemoPropertiesDialog( ArrangementDemoWindow* parent_ = nullptr );
  QVariant property( int index );

protected:
  void setupUi( );
  void updateUi( );

  ArrangementDemoWindow* parent;
  Ui::ArrangementDemoPropertiesDialog* ui;
}; // class ArrangementDemoPropertiesDialog

#endif // ARRANGEMENT_DEMO_PROPERTIES_DIALOG_H
