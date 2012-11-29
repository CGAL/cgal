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

#ifndef DELETE_CURVE_MODE_ITEM_EDITOR_H
#define DELETE_CURVE_MODE_ITEM_EDITOR_H

#include <QComboBox>
#include "DeleteCurveMode.h"

class QWidget;

class DeleteCurveModeItemEditor : public QComboBox
{
  Q_OBJECT
  Q_PROPERTY( DeleteCurveMode mode READ mode WRITE setMode USER true )

public:
  DeleteCurveModeItemEditor( QWidget* parent = 0 );

public:
  DeleteCurveMode mode( ) const;
  void setMode( DeleteCurveMode m );

}; // class DeleteCurveModeItemEditor

#endif // DELETE_CURVE_MODE_ITEM_EDITOR_H
