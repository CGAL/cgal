// Copyright (c) 2003-2004  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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
// Author(s)     : Laurent Rineau

#ifndef QT_WIDGET_STYLE_EDITOR_H
#define QT_WIDGET_STYLE_EDITOR_H

#include <qframe.h>
#include "Qt_widget_styled_layer.h"

namespace CGAL {

class Qt_widget_style_editor : public QFrame {
  Q_OBJECT
public:

  typedef Qt_widget_styled_layer::Style Style;

  Qt_widget_style_editor(Style* style,
			 QWidget *parent = 0 , const char *name = 0);

  virtual ~Qt_widget_style_editor() {}

signals:
  void styleChanged();

private slots:
  void map(QColor);
  void map(int);
  void map(bool);

private:
  Style* style;
  QMap<const QObject*, QString> mapper;
}; // end of class Qt_widget_style_editor

} // end namespace CGAL

#endif // QT_WIDGET_STYLE_EDITOR_H
