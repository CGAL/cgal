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

#ifndef QT_WIDGET_STYLED_LAYER_H
#define QT_WIDGET_STYLED_LAYER_H

#include <qvariant.h>
#include <qstring.h>
#include <qmap.h>

#include <CGAL/IO/Qt_widget_layer.h>
#include <CGAL/function_objects.h>

namespace CGAL {

class Qt_widget_style : public QObject {
  Q_OBJECT

  typedef QMap<QString,QVariant> Map;
public:

  typedef Map::const_iterator const_iterator;
  typedef Map::size_type size_type;

  Qt_widget_style()
    : map() {};

public slots:
  void setBool(QString name, bool b);
  void setInt(QString name, int i);
  void setColor(QString name, QColor c);
  void setPointStyle(QString name, PointStyle s);

public:
  bool getBool(QString name);
  int getInt(QString name);
  QColor getColor(QString name);
  PointStyle getPointStyle(QString name);

  const_iterator begin() const
  {
    return map.begin();
  }

  const_iterator end() const
  {
    return map.end();
  }

  size_type size() const
  {
    return map.size();
  }

private:
  Map map;
};

class Qt_widget_styled_layer : public Qt_widget_layer {
  Q_OBJECT
public:
  typedef Qt_widget_style Style;

  Qt_widget_styled_layer(Style* style = 0,
			 QObject * parent=0, const char * name=0);

  ~Qt_widget_styled_layer();

  void setStyle(Style* style);
  Style * style()  { return style_; }
private:
  Style* style_;
  bool private_style;
};

} // namespace CGAL

#endif // QT_WIDGET_STYLED_LAYER_H
