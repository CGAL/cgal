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

#include <CGAL/basic.h>


#include "Qt_widget_styled_layer.h"

namespace CGAL {

  Qt_widget_styled_layer::Qt_widget_styled_layer(
	    Qt_widget_styled_layer::Style* style,
	    QObject * parent,
	    const char * name)
    : Qt_widget_layer(parent, name),
      style_(style), private_style(false)
  {
    if( style == 0 )
      {
	this->style_ = new Style();
	private_style = true;
      }
  }

  Qt_widget_styled_layer::~Qt_widget_styled_layer()
  {
    if( private_style )
      delete style_;
  }

  void Qt_widget_styled_layer::setStyle(Qt_widget_styled_layer::Style* s)
  {
    if( private_style )
      delete style_;
    private_style = false;
    style_ = s;
  }

  void Qt_widget_style::setBool(QString name, bool b)
  {
    map[name] = b;
  }

  void Qt_widget_style::setInt(QString name, int i)
  {
    map[name] = i;
  }

  void Qt_widget_style::setColor(QString name, QColor c)
  {
    map[name] = c;
  }

  void Qt_widget_style::setPointStyle(QString name, PointStyle s)
  {
    map[name] = static_cast<uint>(s);
    map[name].cast(QVariant::UInt);
  }

  bool Qt_widget_style::getBool(QString name)
  {
    if( ! map.contains(name) )
      return false;
    else
      {
	CGAL_assertion( map[name].type() == QVariant::Bool );
	return map[name].asBool();
      }
  }

  int Qt_widget_style::getInt(QString name)
  {
    if( ! map.contains(name) )
      return 0;
    else
      {
	CGAL_assertion( map[name].type() == QVariant::Int );
	return map[name].asInt();
      }
  }

  QColor Qt_widget_style::getColor(QString name)
  {
    if( ! map.contains(name) )
      return QColor();
    else
      {
	CGAL_assertion( map[name].type() == QVariant::Color );
	return map[name].asColor();
      }
  }

  ::CGAL::PointStyle Qt_widget_style::getPointStyle(QString name)
  {
    if( ! map.contains(name) )
      return PointStyle();
    else
      {
	CGAL_assertion( map[name].type() == QVariant::UInt );
	return PointStyle(map[name].asUInt());
      }
  }

} // namespace CGAL

// moc_source_file: Qt_widget_styled_layer.h
#include "Qt_widget_styled_layer.moc"

