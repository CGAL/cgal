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
// Author(s) : Saar Katz <kats.saar@gmail.com>

#ifndef CGAL_POLYGONLISTITEM_H
#define CGAL_POLYGONLISTITEM_H

#include <QColor>

#include "Typedefs.h"
#include "PolygonItem.h"

class PolygonListItem {
public:
  PolygonListItem(QString name);
  PolygonListItem(QString name, QColor color, bool visibility, 
                  PolygonItem* polygonItem);

  QString getName();
  QColor getColor();
  bool isVisible();
  PolygonItem* getPolygonItem();

  void setName(QString name);
  void setColor(QColor color);
  void setVisible(bool visibility);
  void setPolygonItem(PolygonItem* polygonItem);

private:
	QString m_name;
	QColor m_color;
	bool m_visibility;
	PolygonItem* m_polygonItem;
};

#endif // CGAL_POLYGONLISTITEM_H