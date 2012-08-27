
// Copyright (c) 2008  GeometryFactory Sarl (France).
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
// $URL$
// $Id$
// 
//
// Author(s)     : Andreas Fabri <Andreas.Fabri@geometryfactory.com>
//                 Laurent Rineau <Laurent.Rineau@geometryfactory.com>

#ifndef CGAL_QT_GRAPHICS_VIEW_INPUT_H
#define CGAL_QT_GRAPHICS_VIEW_INPUT_H
#include <CGAL/export/Qt4.h>
#include <CGAL/auto_link/Qt4.h>
#ifndef Q_MOC_RUN
#  include <CGAL/Object.h>
#endif
#include <QObject>

namespace CGAL {
namespace Qt {
class CGAL_QT4_EXPORT GraphicsViewInput  : public QObject
{
  Q_OBJECT

public:
  GraphicsViewInput(QObject* parent) 
    : QObject(parent)
  {}

signals:
  void generate(CGAL::Object o);
  void modelChanged();

public slots:

  virtual void processInput(CGAL::Object /*o*/) {}

};

} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_GRAPHICS_VIEW_INPUT_H
