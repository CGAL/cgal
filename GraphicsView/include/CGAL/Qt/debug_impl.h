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
// SPDX-License-Identifier: GPL-3.0+
// 
//
// Author(s)     : Andreas Fabri <Andreas.Fabri@geometryfactory.com>
//                 Laurent Rineau <Laurent.Rineau@geometryfactory.com>
   
#ifdef CGAL_HEADER_ONLY
#define CGAL_INLINE_FUNCTION inline

#include <CGAL/license/GraphicsView.h>

#else
#define CGAL_INLINE_FUNCTION
#endif

#include <CGAL/Qt/debug.h>
#include <QDir>
#include <iostream>
#include <QtOpenGL/qgl.h>
#include <qopenglfunctions.h>
namespace CGAL {
namespace Qt {


CGAL_INLINE_FUNCTION
void traverse_resources(const QString& name, const QString& dirname, int indent)
{
  std::cerr << qPrintable(QString(indent, ' '))
            << qPrintable(name);
  QString fullname = 
    dirname.isEmpty() ?
    name :
    dirname + "/" + name;
  QDir dir(fullname);
  if(dir.exists()) {
    std::cerr << "/\n";
    Q_FOREACH(QString path, dir.entryList())
    {
      traverse_resources(path, fullname, indent + 2);
    }
  }
  else {
    std::cerr << "\n";
  }
}

} // namesapce Qt
} // namespace CGAL
