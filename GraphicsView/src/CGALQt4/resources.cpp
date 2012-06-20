// Copyright (c) 2011  GeometryFactory Sarl (France).
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
// Author(s)     : Laurent Rineau <Laurent.Rineau@geometryfactory.com>

#include <QDir>
#include <CGAL/Qt/resources.h>

// cannot use namespaces because of the Q_INIT_RESOURCE macro
void CGAL_Qt4_init_resources() {
  Q_INIT_RESOURCE(File);
  Q_INIT_RESOURCE(Triangulation_2); 
  Q_INIT_RESOURCE(Input);
  Q_INIT_RESOURCE(CGAL);
}
