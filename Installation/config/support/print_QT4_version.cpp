// Copyright (c) 2008 GeometryFactory, Sophia Antipolis (France)
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

// Tests if QT5 is available and prints its version string.

#include <iostream>
#include <QtCore/QtGlobal>

int main()
{
    std::cout << "version=" << QT_VERSION_STR << std::endl;

    return 0;
}
