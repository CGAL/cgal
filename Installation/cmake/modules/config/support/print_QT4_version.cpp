// Copyright (c) 2008 GeometryFactory, Sophia Antipolis (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial

// Tests if QT5 is available and prints its version string.

#include <iostream>
#include <QtCore/QtGlobal>

int main()
{
    std::cout << "version=" << QT_VERSION_STR << std::endl;

    return 0;
}
