// Copyright (c) 2006
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Sylvain Pion

// Tests if OPEN GL and GLU are available.

#if defined(_MSC_VER)
#include <wtypes.h>
#include <wingdi.h>
#endif

#ifdef __APPLE__
#  include <OpenGL/gl.h>
#  include <OpenGL/glu.h>
#else
#  include <GL/gl.h>
#  include <GL/glu.h>
#endif // __APPLE__
#include <iostream>

int main()
{
    std::cout << "version=" <<
#if   defined GL_VERSION_2_1
    "2.1"
#elif defined GL_VERSION_2_0
    "2.0"
#elif defined GL_VERSION_1_5
    "1.5"
#elif defined GL_VERSION_1_4
    "1.4"
#elif defined GL_VERSION_1_3
    "1.3"
#elif defined GL_VERSION_1_2
    "1.2"
#elif defined GL_VERSION_1_1
    "1.1"
#elif defined GL_VERSION_1_0
    "1.0"
#else
    "unknown"
#endif
    << std::endl;

    return 0;
}
