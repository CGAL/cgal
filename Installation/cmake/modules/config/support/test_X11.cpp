// Copyright (c) 2005
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
//
// Author(s)     : Michael Hoffmann

// Tests if X11 is available.

#include <iostream>
#include <X11/Xlib.h>

int main()
{
  Display* x = 0;
  std::cout << "version=" << X_PROTOCOL << "." << X_PROTOCOL_REVISION
            << std::endl;
  return 0;
}
