// Copyright (c) 1999  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Stefan Schirra
 
#ifndef CGAL_IO_FORWARD_DECL_WINDOW_STREAM_H
#define CGAL_IO_FORWARD_DECL_WINDOW_STREAM_H

#if defined(CGAL_USE_CGAL_WINDOW)
CGAL_BEGIN_NAMESPACE
class window;
CGAL_END_NAMESPACE
#else
#include <LEDA/window.h>
class leda_window;
#endif

CGAL_BEGIN_NAMESPACE

#if defined(CGAL_USE_CGAL_WINDOW)
typedef   CGAL::window    Window_stream;
#else
typedef   leda_window    Window_stream;
#endif

CGAL_END_NAMESPACE

#endif // CGAL_IO_FORWARD_DECL_WINDOW_STREAM_H
