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
 
// use the port of the LEDA window (CGAL::window) or the "normal" LEDA window ... 
 
#if defined(CGAL_USE_CGAL_WINDOW)
#include <CGAL/IO/cgal_window.h>
#else
#include <CGAL/IO/leda_window.h>
#endif

#ifndef IO_TRIANGULATION_WINDOW_STREAM_H
#include <CGAL/IO/triangulation_Window_stream.h>
#endif  // IO_TRIANGULATION_WINDOW_STREAM_H
#ifndef IO_OPTIMISATION_WINDOW_STREAM_H
#include <CGAL/IO/optimisation_Window_stream.h>
#endif // IO_OPTIMISATION_WINDOW_STREAM_H
#ifndef IO_POLYGON_WINDOW_STREAM_H
#include <CGAL/IO/polygon_Window_stream.h>
#endif // IO_POLYGON_WINDOW_STREAM_H

