// Copyright (c) 1997-2000  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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
// Author(s)     : Michael Seel <seel@mpi-sb.mpg.de>
#ifndef FILTERED_EXTENDED_HOMOGENEOUS_WINDOW_STREAM_H
#define FILTERED_EXTENDED_HOMOGENEOUS_WINDOW_STREAM_H

#ifdef CGAL_USE_LEDA
#include <CGAL/Filtered_extended_homogeneous.h>
#include <CGAL/IO/Window_stream.h>

CGAL_BEGIN_NAMESPACE

template <class RT> 
CGAL::Window_stream& operator<<(CGAL::Window_stream& w, 
                                const Extended_point<RT>& p)
{ w.draw_filled_node(CGAL::to_double(p.x()),CGAL::to_double(p.y())); 
  return w;
}

template <class RT> 
CGAL::Window_stream& operator<<(CGAL::Window_stream& w, 
                                const Extended_segment<RT>& s)
{ w.draw_segment(CGAL::to_double(s.source().x()),
                 CGAL::to_double(s.source().y()),
                 CGAL::to_double(s.target().x()),
                 CGAL::to_double(s.target().y())); 
  return w;
}

template <class RT> 
leda_point pnt(const Extended_point<RT>& p)
{ return leda_point(CGAL::to_double(p.x()),CGAL::to_double(p.y())); }

CGAL_END_NAMESPACE
#endif // CGAL_USE_LEDA
#endif // FILTERED_EXTENDED_HOMOGENEOUS_WINDOW_STREAM_H

