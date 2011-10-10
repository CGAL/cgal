// Copyright (c) 2004,2005  INRIA Sophia-Antipolis (France).
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
// $URL$
// $Id$
//
//
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>

#ifndef SDG_DRAW_H
#define SDG_DRAW_H



template<class T, class Widget>
void draw_diagram(Widget& widget, const T& sdg)
{
  widget << CGAL::BLUE;
#if !defined (__POWERPC__)
  widget << CGAL::PointSize(3);
  widget << CGAL::LineWidth(3);
#endif
  sdg.draw_dual(widget);
}




#endif // SDG_DRAW_H
