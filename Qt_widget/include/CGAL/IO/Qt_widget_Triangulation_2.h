// Copyright (c) 2002-2004  INRIA Sophia-Antipolis (France).
// All rights reserved.
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
// 
//
// Author(s)     : Radu Ursu

#ifndef CGAL_QT_WIDGET_TRIANGULATION_2_H
#define CGAL_QT_WIDGET_TRIANGULATION_2_H

#include <CGAL/IO/Qt_widget.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/apply_to_range.h>


namespace CGAL {

template <class Gt, class Tds>
class Draw_triangulation {
private:
  const Triangulation_2<Gt, Tds>& t;
  Qt_widget& w;
public:
  Draw_triangulation(const Triangulation_2<Gt, Tds>& _t, Qt_widget& _w)
    : t(_t), w(_w)
  {}
  void operator()(typename Triangulation_2<Gt, Tds>::Face_handle fh)
  {
    for (int i=0; i<3; i++)
      if (fh < fh->neighbor(i) || t.is_infinite(fh->neighbor(i)))
        w << t.segment(fh,i);
  }
};


template < class Gt, class Tds>
Qt_widget&
operator<<(Qt_widget& w,  const Triangulation_2<Gt, Tds> &t)
{
  if (t.dimension()<2) {
    t.draw_triangulation(w);
    return w;
  }
  typedef typename Triangulation_2<Gt, Tds>::Point OpPoint;
  w.lock();
  Draw_triangulation<Gt, Tds> draw(t, w);
  apply_to_range(t, OpPoint(w.x_min(), w.y_max()), 
                 OpPoint(w.x_max(), w.y_min()), draw);
  w.unlock();
  return w;
}


}// end namespace CGAL

#endif
