// ============================================================================
//
// Copyright (c) 1997-2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// file          : include/CGAL/IO/Qt_widget_Triangulation_2.h
// package       : Qt_widget
// author(s)     : Radu Ursu <rursu@sophia.inria.fr>
// release       : 
// release_date  : 
//
// coordinator   : Laurent Rineau <rineau@clipper.ens.fr>
//
// ============================================================================


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
      if ((*fh).neighbor(i) > fh || t.is_infinite((*fh).neighbor(i)))
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
