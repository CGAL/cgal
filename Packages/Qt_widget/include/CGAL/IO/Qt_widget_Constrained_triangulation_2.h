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
// file          : include/CGAL/IO/Qt_widget_Constrained_triangulation_2.h
// package       : Qt_widget
// author(s)     : Radu Ursu
// release       : 
// release_date  : 
//
// coordinator   : Laurent Rineau <rineau@clipper.ens.fr>
//
// ============================================================================


#ifndef CGAL_QT_WIDGET_CONSTRAINED_TRIANGULATION_2_H
#define CGAL_QT_WIDGET_CONSTRAINED_TRIANGULATION_2_H

#include <CGAL/IO/Qt_widget.h>
#include <CGAL/Constrained_triangulation_2.h>

namespace CGAL{

template < class Gt, class Tds>
Qt_widget&
operator<<(Qt_widget& w,  const Constrained_triangulation_2<Gt,Tds> &t)
{
  w.lock();
  t.draw_triangulation(w);
  w.unlock();
  return w;
}

}//end namespace CGAL

#endif
