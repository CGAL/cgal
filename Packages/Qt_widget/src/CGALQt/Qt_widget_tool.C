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
// file          : src/Qt_widget_tool.C
// package       : Qt_widget
// author(s)     : Laurent Rineau & Radu Ursu
// release       : 
// release_date  : 
//
// coordinator   : Laurent Rineau <rineau@clipper.ens.fr>
//
// ============================================================================

#ifdef CGAL_USE_QT

#include <CGAL/Bbox_2.h>
#include <CGAL/IO/Qt_widget_tool.h>
#include <CGAL/IO/Qt_widget.h>

namespace CGAL {

Qt_widget_tool::Qt_widget_tool() : widget(0) {};

void Qt_widget_tool::attach(Qt_widget *w)
{
  widget=w;
  attaching();
}

void Qt_widget_tool::detach()
{
  detaching();
  widget=0;
}


} // namespace CGAL

#include "Qt_widget_tool.moc"

#endif
