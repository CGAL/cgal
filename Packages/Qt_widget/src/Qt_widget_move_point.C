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
// file          : src/Qt_Window_MovePoint.C
// package       : QT_window
// author(s)     : Laurent Rineau
// release       : 
// release_date  : 
//
// coordinator   : Laurent Rineau <rineau@clipper.ens.fr>
//
// ============================================================================

#ifdef CGAL_USE_QT

#include <CGAL/IO/Qt_Widget_MovePoint.h>

namespace CGAL {

void Qt_widget_movepoint_helper::delete_point() { delete_pointi(); };
void Qt_widget_movepoint_helper::move_point() { move_pointi(); };

} // namespace CGAL

#include "Qt_Widget_MovePoint.moc"

#endif
