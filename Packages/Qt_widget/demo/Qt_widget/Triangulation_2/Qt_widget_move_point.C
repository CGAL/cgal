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
// file          : src/Qt_widget_MovePoint.C
// package       : Qt_widget
// author(s)     : Laurent Rineau
// release       : 
// release_date  : 
//
// coordinator   : Laurent Rineau <rineau@clipper.ens.fr>
//
// ============================================================================

#ifdef CGAL_USE_QT

#include "Qt_widget_move_point.h"

void Qt_widget_movepoint_helper::delete_vertex() { delete_vertexi(); };
void Qt_widget_movepoint_helper::move_vertex() { move_vertexi(); };
void Qt_widget_movepoint_helper::change_weight() { change_weighti(); };
void Qt_widget_movepoint_helper::stateChanged(int i){
  if(i==2)
    activate();
  else if(i == 0)
    deactivate();
}

#include "Qt_widget_move_point.moc"

#endif
