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
// file          : triangulation_2_edit_vertex.C
// package       : Qt_widget
// author(s)     : Laurent Rineau
// release       : 
// release_date  : 
//
// coordinator   : Laurent Rineau <rineau@clipper.ens.fr>
//
// ============================================================================

#ifdef CGAL_USE_QT

#include "triangulation_2_edit_vertex.h"

void triangulation_2_edit_vertex_helper::delete_vertex()
{ 
  delete_vertexi();
  emit(triangulation_changed());
};
void triangulation_2_edit_vertex_helper::move_vertex() { move_vertexi(); };
void triangulation_2_edit_vertex_helper::change_weight() { change_weighti(); };
void triangulation_2_edit_vertex_helper::stateChanged(int i){
  if(i==2)
    activate();
  else if(i == 0)
    deactivate();
}

#include "triangulation_2_edit_vertex.moc"

#endif
