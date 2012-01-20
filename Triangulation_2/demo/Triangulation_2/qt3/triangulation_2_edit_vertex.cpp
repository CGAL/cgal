// Copyright (c) 1997-2002  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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
// Author(s)     : Laurent Rineau

#include <CGAL/basic.h>


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

