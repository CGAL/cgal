// Copyright (c) 1997  Tel-Aviv University (Israel).
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
// Author(s)     : Eti Ezra <estere@post.tau.ac.il>

#ifndef CGAL_IO_WRITE_PM_H
#define CGAL_IO_WRITE_PM_H 

#ifndef CGAL_BASIC_H
#include <CGAL/basic.h>
#endif

#include <iostream>

#ifndef CGAL_PLANAR_MAP_2_H
#include <CGAL/Planar_map_2.h>
#endif

#ifndef CGAL_INVERSE_INDEX_H
#include <CGAL/Inverse_index.h>
#endif

CGAL_BEGIN_NAMESPACE

template <class PM, class Writer>
void write_pm(const PM & pm, Writer & writer, std::ostream &)
{
  // Print header. write #vertices, #halfedges, #faces.
  writer.write_title("Begin Planar Map");
  writer.write_comment("Number of vertices halfedges and faces in Planar map");
  writer.write_pm_vhf_sizes(pm.number_of_vertices(),
                            pm.number_of_halfedges(),
                            pm.number_of_faces());

  writer.write_comment("vertices", pm.number_of_vertices());
  writer.write_vertices(pm.vertices_begin(), pm.vertices_end());
  
  writer.write_comment("halfedges", pm.number_of_halfedges());
  writer.write_halfedges(pm.halfedges_begin(), pm.halfedges_end());

  writer.write_comment("faces", pm.number_of_faces());
  writer.write_faces(pm.faces_begin(), pm.faces_end());
  
  writer.write_title("End Planar Map");
  //writer.write_footer();
}

CGAL_END_NAMESPACE

#endif
