// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-2.3-I-44 $
// release_date  : $CGAL_Date: 2001/03/09 $
//
// file          : include/CGAL/IO/write_pm.h
// package       : pm (5.45)
// maintainer    : Eyal Flato <flato@math.tau.ac.il>
// source        : 
// revision      : 
// revision_date : 
// author(s)     : Eti Ezra <estere@post.tau.ac.il>
//
//
// coordinator   : Tel-Aviv University (Dan Halperin <halperin@math.tau.ac.il>)
//
// Chapter       : 
// ======================================================================

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
