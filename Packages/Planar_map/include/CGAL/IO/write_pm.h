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
  //typedef Planar_map_2<Dcel,Traits>                     PM;
  typedef typename PM::Halfedge_const_iterator          HCI;
  typedef typename PM::Vertex_const_iterator            VCI;
  typedef typename PM::Face_const_iterator              FCI;
  typedef Inverse_index<HCI>                            H_index;
  typedef Inverse_index<VCI>                            V_index;
  
  //H_index h_index(pm.halfedges_begin(), pm.halfedges_end()); 
  //V_index v_index(pm.vertices_begin(), pm.vertices_end()); 
  
  // Print header. write #vertices, #halfedges, #faces.
  writer.write_title("Printing Planar map" /*, pm.number_of_vertices(), pm.number_of_halfedges(), pm.number_of_faces()*/);
  writer.write_comment("Printing number of vertices halfedges and faces in Planar map");
  writer.write_pm_vhf_sizes(pm.number_of_vertices(),
                            pm.number_of_halfedges(),
                            pm.number_of_faces());

  writer.write_comment("vertices", pm.number_of_vertices());
  //writer.write_vertices(pm.vertices_begin(), pm.vertices_end());
  VCI v_iter;
  for (v_iter = pm.vertices_begin(); v_iter != pm.vertices_end(); ++v_iter)
    writer.write_vertex(v_iter);
  
  writer.write_comment("halfedges", pm.number_of_halfedges());
  //writer.write_halfedges_header();
  
  // writer.write_index( v_index[ VCI(ei->target())] );*/
  // writer.write_halfedges(pm.halfedges_begin(), pm.halfedges_end());
  V_index v_index(pm.vertices_begin(), pm.vertices_end());
  HCI h_iter;
  for (h_iter = pm.halfedges_begin(); h_iter != pm.halfedges_end(); ++h_iter)
    writer.write_halfedge(h_iter, v_index);

  writer.write_comment("faces", pm.number_of_faces());
  //writer.write_faces_header();
  
  //typename Planar_map_2<Dcel,Traits>::Holes_const_iterator iccbit;
  //typename Planar_map_2<Dcel,Traits>::Ccb_halfedge_const_circulator ccb_circ;
  
  //int ck;
  // writer.write_faces(pm.faces_begin(), pm.faces_end());
  H_index h_index(pm.halfedges_begin(), pm.halfedges_end());
  FCI f_iter;
  for (f_iter = pm.faces_begin(); f_iter != pm.faces_end(); ++f_iter)
    writer.write_face(f_iter, h_index);
  
  writer.write_title("End of Planar map");
  //writer.write_footer();
}


CGAL_END_NAMESPACE

#endif
