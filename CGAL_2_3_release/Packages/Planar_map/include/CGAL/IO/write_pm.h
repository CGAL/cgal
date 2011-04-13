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

#ifndef CGAL_PLANAR_MAP_2_H
#include <CGAL/Planar_map_2.h>
#endif

#ifndef CGAL_INVERSE_INDEX_H
#include <CGAL/Inverse_index.h>
#endif

#include <iostream>

CGAL_BEGIN_NAMESPACE


template <class PM, class Writer>
void write_pm( const PM& pm,
               Writer& writer, 
               std::ostream& o) 
{

  //typedef Planar_map_2<Dcel,Traits>                     PM;
  typedef typename PM::Halfedge_const_iterator          HCI;
  typedef typename PM::Vertex_const_iterator            VCI;
  //typedef Inverse_index<HCI>                            H_index;
  //typedef Inverse_index<VCI>                            V_index;
  
  //H_index h_index(pm.halfedges_begin(), pm.halfedges_end()); 
  //V_index v_index(pm.vertices_begin(), pm.vertices_end()); 
  
  // Print header. write #vertices, #halfedges, #faces.
  writer.write_title("Printing Planar map" /*, pm.number_of_vertices(), pm.number_of_halfedges(), pm.number_of_faces()*/);
  writer.write_comment("Printing number of vertices halfedges and faces in Planar map");
  writer.write_pm_vhf_sizes();

  writer.write_comment("vertices", pm.number_of_vertices());
    
  typename PM::Vertex_const_iterator vi;
  for (vi = pm.vertices_begin(); vi != pm.vertices_end(); ++vi) {
    writer.write_vertex(vi);
  }
  
  writer.write_comment("halfedges", pm.number_of_halfedges());
  //writer.write_halfedges_header();
  
  typename PM::Halfedge_const_iterator ei;
  for (ei = pm.halfedges_begin(); ei != pm.halfedges_end(); ++ei) {
    /* if (writer.verbose()) 
       writer.write_vertex(ei->target());
       else
       writer.write_index( v_index[ VCI(ei->target())] );*/
    
    writer.write_halfedge(ei);
  } 
  
  writer.write_comment("faces", pm.number_of_faces());
  //writer.write_faces_header();
  
  typename PM::Face_const_iterator fi;
  //typename Planar_map_2<Dcel,Traits>::Holes_const_iterator iccbit;
  //typename Planar_map_2<Dcel,Traits>::Ccb_halfedge_const_circulator ccb_circ;
  
  //int ck;
  int fc_counter;
  for (fc_counter=0,fi = pm.faces_begin(); fi != pm.faces_end(); ++fi,fc_counter++) { 
    writer.write_face(fi);
    
    /*writer.write_face_begin();
      
      if (fi->is_unbounded()){
      writer.write_comment("UNBOUNDED");
      }
      else {
      
      writer.write_comment("outer ccb");
      
      typename Planar_map_2<Dcel,Traits>::Ccb_halfedge_const_circulator first = fi->outer_ccb(), iter = first;
      
      std::size_t n = circulator_size(first);
      writer.write_number_halfedges_on_boundary(n);
      
      do{
      writer.write_index( h_index[HCI(iter)] );
      ++iter;
      } while (iter != first);
      
      writer.write_component_face_end();
      }
      
      for (ck=0, iccbit = (*fi).holes_begin(); iccbit != (*fi).holes_end(); ++iccbit, ck++);
      writer.write_comment("number of holes");
      writer.write_number_of_holes(ck);
      
      for (iccbit = (*fi).holes_begin(); iccbit != (*fi).holes_end(); ++iccbit){
      writer.write_comment("inner ccb");
      
      ccb_circ = (*iccbit);
      typename Planar_map_2<Dcel,Traits>::Ccb_halfedge_const_circulator 
      //typename Planar_map_2<Dcel,Traits>::Ccb_halfedge_circulator 
      iter = ccb_circ;
      
      std::size_t n = circulator_size(iter);
      writer.write_number_halfedges_on_boundary(n);
      
      do {
      writer.write_index( h_index[ HCI(iter)] );
      ++iter;
      } while (iter != ccb_circ);
      
      writer.write_face_end();
      }
      
      writer.write_component_face_end();*/
  }
  
  writer.write_title("End of Planar map");
  //writer.write_footer();
}


CGAL_END_NAMESPACE

#endif











