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
// release       : $CGAL_Revision: CGAL-2.5-I-11 $
// release_date  : $CGAL_Date: 2002/08/04 $
//
// file          : include/CGAL/Map_overlay_incremental.h
// package       : Map_overlay (1.12)
// maintainer    : Efi Fogel <efif@math.tau.ac.il>
// source        : 
// revision      : 
// revision_date : 
// author(s)     : Eti Ezra          <estere@post.tau.ac.il>
//
// coordinator   : Tel-Aviv University (Dan Halperin <halperin@math.tau.ac.il>)
//
// Chapter       : 
// ======================================================================
#ifndef CGAL_MAP_OVERLAY_INCREMENTAL_H
#define CGAL_MAP_OVERLAY_INCREMENTAL_H

#include <vector>
#include <list>

#ifndef CGAL_MAP_OVERLAY_BASE_H
#include <CGAL/Map_overlay_base.h>
#endif

#ifndef CGAL_MAP_OVERLAY_DEFAULT_NOTIFIER_H
#include <CGAL/Map_overlay_default_notifier.h>
#endif

//#include <CGAL/IO/leda_window.h>  //used for visualization -
//should be placed after the other CGAL files.

CGAL_BEGIN_NAMESPACE

template  <class Arrangement_, 
  class Map_overlay_ChangeNotification_ = 
        Map_overlay_default_notifier<Arrangement_> > 
class  Map_overlay_incremental : 
  public Map_overlay_base<Arrangement_, Map_overlay_ChangeNotification_> {
public:  
  
  typedef Arrangement_                          Arrangement;
  typedef Map_overlay_ChangeNotification_       Map_overlay_ChangeNotification;
  //typedef typename Arrangement::Planar_map                      PM;
  //typedef typename Arrangement::Pmwx                          Pmwx;
  typedef typename Arrangement::Vertex            Vertex;
  typedef typename Arrangement::Face              Face;
  typedef typename Arrangement::Halfedge          Halfedge;
  typedef typename Arrangement::Vertex_handle     Vertex_handle;
  typedef typename Arrangement::Halfedge_handle          Halfedge_handle;
  typedef typename Arrangement::Face_handle              Face_handle;
  typedef typename Arrangement::Vertex_const_handle      Vertex_const_handle;
  typedef typename Arrangement::Halfedge_const_handle    Halfedge_const_handle;
  typedef typename Arrangement::Face_const_handle        Face_const_handle;
  typedef typename Arrangement::Vertex_iterator          Vertex_iterator;
  typedef typename Arrangement::Halfedge_iterator        Halfedge_iterator;
  typedef typename Arrangement::Face_iterator            Face_iterator;
  typedef typename Arrangement::Vertex_const_iterator    Vertex_const_iterator;
  typedef typename Arrangement::Halfedge_const_iterator   
                                                       Halfedge_const_iterator;
  typedef typename Arrangement::Face_const_iterator    Face_const_iterator;
  typedef typename Arrangement::Ccb_halfedge_circulator  
                                                    Ccb_halfedge_circulator;
  typedef typename Arrangement::Ccb_halfedge_const_circulator   
                                               Ccb_halfedge_const_circulator;
  typedef typename Arrangement::Holes_iterator         Holes_iterator;
  typedef typename Arrangement::Holes_const_iterator   Holes_const_iterator;
  typedef typename Arrangement::Locate_type            Locate_type;

  typedef typename  Arrangement::Traits               Traits;
  typedef typename  Traits::X_curve                   X_curve;
  typedef typename  Traits::Point                     Point;
  
private:
  typedef Map_overlay_incremental<Arrangement,Map_overlay_ChangeNotification> 
                                                                         Self;
public:

  Map_overlay_incremental() {}

  ~Map_overlay_incremental() {}
  
  void  map_overlay(const Arrangement &a1, 
                    const Arrangement &a2, 
                    Map_overlay_ChangeNotification *change_notf, 
                    Arrangement &result)
  {
    Halfedge_const_iterator  halfedge_iter;
    
    // Creaing the new arrangement induced by both subdivisions. 
    // inserting the curves of every halfedge one by one.
    
    for (halfedge_iter = a1.halfedges_begin(); 
         halfedge_iter != a1.halfedges_end(); 
         ++halfedge_iter, ++halfedge_iter){
      change_notf->set_curve_attributes(halfedge_iter->curve(), 
                                        halfedge_iter, 
                                        true);
      result.insert(halfedge_iter->curve(), change_notf);
    }
    
    for (halfedge_iter = a2.halfedges_begin(); 
         halfedge_iter != a2.halfedges_end(); 
         ++halfedge_iter, ++halfedge_iter){
      change_notf->set_curve_attributes(halfedge_iter->curve(), 
                                        halfedge_iter, 
                                        false); 
      result.insert(halfedge_iter->curve(), change_notf);
    }
    
    //cout<<"After creating Arr\n";
    change_notf->update_all_faces(result);
  }

};

CGAL_END_NAMESPACE

#endif











