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
// file          : include/CGAL/Holes_split_notifier.h
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
#ifndef CGAL_HOLES_SPLIT_NOTIFIER_H
#define CGAL_HOLES_SPLIT_NOTIFIER_H   


CGAL_BEGIN_NAMESPACE         

template <class Planar_map_>
class Holes_split_notifier : public Planar_map_::Change_notification
{
public:
  typedef Planar_map_                                   Planar_map;
  typedef typename Planar_map::Halfedge_handle          Halfedge_handle;
  typedef typename Planar_map::Halfedge_const_handle    Halfedge_const_handle;
  
  typedef typename Planar_map::Traits                   Traits;
  typedef typename Traits::Point                        Point;
  typedef typename Traits::X_curve                      X_curve;
  
  Holes_split_notifier() : decomposing_(false) {}
  
  virtual ~Holes_split_notifier() {}
  
  // We assume the curves comes from the original triangles are
  // oriented counter clockwise.
  void add_edge(const X_curve& cv, 
                Halfedge_handle e, 
                bool left_to_right, 
                bool overlap = false)
  { 
    if (decomposing_){
      e->set_decomposing(true);
      e->twin()->set_decomposing(true);
    }
  }
 
  
  //---------------------------------------- new functions  
  
  void  set_decomposing_edge(bool decomposing) 
  { 
    decomposing_ = decomposing; 
  }

  bool  is_decomposing_edge() const { 
    return decomposing_; 
  }
  
private:
  bool  decomposing_;
};

CGAL_END_NAMESPACE

#endif




