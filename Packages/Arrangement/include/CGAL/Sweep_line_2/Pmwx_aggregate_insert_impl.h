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
// release       : $CGAL_Revision: CGAL-2.3-I-81 $
// release_date  : $CGAL_Date: 2001/07/10 $
//
// file          : include/CGAL/Sweep_line_2/Sweep_curves_to_planar_map_2.h
// package       : Arrangement (2.07)
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

#ifndef CGAL_PMWX_AGGREGATE_INSERT_IMPL_H
#define CGAL_PMWX_AGGREGATE_INSERT_IMPL_H

#include <list>

#include <CGAL/Sweep_line_base_2.h>


CGAL_BEGIN_NAMESPACE


template <class CurveInputIterator, class SweepLineTraits_2, 
          class PM_, class Change_notification_>
class Pmwx_aggregate_insert_impl :
               public Sweep_line_base_2<CurveInputIterator, SweepLineTraits_2>
{
public:
  typedef SweepLineTraits_2 Traits;
  typedef typename Traits::X_curve_2 X_curve_2;

  typedef PM_ PM;
  typedef typename  PM::Halfedge_iterator  Halfedge_iterator; 
  
  typedef Change_notification_ Change_notification;
  typedef Sweep_line_base_2<CurveInputIterator, Traits> Base;

  Pmwx_aggregate_insert_impl() : 
    Base() {}
  
  Pmwx_aggregate_insert_impl(Traits *traits_) : 
    Base(traits_) {} 
  
  virtual ~Pmwx_aggregate_insert_impl() {}
  
  void Init(CurveInputIterator begin, CurveInputIterator end, PM &pm);

  void insert_curves(CurveInputIterator begin, 
                     CurveInputIterator end, 
                     PM &planarMap,
                     Change_notification* change_notification)
  {
    std::vector<X_curve_2> subcurves;
    Init(begin, end, planarMap);
    PerformIntersection(std::back_inserter(subcurves));
    planarMap.non_intersecting_insert(subcurves.begin(), subcurves.end());
  }
  
protected:



private:
};

/*! Initializes the data structures to work with:
  - x-monotonize the inf\put curves
  - for each end point of each curve create an event
  - for each curve in the planarmap do the same
  - initialize the event queue
  -
*/
template <class CurveInputIterator,  class SweepLineTraits_2,
          class PM_, class Change_notification_>
inline void 
Pmwx_aggregate_insert_impl<CurveInputIterator, SweepLineTraits_2, PM_, 
Change_notification_>::
Init(CurveInputIterator begin, CurveInputIterator end, PM &pm)
{
  Sweep_line_base_2<CurveInputIterator, Traits>::Init(begin, end);

  Halfedge_iterator eit;
  for (eit = pm.halfedges_begin(); eit != pm.halfedges_end(); ++eit, ++eit) 
  {
    InitCurve(eit->curve());
  }
  pm.clear();
}

CGAL_END_NAMESPACE

#endif
