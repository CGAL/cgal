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
// file          : include/CGAL/Pmwx_sweep_line_curve.h
// package       : arr (1.87)
// maintainer    : Tali Zvi <talizvi@post.tau.ac.il>
// source        : 
// revision      : 
// revision_date : 
// author(s)     : Tali Zvi <talizvi@post.tau.ac.il>
//
//
// coordinator   : Tel-Aviv University (Dan Halperin <halperin@math.tau.ac.il>)
//
// Chapter       : 
// ======================================================================
#ifndef CGAL_PMWX_SWEEP_LINE_CURVE_H
#define CGAL_PMWX_SWEEP_LINE_CURVE_H

#include <CGAL/Sweep_line_2/Sweep_line_subcurve.h>
#include <CGAL/Sweep_line_2/Pmwx_insert_info.h>
#include <CGAL/Sweep_line_2/Pmwx_sweep_line_event.h>

CGAL_BEGIN_NAMESPACE

template<class SweepLineTraits_2, class VertexHandle, class HalfedgeHandle>
class Pmwx_sweep_line_curve : public Sweep_line_subcurve<SweepLineTraits_2>
{
public:
  typedef SweepLineTraits_2 Traits;
  typedef typename Traits::Point_2 Point_2;
  typedef typename Traits::Curve_2 Curve_2;
  typedef typename Traits::X_curve_2 X_curve_2;
  typedef Sweep_line_subcurve<SweepLineTraits_2> Base;
  typedef Pmwx_sweep_line_curve<Traits, VertexHandle, HalfedgeHandle> Self;

  typedef Pmwx_insert_info<VertexHandle,HalfedgeHandle> PmwxInsertInfo;
  typedef Pmwx_sweep_line_event<Traits, Self> Event;

  Pmwx_sweep_line_curve(int id, X_curve_2 &curve, Point_2 *reference, 
			SweepLineTraits_2 *traits) : 
    Base(id, curve, reference, traits) , m_insertInfo(0), m_lastEvent(0)
  {
  }

  void setInsertInfo(PmwxInsertInfo *insertInfo) {
    m_insertInfo = insertInfo;
  }

  PmwxInsertInfo *getInsertInfo() const {
    return m_insertInfo;
  }

  void setLastEvent(Event *e) {
    m_lastEvent = e;
  }

  Event *getLastEvent() const {
    return m_lastEvent;
  }

private:

  PmwxInsertInfo *m_insertInfo;
  Event *m_lastEvent;
  
};

CGAL_END_NAMESPACE

#endif // CGAL_PMWX_SWEEP_LINE_CURVE_H

