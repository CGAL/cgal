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
// file          : include/CGAL/Pmwx_sweep_line_event.h
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
#ifndef CGAL_PMWX_SWEEP_LINE_EVENT_H
#define CGAL_PMWX_SWEEP_LINE_EVENT_H

#include <CGAL/Sweep_line_2/Sweep_line_event.h>
//#include <CGAL/Sweep_line_2/Pmwx_sweep_line_curve.h>
#include <CGAL/Sweep_line_2/Pmwx_insert_info.h>

CGAL_BEGIN_NAMESPACE

template<class SweepLineTraits_2, class CurveWrap>
class Pmwx_sweep_line_event : 
  public Sweep_line_event<SweepLineTraits_2, CurveWrap>
{
public:
  typedef SweepLineTraits_2 Traits;
  typedef typename Traits::X_curve_2 X_curve_2;
  typedef typename Traits::Point_2 Point_2;

  typedef Sweep_line_event<SweepLineTraits_2, CurveWrap> Base;
  typedef Pmwx_sweep_line_event<Traits, CurveWrap> Self;

  typedef typename CurveWrap::PmwxInsertInfo PmwxInsertInfo;

  typedef std::list<Self *> VerticalXEventList;
  typedef VerticalXEventList::iterator VerticalXEventListIter; 

  /*! Constructor */
  Pmwx_sweep_line_event(const Point_2 &point, Traits *traits) :
    Base(point, traits)
    {}

  PmwxInsertInfo *getInsertInfo() {
    return &m_insertInfo;
  }

  /*! Insert a new intersection point on any of the vertical curves.
      The list of points is sorted by their y values.
      If the requireSort flag is true, the appripriate place in the list 
      is searched for. If not, the point is assumed to have the largest y 
      value, and is inserted at the end of the list. 
      If the pioint already exists, the point is nott inserted again.
      @param p a reference to the point
      @param requireSort false if the point is to be added at the end
      of the list.
  */
  void addVerticalCurveXEvent(Self *e, bool requireSort=false) 
  {
    if ( m_verticalCurveXEvents.empty() ) 
    {
      m_verticalCurveXEvents.push_back(e); 
      return;
    }
    
    if ( !requireSort ) 
    {
      if ( m_verticalCurveXEvents.back() != e ) {
	m_verticalCurveXEvents.push_back(e);
      }
    } else
    {
      VerticalXEventListIter iter = m_verticalCurveXEvents.begin();
      while ( iter != m_verticalCurveXEvents.end() )
      {
	if ( m_traits->compare_y((*iter)->getPoint(), e->getPoint()) 
	     == SMALLER ) {
	  ++iter; 
	}
	else
	  break;
      }
      if ( iter == m_verticalCurveXEvents.end() )
	m_verticalCurveXEvents.push_back(e);
      else if (m_verticalCurveXEvents.back() != e) {
	m_verticalCurveXEvents.insert(iter, e);
      }
    }
  }

  VerticalXEventList &getVerticalXEventList() {
    return m_verticalCurveXEvents;
  }
  
private:
  PmwxInsertInfo m_insertInfo;

  VerticalXEventList m_verticalCurveXEvents;
};

CGAL_END_NAMESPACE

#endif // CGAL_PMWX_SWEEP_LINE_EVENT_H
