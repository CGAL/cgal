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
// Author(s)     : Tali Zvi <talizvi@post.tau.ac.il>
//                 Baruch Zukerman <baruchzu@post.tau.ac.il>

#ifndef CGAL_PMWX_SWEEP_LINE_EVENT_H
#define CGAL_PMWX_SWEEP_LINE_EVENT_H

#include <CGAL/Sweep_line_2/Sweep_line_event.h>
#include <CGAL/Sweep_line_2/Pmwx_insert_info.h>
#include <CGAL/assertions.h>

CGAL_BEGIN_NAMESPACE

/*! @class Pmwx_sweep_line_event
 *
 * Stores the data associated with an event.
 * In addition to the information stored in Sweep_line_event, when 
 * constructing a * planar map, additional information is kept, in 
 * order to speed insertion of curves into the planar map.
 *
 * The additional infomation contains the following:
 * - among the left curves of the event, we keep the highest halfedge that 
 *   was inserted into the planar map at any given time.
 *
 * Inherits from Sweep_line_event.
 * \sa Sweep_line_event
 */

template<class SweepLineTraits_2, class CurveWrap>
class Pmwx_sweep_line_event : 
  public Sweep_line_event<SweepLineTraits_2, CurveWrap>
{
public:
  typedef SweepLineTraits_2 Traits;
  typedef typename Traits::X_monotone_curve_2 X_monotone_curve_2;
  typedef typename Traits::Point_2 Point_2;

  typedef Sweep_line_event<SweepLineTraits_2, CurveWrap> Base;
  typedef Pmwx_sweep_line_event<Traits, CurveWrap> Self;

  typedef CurveWrap SubCurve;
  typedef std::list<SubCurve *> SubcurveContainer;
  typedef typename SubcurveContainer::iterator SubCurveIter;
  typedef typename SubcurveContainer::reverse_iterator SubCurveRevIter;

    
  typedef typename SubCurve::PmwxInsertInfo        PmwxInsertInfo;
  typedef typename PmwxInsertInfo::Halfedge_handle Halfedge_handle;

  using Base::m_rightCurves;


  Pmwx_sweep_line_event()
  {}


  /*! destructor */
  ~Pmwx_sweep_line_event()
  {}


  // TODO - remove ?
  void init(const Point_2 &point)
  {
    Base::init(point);
  }

  PmwxInsertInfo *get_insert_info()
  {
    return &m_insertInfo;
  }

  std::pair<bool, SubCurveIter> add_curve_to_right(SubCurve *curve)
  {
    std::pair<bool,SubCurveIter> res = Base::add_curve_to_right(curve);
    if(res.second != m_rightCurves.end() && res.first == false )
      m_insertInfo.inc_right_curves_counter();
    return res;
  }


  int get_halfedge_jump_count(CurveWrap *curve)
  {
    int i = 0;
    int skip = 0;
   
    SubCurveIter iter = m_rightCurves.begin();
    for(; iter!=m_rightCurves.end(); ++iter)
    {
      if((*iter) == NULL)
        skip++;
    }
    skip--;  // now 'skip' holds the amount of the right curves of the event
		         // that are already inserted to the planar map  - 1 (minus 1)

    iter = m_rightCurves.end();
    --iter;
     
    unsigned int num_left_curves = this->get_num_left_curves();
    for ( ; iter != m_rightCurves.begin() ; --iter )
    {
      if(curve == (*iter))
      {
        (*iter) = NULL; 
        if (( i == 0 ) && ( num_left_curves == 0 )) 
        {
          return skip;
        }
        if ( num_left_curves == 0 ) 
	      {   
          return i-1;
        }
        return i;
      }
      if ( (*iter) == NULL )
        i++;
    }

    CGAL_assertion(curve == (*iter));
    (*iter) = NULL; 
    
    if ( num_left_curves == 0 )
      i--;
    return i;
  }

 

  bool is_curve_largest(CurveWrap *curve)
  {
    for( SubCurveRevIter rev_iter = m_rightCurves.rbegin();
         rev_iter != m_rightCurves.rend() && curve != (*rev_iter) ;
         ++rev_iter)
    {
      if((*rev_iter) == NULL)
         return false;
    }
    return true;
  }
  

protected:
  PmwxInsertInfo m_insertInfo;
};

CGAL_END_NAMESPACE

#endif // CGAL_PMWX_SWEEP_LINE_EVENT_H
