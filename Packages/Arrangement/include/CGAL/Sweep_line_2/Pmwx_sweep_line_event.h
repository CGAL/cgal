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
#ifndef CGAL_PMWX_SWEEP_LINE_EVENT_H
#define CGAL_PMWX_SWEEP_LINE_EVENT_H

#include <CGAL/Sweep_line_2/Sweep_line_event.h>
#include <CGAL/Sweep_line_2/Pmwx_insert_info.h>
#include <CGAL/Sweep_line_2/Pmwx_sweep_line_functors.h>
#include <CGAL/assertions.h>

#include<functional>

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
 * - an array of booleans that indicates for each curve to the right of 
 *   the event, whether it is already in the planar map or not. This is 
 *   used to speed insertions of curves into the planar map.
 * - an array of events that occur on the vertical curves that go through 
 *   this event. This is used instead of the array of points that is kept 
 *   in the base class.
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

  typedef typename CurveWrap::PmwxInsertInfo PmwxInsertInfo;

	typedef Event_less_functor<Self ,SweepLineTraits_2>  EventLess;
	typedef std::set<Self * , EventLess > VerticalXEventSet;
  
  typedef typename PmwxInsertInfo::Halfedge_handle Halfedge_handle;

  Pmwx_sweep_line_event():m_isCurveInPm(1)
  {}

  /*! Constructor */
  Pmwx_sweep_line_event(const Point_2 &point, Traits *traits) :
    Base(point, traits)
  {}


  /*! destructor */
  ~Pmwx_sweep_line_event()
  {}


  void init(const Point_2 &point, Traits *traits)
  {
    Base::init(point,traits);
  }

  PmwxInsertInfo *get_insert_info()
  {
    return &m_insertInfo;
  }

  void add_curve_to_right(SubCurve *curve)
  {
    if(Base::add_curve_to_right(curve))
      m_insertInfo.inc_right_curves_counter();
  }


  /*! Initialize the array that indicates wheter a curve to the right of the
   * event was already inserted into the planar map.
   */
  void init_right_curves()
  {
    unsigned int num_right_curves = get_num_right_curves();
    m_isCurveInPm.resize(num_right_curves);
    for ( unsigned int i = 0 ; i < num_right_curves ; i++ )
      m_isCurveInPm[i] = false;
  }
  
  /*! Caculates the number of halfedges in the planar map between the highest
   *  halfedge to the left of the event (which is stored in the insertInfo 
   *  member) and the position of the the specified curve around the vertex 
   *  in the planar map.
   *
   * @param curve a pointer to a curve that is going to be inserted 
   * @return the number of halfedges to skip before inserting the curve
   */
  int get_halfedge_jump_count(CurveWrap *curve)
  {
    int i = 0;
    int counter = 0;
    int skip = 0;
    unsigned int _size = m_isCurveInPm.size();
   
    for (unsigned int j = 0 ; j < _size ; j++ ) {
      if ( m_isCurveInPm[j] == true ) {
        skip++;
      }
    }
    skip--;  // now 'skip' holds the amount of the right curves of the event
		         // that are already inserted to the planar map  - 1 (minus 1)

    SubCurveIter iter = this->m_rightCurves.end();
    --iter;
    

	 
    // a boolean indictes if there is a vertical that was inserted to the map with endpoint as the
    // event point.
    /*bool exist_vertical = m_insertInfo.get_vertical_below_event_flag() ||
                          m_insertInfo.get_vertical_above_event_flag();*/
	   
    unsigned int num_left_curves = get_num_left_curves();
    for ( ; iter != m_rightCurves.begin() ; --iter )
    {
      if(curve == (*iter))
      {
        m_isCurveInPm[counter] = true;
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
      if ( m_isCurveInPm[counter] == true )
        i++;
      counter++;
    }

    CGAL_assertion(curve == (*iter));
    m_isCurveInPm[counter] = true;
    
    if ( num_left_curves == 0 )
      i--;
    return i;
  }

  /*! Returns true if the curve is the highest one among the right curves 
   *  that were already inserted into the planar map.
   */
  bool is_curve_largest(CurveWrap *curve)
  {
    int counter = 0;
    SubCurveIter iter = this->m_rightCurves.end();
    --iter;
   
    while(curve != (*iter))
    {
      if ( m_isCurveInPm[counter] == true )
        return false;
      counter++;
      --iter;
    }
    return true;
  }

private:
  PmwxInsertInfo m_insertInfo;
  std::vector<bool> m_isCurveInPm;
};

CGAL_END_NAMESPACE

#endif // CGAL_PMWX_SWEEP_LINE_EVENT_H
