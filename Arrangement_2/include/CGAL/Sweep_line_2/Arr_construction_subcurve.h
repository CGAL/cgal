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
// $URL: svn+ssh://golubevs@scm.gforge.inria.fr/svn/cgal/trunk/Arrangement_2/include/CGAL/Sweep_line_2/Arr_construction_curve.h $
// $Id: Arr_construction_curve.h 31509 2006-06-11 12:02:54Z baruchzu $
// 
//
// Author(s)     : Tali Zvi <talizvi@post.tau.ac.il>
//                 Baruch Zukerman <baruchzu@post.tau.ac.il>

#ifndef CGAL_ARR_CONSTRUCTION_SUBCURVE_H
#define CGAL_ARR_CONSTRUCTION_SUBCURVE_H

#include <CGAL/Sweep_line_2/Sweep_line_subcurve.h>

CGAL_BEGIN_NAMESPACE

/*! @class Arr_construction_subcurve 
 *  
 * a class that holds information about a curve that is added to 
 * the arrangement.
 * In addition to the information that is contained in Sweep_line_subcurve,
 * when an arrangement is constructed, a pointer to the last handled event  
 * on the curve is stored. This information is used to retrieve
 * hints when a subcurve of this curve is inserted into the planar map.
 *
 * Inherits from Sweep_line_subcurve
 * \sa Sweep_line_subcurve
 */

template<class SweepLineTraits_2>
class Arr_construction_subcurve : public Sweep_line_subcurve<SweepLineTraits_2>
{
public:
  typedef SweepLineTraits_2                        Traits;
  typedef typename Traits::Point_2                 Point_2;
  typedef typename Traits::X_monotone_curve_2      X_monotone_curve_2;

  typedef Sweep_line_subcurve<SweepLineTraits_2>   Base;
  typedef Arr_construction_subcurve<Traits>        Self;

  typedef typename Base::StatusLineIter            StatusLineIter;
  typedef void*                                    Event_ptr;


  Arr_construction_subcurve(): Base(),
                            m_lastEvent(0),
                            m_index(0)
  {}

  Arr_construction_subcurve(X_monotone_curve_2 &curve): Base( curve),
                                                     m_lastEvent(0),
                                                     m_index(0)
  {}

  
  void init(const X_monotone_curve_2 &curve)
  {
    Base::init(curve);
  }

  template<class SweepEvent>
  void set_left_event(SweepEvent* left)
  {
    Base::set_left_event(left);
    m_lastEvent = left;
  }

  void set_last_event(Event_ptr e) {
    m_lastEvent = e;
  }

  Event_ptr last_event() const {
    return m_lastEvent;
  }

  unsigned int index() const
  {
    return (m_index);
  }

  void set_index(unsigned int i)
  {
    m_index = i;
  }

  bool has_valid_index() const
  {
    return (m_index != 0);
  }

  void push_back_halfedge_index(unsigned int i)
  {
    m_haldedges_indexes.push_back(i);
  }

  void clear_haldedges_indexes()
  {
    m_haldedges_indexes.clear();
  }

  bool has_haldedges_indexes() const
  {
    return (!m_haldedges_indexes.empty());
  }

  std::list<unsigned int>& get_haldedges_indexes_list()
  {
    return (m_haldedges_indexes);
  }

protected:

  /*! the last event that was handled on the curve */
  Event_ptr  m_lastEvent;

  /*! index for a subcurve that may represent a hole (emarge from the left
   * most vertex of a hole, and its the upper most curve). other subcurves
   * will have 0 value  (invalid index)
   */
  unsigned int m_index;

  /*! indxes of all haldedges below that may represent a hole */
  std::list<unsigned int>  m_haldedges_indexes;
  
};


CGAL_END_NAMESPACE

#endif 
