// Copyright (c) 2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
//
// Author(s)     : Tali Zvi <talizvi@post.tau.ac.il>
//                 Baruch Zukerman <baruchzu@post.tau.ac.il>
//                 Efi Fogel <efif@post.tau.ac.il>

#ifndef CGAL_ARR_CONSTRUCTION_EVENT_H
#define CGAL_ARR_CONSTRUCTION_EVENT_H

/*! \file
 * Definition of the Arr_construction_event class-template.
 */

#include <CGAL/Sweep_line_2/Sweep_line_event.h>
#include <CGAL/assertions.h>
#include <vector>

namespace CGAL {

/*! \class Arr_construction_event
 *
 * Stores the data associated with an event.
 * In addition to the information stored in Sweep_line_event, when
 * constructing an arrangement, additional information is kept, in
 * order to speed insertion of curves into the planar map.
 *
 * The additional infomation contains the following:
 * - among the left curves of the event, we keep the highest halfedge that
 *   was inserted into the arrangement at any given time and when there no
 *   left curves, we keep the highest halfedge that was inseted to the right.
 *
 * Inherits from `Sweep_line_event`.
 * \sa `Sweep_line_event`
 */

template<class Traits_, class Subcurve_, class Arrangement_>
class Arr_construction_event :
  public Sweep_line_event<Traits_, Subcurve_>
{
public:

  typedef Traits_                                         Traits_2;
  typedef Subcurve_                                       Subcurve;
  typedef Arrangement_                                    Arrangement_2;
  typedef typename Arrangement_2::Vertex_handle           Vertex_handle;
  typedef typename Arrangement_2::Halfedge_handle         Halfedge_handle;

  typedef typename Traits_2::X_monotone_curve_2           X_monotone_curve_2;
  typedef typename Traits_2::Point_2                      Point_2;

  typedef Sweep_line_event<Traits_2,
                           Subcurve>                      Base;

  typedef Arr_construction_event<Traits_2,
                                 Subcurve,
                                 Halfedge_handle>         Self;

  typedef typename Base::Subcurve_container         Subcurve_container;
  typedef typename Base::Subcurve_iterator          Subcurve_iterator;
  typedef typename Base::Subcurve_reverse_iterator  Subcurve_reverse_iterator;

protected:

  // Data members:
  std::vector<bool>  m_isCurveInArr;          // Stores for each incident
                                              // subcurve whether it has been
                                              // inserted into the arrangement.

  Halfedge_handle    m_halfedge;              // A halfedge handle.
  Vertex_handle      m_vertex;                // A vertex handle.

  unsigned int       m_right_curves_counter;  // Number of subcurves defined
                                              // to the event's right that
                                              // haven't been added to the
                                              // arrangement, when that counter
                                              // is zero, we can deallocate the
                                              // event.

public:

  /*! Default constructor. */
  Arr_construction_event():
    m_halfedge(),
    m_vertex(),
    m_right_curves_counter(0)
  {}

  /*! Destructor */
  ~Arr_construction_event()
  {}

  /*! Add a curve to the right of the event. */
  std::pair<bool, Subcurve_iterator>
  add_curve_to_right (Subcurve *curve,
                      const Traits_2 * tr)
  {
    std::pair<bool,Subcurve_iterator> res =
      Base::add_curve_to_right(curve, tr);

    if(res.second != this->m_rightCurves.end() && res.first == false)
      ++m_right_curves_counter;

    return res;
  }

  /*! Add a curve pair to the right of the event. */
  std::pair<bool, Subcurve_iterator>
  add_curve_pair_to_right (Subcurve *sc1, Subcurve *sc2)
  {
    //increment twice the counter of right curves
    m_right_curves_counter+=2;
    return (Base::add_curve_pair_to_right(sc1, sc2));
  }

  /*! using the additional data that we store at the event, we compute
   *  how much we have to jump (he = he->next()->twin()) from the halfedge
   *  that is stored in the event, to the halefge that is previous to 'curve'
   *  that is about to be inserted into the arrangement.
   */
  int compute_halfedge_jump_count(Subcurve *curve)
  {
    int          i = 0;
    int          skip = 0;
    int          counter = 0;
    unsigned int j;

    for (j = 0 ; j < m_isCurveInArr.size() ; j++ )
    {
      if (m_isCurveInArr[j])
        skip++;
    }
    skip--;  // now 'skip' holds the amount of the right curves of the event
             // that are already inserted to the planar map  - 1 (minus 1)

    Subcurve_iterator  iter = this->m_rightCurves.end();
    size_t num_left_curves = this->number_of_left_curves();

    for (--iter; iter != this->m_rightCurves.begin() ; --iter, ++counter)
    {
      if (curve == (*iter))
      {
        m_isCurveInArr[counter] = true;

        if ((i == 0) && (num_left_curves == 0))
          return (skip);
        if (num_left_curves == 0)
          return (i - 1);

        return (i);
      }

      if (m_isCurveInArr[counter])
        i++;
    }

    CGAL_assertion(curve == (*iter));
    m_isCurveInArr[counter] = true;

    if (num_left_curves == 0)
      i--;

    return (i);
  }

  /*! return true iff 'curve' is the toppest curve among the halfedges
   *  to the right fo the event that were already were inserted to the
   * arrangement.
   */
  bool is_curve_largest (Subcurve *curve)
  {
    int counter = 0;

    Subcurve_reverse_iterator  rev_iter;
    for (rev_iter = this->m_rightCurves.rbegin();
         rev_iter != this->m_rightCurves.rend() && curve != (*rev_iter) ;
         ++rev_iter, ++ counter)
    {
      if(m_isCurveInArr[counter] == true)
         return false;
    }
    return true;
  }

  /*!
   * Resize the bit-vector indicating whether the incident curves are already
   * in the arrangement, and set all flags to false.
   */
  void init_subcurve_in_arrangement_flags (size_t n)
  {
    m_isCurveInArr.resize (n, false);
    return;
  }

  /*! Check if the i'th subcurve is in the arrangement. */
  bool is_subcurve_in_arrangement (unsigned int i) const
  {
    return (m_isCurveInArr[i]);
  }

  /*!
   * Set the flag indicating whether the i'th subcurve is in the arrangement.
   */
  void set_subcurve_in_arrangement (unsigned int i, bool flag)
  {
    m_isCurveInArr[i] = flag;
    return;
  }

  /*! Set the halfedge handle. */
  void set_halfedge_handle (Halfedge_handle h)
  {
    m_halfedge = h;
  }

  /*! Get the halfedge handle. */
  Halfedge_handle halfedge_handle() const
  {
    return m_halfedge;
  }

  /*! Set the vertex handle. */
  void set_vertex_handle (Vertex_handle v)
  {
    m_vertex = v;
  }

  /*! Get the vertex handle. */
  Vertex_handle vertex_handle() const
  {
    return m_vertex;
  }

  /*! Decrement the count of curves to the right that we have't done yet with
   *  (haven't been inserted to the arrangement). */
  unsigned int dec_right_curves_counter()
  {
    return (--m_right_curves_counter);
  }

  /*! Get the number of subcurves to the right of the event that we have't
   * done yet with (haven't been inserted to the arrangement).
   */
  unsigned int right_curves_counter() const
  {
    return (m_right_curves_counter);
  }

};

} //namespace CGAL

#endif
