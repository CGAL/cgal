// Copyright (c) 2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Tali Zvi <talizvi@post.tau.ac.il>
//                 Baruch Zukerman <baruchzu@post.tau.ac.il>
//                 Efi Fogel <efif@post.tau.ac.il>

#ifndef CGAL_ARR_CONSTRUCTION_EVENT_BASE_H
#define CGAL_ARR_CONSTRUCTION_EVENT_BASE_H

#include <CGAL/license/Arrangement_on_surface_2.h>

/*! \file
 *
 * Definition of the Arr_construction_event_base class-template.
 */

#include <vector>

#include <CGAL/Surface_sweep_2/Default_event_base.h>
#include <CGAL/assertions.h>

namespace CGAL {

namespace Ss2 = Surface_sweep_2;

/*! \class Arr_construction_event_base
 *
 * This template represents an event used by the surface-sweep framework.  It
 * inherits either from `Default_event_base` (the default) or
 * 'No_overlap_event_base' depending on whether the curve may overlap or not.
 * It stores the data associated with an event in addition to the information
 * stored in its base class. When constructing an arrangement, additional
 * information is stored, in order to expedite the insertion of curves into the
 * arrangement.
 *
 * The additional infomation contains the following:
 * - among the left curves of the event, we keep the highest halfedge that
 *   was inserted into the arrangement at any given time and when there are no
 *   left curves, we keep the highest halfedge that was inseted to the right.
 *
 * \tparam GeometryTraits_2 the geometry traits.
 * \tparam Allocator_ a type of an element that is used to acquire/release
 *                    memory for elements of the event queue and the status
 *                    structure, and to construct/destroy the elements in that
 *                    memory. The type must meet the requirements of Allocator.
 * \tparam SurfaceSweepEvent a template, an instance of which is used as the
 *                           base class.
 *
 * \sa `Default_event_base`
 * \sa `No_overlap_event_base`
 */
template <typename GeometryTraits_2, typename Subcurve_, typename Arrangement_,
          template <typename, typename>
          class SurfaceSweepEvent = Ss2::Default_event_base>
class Arr_construction_event_base :
  public SurfaceSweepEvent<GeometryTraits_2, Subcurve_>
{
public:
  typedef GeometryTraits_2                              Geometry_traits_2;
  typedef Subcurve_                                     Subcurve;
  typedef Arrangement_                                  Arrangement_2;

  typedef typename Arrangement_2::Vertex_handle         Vertex_handle;
  typedef typename Arrangement_2::Halfedge_handle       Halfedge_handle;

private:
  typedef Geometry_traits_2                             Gt2;
  typedef SurfaceSweepEvent<Gt2, Subcurve>              Base;
  typedef Arr_construction_event_base<Gt2, Subcurve, Halfedge_handle,
                                 SurfaceSweepEvent>     Self;

public:
  typedef typename Gt2::X_monotone_curve_2              X_monotone_curve_2;
  typedef typename Gt2::Point_2                         Point_2;

  typedef typename Base::Subcurve_container             Subcurve_container;
  typedef typename Base::Subcurve_iterator              Subcurve_iterator;
  typedef typename Base::Subcurve_reverse_iterator    Subcurve_reverse_iterator;

protected:
  // Data members:
  std::vector<bool> m_isCurveInArr;           // Stores for each incident
                                              // subcurve whether it has been
                                              // inserted into the arrangement.

  Halfedge_handle m_halfedge;                 // A halfedge handle.
  Vertex_handle m_vertex;                     // A vertex handle.

  unsigned int m_right_curves_counter;        // Number of subcurves defined
                                              // to the event's right that
                                              // haven't been added to the
                                              // arrangement, when that counter
                                              // is zero, we can deallocate the
                                              // event.

public:
  /*! Default constructor. */
  Arr_construction_event_base():
    m_halfedge(),
    m_vertex(),
    m_right_curves_counter(0)
  {}

  /*! Destructor */
  ~Arr_construction_event_base() {}

  /*! Add a curve to the right of the event. */
  std::pair<bool, Subcurve_iterator>
  add_curve_to_right(Subcurve* curve, const Gt2* tr)
  {
    std::pair<bool,Subcurve_iterator> res = Base::add_curve_to_right(curve, tr);

    if (res.second != this->right_curves_end() && res.first == false)
      ++m_right_curves_counter;

    return res;
  }

  /*! Add a curve pair to the right of the event. */
  std::pair<bool, Subcurve_iterator>
  add_curve_pair_to_right(Subcurve* sc1, Subcurve* sc2)
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
  int compute_halfedge_jump_count(Subcurve* curve)
  {
    int i = 0;
    int skip = 0;
    int counter = 0;
    unsigned int j;

    for (j = 0 ; j < m_isCurveInArr.size() ; j++) {
      if (m_isCurveInArr[j]) skip++;
    }
    skip--;  // now 'skip' holds the amount of the right curves of the event
             // that are already inserted to the planar map  - 1 (minus 1)

    Subcurve_iterator iter = this->right_curves_end();
    size_t num_left_curves = this->number_of_left_curves();

    for (--iter; iter != this->right_curves_begin() ; --iter, ++counter) {
      if (curve == (*iter)) {
        m_isCurveInArr[counter] = true;

        if ((i == 0) && (num_left_curves == 0)) return (skip);
        if (num_left_curves == 0) return (i - 1);
        return (i);
      }

      if (m_isCurveInArr[counter]) i++;
    }

    CGAL_assertion(curve == (*iter));
    m_isCurveInArr[counter] = true;

    if (num_left_curves == 0) i--;

    return (i);
  }

  /*! Return true iff 'curve' is the toppest curve among the halfedges
   *  to the right fo the event that were already were inserted to the
   * arrangement.
   */
  bool is_curve_largest(Subcurve *curve)
  {
    int counter = 0;

    Subcurve_reverse_iterator  rev_iter;
    for (rev_iter = this->right_curves_rbegin();
         rev_iter != this->right_curves_rend() && curve != (*rev_iter) ;
         ++rev_iter, ++ counter)
    {
      if (m_isCurveInArr[counter] == true) return false;
    }
    return true;
  }

  /*! Resize the bit-vector indicating whether the incident curves are already
   * in the arrangement, and set all flags to false.
   */
  void init_subcurve_in_arrangement_flags(size_t n)
  { m_isCurveInArr.resize(n, false); }

  /*! Check if the i'th subcurve is in the arrangement. */
  bool is_subcurve_in_arrangement(unsigned int i) const
  { return (m_isCurveInArr[i]); }

  /*! Set the flag indicating whether the i'th subcurve is in the arrangement.
   */
  void set_subcurve_in_arrangement(unsigned int i, bool flag)
  { m_isCurveInArr[i] = flag; }

  /*! Set the halfedge handle. */
  void set_halfedge_handle(Halfedge_handle h) { m_halfedge = h; }

  /*! Get the halfedge handle. */
  Halfedge_handle halfedge_handle() const { return m_halfedge; }

  /*! Set the vertex handle. */
  void set_vertex_handle(Vertex_handle v) { m_vertex = v; }

  /*! Get the vertex handle. */
  Vertex_handle vertex_handle() const { return m_vertex; }

  /*! Decrement the count of curves to the right that we have't done yet with
   *  (haven't been inserted to the arrangement). */
  unsigned int dec_right_curves_counter() { return (--m_right_curves_counter); }

  /*! Get the number of subcurves to the right of the event that we have't
   * done yet with (haven't been inserted to the arrangement).
   */
  unsigned int right_curves_counter() const { return (m_right_curves_counter); }
};

} // namespace CGAL

#endif
