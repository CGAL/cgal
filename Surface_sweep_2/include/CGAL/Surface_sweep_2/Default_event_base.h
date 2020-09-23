// Copyright (c) 2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Tali Zvi <talizvi@post.tau.ac.il>,
//             Baruch Zukerman <baruchzu@post.tau.ac.il>
//             Ron Wein <wein@post.tau.ac.il>
//             Efi Fogel <efif@gmail.com>

#ifndef CGAL_SURFACE_SWEEP_2_DEFAULT_EVENT_BASE_H
#define CGAL_SURFACE_SWEEP_2_DEFAULT_EVENT_BASE_H

#include <CGAL/license/Surface_sweep_2.h>

/*! \file
 *
 * Defintion of the Default_event_base class.
 */

#include <CGAL/Surface_sweep_2/No_overlap_event_base.h>

namespace CGAL {
namespace Surface_sweep_2 {

/*! \class Default_event_base
 *
 * A class associated with an event in a sweep line algorithm.
 * An intersection point in the sweep line algorithm is refered to as an event.
 * This class contains the information that is associated with any given
 * event point. This information contains the following:
 * - the actual point
 * - a list of curves that pass through the event point and defined to
 *   the left of the event point.
 * - a list of curves that pass through the event point and defined to
 *   the right of the event point.
 *
 * The class mostly exists to store information and does not have any
 * significant functionality otherwise.
 *
 */
template <typename GeometryTraits_2, typename Subcurve_>
class Default_event_base :
  public No_overlap_event_base<GeometryTraits_2, Subcurve_>
{
public:
  typedef GeometryTraits_2                              Geometry_traits_2;
  typedef Subcurve_                                     Subcurve;

private:
  typedef Geometry_traits_2                             Gt2;
  typedef No_overlap_event_base<Gt2, Subcurve>          Base;

public:
  typedef typename Base::X_monotone_curve_2             X_monotone_curve_2;
  typedef typename Base::Point_2                        Point_2;

  typedef typename Base::Left_side_category             Left_side_category;
  typedef typename Base::Bottom_side_category           Bottom_side_category;
  typedef typename Base::Top_side_category              Top_side_category;
  typedef typename Base::Right_side_category            Right_side_category;

  typedef typename Base::Subcurve_container             Subcurve_container;
  typedef typename Base::Subcurve_iterator              Subcurve_iterator;
  typedef typename Base::Subcurve_const_iterator     Subcurve_const_iterator;
  typedef typename Base::Subcurve_reverse_iterator   Subcurve_reverse_iterator;

public:
  /*! Default constructor. */
  Default_event_base() {}

  /*! Add a subcurve to the container of left curves. */
  void add_curve_to_left(Subcurve* curve)
  {
    // Look for the subcurve.
    bool curve_added = false;
    std::vector<Subcurve_iterator> left_curves_to_remove;
    for (Subcurve_iterator iter = this->left_curves_begin();
         iter != this->left_curves_end(); ++iter)
    {
      // Do nothing if the curve exists.
      if ((curve == *iter) || (*iter)->is_inner_node(curve)) return;

      // Replace the existing curve in case of overlap, only if the set of
      // ancesters of curve contains the set of ancesters of *iter
      if (curve->has_common_leaf(*iter)) {
        if (curve->number_of_original_curves() >
            (*iter)->number_of_original_curves())
        {
          if (curve->are_all_leaves_contained(*iter)) {
            if (curve_added) {
              left_curves_to_remove.push_back(iter);
              continue;
            }

            *iter = curve;
            curve_added = true;
          }
        }
        else {
          if ((*iter)->are_all_leaves_contained(curve)) {
            CGAL_assertion(!curve_added);
            return;
          }
        }
      }
    }

    for (typename std::vector<Subcurve_iterator>::iterator it =
         left_curves_to_remove.begin(); it != left_curves_to_remove.end(); ++it)
      this->left_curves_erase(*it);

    if (curve_added) return;

    // The curve does not exist; insert it to the container.
    this->push_back_curve_to_left(curve);
  }

  /*! Add a subcurve to the container of right curves. */
  std::pair<bool, Subcurve_iterator>
  add_curve_to_right(Subcurve* curve, const Gt2* tr)
  {
    if (! this->has_right_curves()) {
      this->push_back_curve_to_right(curve);
      return (std::make_pair(false, this->right_curves_begin()));
    }

    // Check if its an event at open boundary,
    // and if so then there is no overlap
    //(there cannot be two non-overlap curves at the same event at open
    // boundary).
    if (!this->is_closed())
      return (std::make_pair(true, this->right_curves_begin()));

    Subcurve_iterator iter = this->right_curves_begin();
    Comparison_result res;

    while ((res = tr->compare_y_at_x_right_2_object()
            (curve->last_curve(), (*iter)->last_curve(), this->point())) ==
           LARGER)
    {
      ++iter;
      if (iter == this->right_curves_end()) {
        this->m_right_curves.insert(iter, curve);
        return std::make_pair(false, --iter);
      }
    }

    //overlap !!
    if (res == EQUAL) return std::make_pair(true, iter);

    this->m_right_curves.insert(iter, curve);
    return std::make_pair(false, --iter);
  }

  Subcurve_iterator
  get_curve_after_on_right(Subcurve* curve)
  {
    Subcurve_iterator iter = this->right_curves_begin();
    for (Subcurve_iterator end = this->right_curves_end(); iter!=end; ++iter)
    {
      // TODO refine the condition
      if ( (*iter)->is_leaf(curve) || curve->is_leaf(*iter) || curve->has_common_leaf(*iter) )
        break;
    }
    CGAL_assertion( iter!=this->right_curves_end() );
    ++iter;
    return iter;
  }

  /*! Remove a curve from the set of left curves. */
  void remove_curve_from_left(Subcurve* curve)
  {
    for (Subcurve_iterator iter = this->left_curves_begin();
         iter != this->left_curves_end(); ++iter)
    {
      if (curve == *iter) {
        this->left_curves_erase(iter);
        return;
      }
    }
  }

  /*! Remove a curve from the set of right curves.
   */
  void remove_curve_from_right(Subcurve* curve)
  {
    for (Subcurve_iterator iter = this->right_curves_begin();
         iter != this->right_curves_end(); ++iter)
    {
      if ((curve == *iter) || curve->are_all_leaves_contained(*iter)) {
        this->right_curves_erase(iter);
        return;
      }
    }
  }

  bool is_right_curve_bigger(Subcurve* c1, Subcurve* c2, const Gt2* tr)
  {
    bool found_c1 = false;
    bool found_c2 = false;
    for (Subcurve_iterator iter = this->right_curves_begin();
         iter != this->right_curves_end(); ++iter)
    {
      if (!found_c1 && ((*iter == c1) || (*iter)->are_all_leaves_contained(c1)))
      {
        if (found_c2) return true;
        else found_c1 = true;
      }

      if (!found_c2 && ((*iter == c2) || (*iter)->are_all_leaves_contained(c2)))
      {
        if (found_c1) return false;
        else found_c2 = true;
      }
    }
    CGAL_assertion(!found_c1 || !found_c2);

    return tr->compare_y_at_x_right_2_object()
      (c1->last_curve(), c2->last_curve(), this->point()) == LARGER;
  }

  std::vector< std::pair<Subcurve*, Subcurve*> > overlaps_on_right;
};

} // namespace Surface_sweep_2
} // namespace CGAL

#endif
