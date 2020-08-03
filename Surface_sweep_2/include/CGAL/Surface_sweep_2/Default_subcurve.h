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
//             Efi Fogel <efifogel@gmail.com>

#ifndef CGAL_SURFACE_SWEEP_DEFAULT_SUBCURVE_H
#define CGAL_SURFACE_SWEEP_DEFAULT_SUBCURVE_H

#include <CGAL/license/Surface_sweep_2.h>

/*! \file
 *
 * Defintion of the Default_subcurve class, which is an extended curve
 * type, referred to as Subcurve, used by the surface-sweep framework.
 *
 * The surface-sweep framework is implemented as a template that is
 * parameterized, among the other, by the Subcurve and Event types. That is,
 * instance types of Subcurve and Event must be available when the
 * surface-sweep template is instantiated.
 *
 * Default_subcurve derives from an instance of the No_overlap_subcurve class
 * template. The user is allowed to introduce new types that derive from an
 * instance of the Default_subcurve class template. However, some of the fields
 * of this template depends on the Subcurve type.  We use the curiously
 * recurring template pattern (CRTP) idiom to force the correct matching of
 * these types.
 */

#include <CGAL/Surface_sweep_2/No_overlap_subcurve.h>
#include <CGAL/Multiset.h>
#include <CGAL/assertions.h>
#include <CGAL/Default.h>
#include <CGAL/Small_unordered_set.h>
#include <set>

namespace CGAL {
namespace Surface_sweep_2 {

/*! \class Default_subcurve_base
 *
 * This is the base class of the Default_subcurve class template used by
 * the (CRTP) idiom.
 * \tparam GeometryTraits_2 the geometry traits.
 * \tparam Event_ the event type.
 * \tparam Allocator_ a type of an element that is used to acquire/release
 *                    memory for elements of the event queue and the status
 *                    structure, and to construct/destroy the elements in that
 *                    memory. The type must meet the requirements of Allocator.
 * \tparam Subcurve_ the subcurve actual type.
 *
 * The information contained in this class is:
 * - two pointers to subcurves that are the originating subcurves in case of
 *   an overlap, otherwise thay are both nullptr.
 */
template <typename GeometryTraits_2, typename Event_, typename Allocator_,
          typename Subcurve_>
class Default_subcurve_base :
  public No_overlap_subcurve<GeometryTraits_2, Event_, Allocator_, Subcurve_>
{
public:
  typedef GeometryTraits_2                              Geometry_traits_2;
  typedef Subcurve_                                     Subcurve;
  typedef Event_                                        Event;

private:
  typedef Geometry_traits_2                             Gt2;
  typedef No_overlap_subcurve<Gt2, Event, Subcurve>     Base;
  typedef Default_subcurve_base<GeometryTraits_2, Event_, Allocator_, Subcurve_>
                                                        Self;
  typedef Small_unordered_set<Self*, 8>                 Intersected_set;

public:
  typedef typename Gt2::X_monotone_curve_2              X_monotone_curve_2;

  /*! Construct default.
   */
  Default_subcurve_base() :
    m_orig_subcurve1(nullptr),
    m_orig_subcurve2(nullptr)
  {}

  /*! Construct from a curve.
   */
  Default_subcurve_base(const X_monotone_curve_2& curve) :
    Base(curve),
    m_orig_subcurve1(nullptr),
    m_orig_subcurve2(nullptr)
  {}

protected:
  Subcurve* m_orig_subcurve1;           // The overlapping hierarchy
  Subcurve* m_orig_subcurve2;           // (relevant only in case of overlaps).
  Intersected_set m_intersected;


public:
  /*! Get the subcurves that originate an overlap. */
  Subcurve* originating_subcurve1() { return m_orig_subcurve1; }

  Subcurve* originating_subcurve2() { return m_orig_subcurve2; }

  const Subcurve* originating_subcurve1() const { return m_orig_subcurve1; }

  const Subcurve* originating_subcurve2() const { return m_orig_subcurve2; }

  bool intersection_exists (Self* other) { return !m_intersected.insert(other); }

  /*! Set the subcurves that originate an overlap. */
  void set_originating_subcurve1(Subcurve* orig_subcurve1)
  { m_orig_subcurve1 = orig_subcurve1; }

  void set_originating_subcurve2(Subcurve* orig_subcurve2)
  { m_orig_subcurve2 = orig_subcurve2; }

  /*! Get all the leaf-nodes in the hierarchy of overlapping subcurves. */
  template <typename OutputIterator>
  OutputIterator all_leaves(OutputIterator oi)
  {
    if (m_orig_subcurve1 == nullptr) {
      *oi++ = reinterpret_cast<Subcurve*>(this);
      return oi;
    }

    oi = m_orig_subcurve1->all_leaves(oi);
    oi = m_orig_subcurve2->all_leaves(oi);
    return oi;
  }

  /*! Check whether the given subcurve is a node in the overlapping hierarchy.
   */
  bool is_inner_node(Subcurve* s)
  {
    if (this == s) return true;
    if (m_orig_subcurve1 == nullptr) return false;
    return (m_orig_subcurve1->is_inner_node(s) ||
            m_orig_subcurve2->is_inner_node(s));
  }

  /*! Check whether the given subcurve is a leaf in the overlapping hierarchy.
   */
  bool is_leaf(Subcurve* s)
  {
    if (m_orig_subcurve1 == nullptr) return (this == s);
    return (m_orig_subcurve1->is_leaf(s) ||
            m_orig_subcurve2->is_leaf(s));
  }

  /*! Check whether other is entirely contained in the hierarchy of subcurves
   * of this.
   */
  bool are_all_leaves_contained(Subcurve* other)
  {
    std::set<Subcurve*> leaves_of_this;
    all_leaves(std::inserter(leaves_of_this, leaves_of_this.begin()));
    std::vector<Subcurve*> leaves_of_other;
    other->all_leaves(std::back_inserter(leaves_of_other));
    if (leaves_of_other.size() > leaves_of_this.size()) return false;
    for (typename std::vector<Subcurve*>::iterator it = leaves_of_other.begin();
         it != leaves_of_other.end(); ++it)
    {
      if (leaves_of_this.count(*it) == 0) return false;
    }

    return true;
  }

  /*! Check whether the two hierarchies contain the same leaf nodes. */
  bool has_same_leaves(Subcurve* s)
  {
    std::list<Subcurve*> my_leaves;
    std::list<Subcurve*> other_leaves;

    all_leaves(std::back_inserter(my_leaves));
    s->all_leaves(std::back_inserter(other_leaves));

    typename std::list<Subcurve*>::iterator iter;
    for (iter = my_leaves.begin(); iter != my_leaves.end(); ++iter) {
      if (std::find(other_leaves.begin(), other_leaves.end(), *iter) ==
          other_leaves.end())
        return false;
    }

    for (iter = other_leaves.begin(); iter != other_leaves.end(); ++iter) {
      if (std::find(my_leaves.begin(), my_leaves.end(), *iter) ==
          my_leaves.end())
        return false;
    }

    return true;
  }

  /*! Check whether the two hierarchies (s1+s2 considered as an overlapping
   * curve not already created) contain the same leaf nodes.
   */
  bool has_same_leaves(Subcurve* s1, Subcurve* s2)
  {
    std::list<Subcurve*> my_leaves;
    std::list<Subcurve*> other_leaves;

    all_leaves(std::back_inserter(my_leaves));
    s1->all_leaves(std::back_inserter(other_leaves));
    s2->all_leaves(std::back_inserter(other_leaves));

    typename std::list<Subcurve*>::iterator iter;
    for (iter = my_leaves.begin(); iter != my_leaves.end(); ++iter) {
      if (std::find(other_leaves.begin(), other_leaves.end(), *iter) ==
          other_leaves.end())
        return false;
    }

    for (iter = other_leaves.begin(); iter != other_leaves.end(); ++iter) {
      if (std::find(my_leaves.begin(), my_leaves.end(), *iter) ==
          my_leaves.end())
        return false;
    }

    return true;
  }

  /*! Check whether the two hierarchies contain a common leaf node. */
  bool has_common_leaf(Subcurve* s)
  {
    std::list<Subcurve*> my_leaves;
    std::list<Subcurve*> other_leaves;

    all_leaves(std::back_inserter(my_leaves));
    s->all_leaves(std::back_inserter(other_leaves));

    typename std::list<Subcurve*>::iterator iter;
    for (iter = my_leaves.begin(); iter != my_leaves.end(); ++iter) {
      if (std::find(other_leaves.begin(), other_leaves.end(), *iter) !=
          other_leaves.end())
        return true;
    }
    return false;
  }

  /*! Get all distinct nodes from the two hierarchies. */
  template <typename OutputIterator>
  OutputIterator distinct_nodes(Subcurve* s, OutputIterator oi)
  {
    if (m_orig_subcurve1 == nullptr) {
      Subcurve* subcurve = reinterpret_cast<Subcurve*>(this);
      if (s->is_leaf(subcurve)) *oi++ = subcurve;
      return oi;
    }

    if (! s->is_inner_node(m_orig_subcurve1)) *oi++ = m_orig_subcurve1;
    else oi++ = m_orig_subcurve1->distinct_nodes(s, oi);

    if (! s->is_inner_node(m_orig_subcurve2)) *oi++ = m_orig_subcurve2;
    else oi++ = m_orig_subcurve2->distinct_nodes(s, oi);

    return oi;
  }

  /*! Get the depth of the overlap hierarchy. */
  unsigned int overlap_depth()
  {
    if (m_orig_subcurve1 == nullptr) return (1);

    unsigned int depth1 = m_orig_subcurve1->overlap_depth();
    unsigned int depth2 = m_orig_subcurve2->overlap_depth();
    if (depth1 > depth2) return (depth1 + 1);
    else return (depth2 + 1);
  }

  /*! Get the number of input curves contributing to the subcurve */
  unsigned int number_of_original_curves() const
  {
    if (m_orig_subcurve1 == nullptr) return 1;
    unsigned int d1 = m_orig_subcurve1->number_of_original_curves();
    unsigned int d2 = m_orig_subcurve2->number_of_original_curves();
    return d1+d2;
  }
};

/*! \class Default_subcurve
 *
 * This is a class template that wraps a traits curve of type
 * X_monotone_curve_2.  It contains data that is used when applying the sweep
 * algorithm on a set of x-monotone curves. This class derives from the
 * No_overlap_subcurve class template.
 *
 * \tparam GeometryTraits_2 the geometry traits.
 * \tparam Event_ the event type.
 * \tparam Allocator_ a type of an element that is used to acquire/release
 *                    memory for elements of the event queue and the status
 *                    structure, and to construct/destroy the elements in that
 *                    memory. The type must meet the requirements of Allocator.
 * \tparam Subcurve_ the type of the subcurve or Default. If the default is not
 *         overriden it implies that the type is
 *         No_overlap_subcurve
 */
template <typename GeometryTraits_2, typename Event_,
          typename Allocator_ = CGAL_ALLOCATOR(int),
          typename Subcurve_ = Default>
class Default_subcurve :
  public Default_subcurve_base<GeometryTraits_2, Event_, Allocator_,
                               typename Default::Get<Subcurve_,
                                                     Default_subcurve<
                                                       GeometryTraits_2, Event_,
                                                       Allocator_,
                                                       Subcurve_> >::type>
{
public:
  typedef GeometryTraits_2                              Geometry_traits_2;
  typedef Event_                                        Event;
  typedef Allocator_                                    Allocator;

private:
  typedef Geometry_traits_2                             Gt2;
  typedef Default_subcurve<Gt2, Event, Allocator, Subcurve_>
                                                        Self;
  typedef typename Default::Get<Subcurve_, Self>::type  Subcurve;
  typedef Default_subcurve_base<Gt2, Event, Allocator, Subcurve>
                                                        Base;

public:
  typedef typename Gt2::X_monotone_curve_2              X_monotone_curve_2;

public:
  /*! Construct default.
   */
  Default_subcurve() {}

  /*! Construct from a curve.
   */
  Default_subcurve(const X_monotone_curve_2& curve) : Base(curve) {}

  /*! Destruct.
   */
  ~Default_subcurve() {}

#ifdef CGAL_SS_VERBOSE
  void Print() const;
#endif
};

#ifdef CGAL_SS_VERBOSE
template <typename Gt2, typename Evt, typename Allocator, typename Scv>
void Default_subcurve<Gt2, Evt, Allocator, Scv>::Print() const
{
  std::cout << "Curve " << this
            << "  (" << this->last_curve() << ") "
            << " [sc1: " << this->originating_subcurve1()
            << ", sc2: " << this->originating_subcurve2() << "]";
}
#endif

} // namespace Surface_sweep_2
} // namespace CGAL

#endif
