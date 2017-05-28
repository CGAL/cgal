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
// Author(s) : Tali Zvi <talizvi@post.tau.ac.il>,
//             Baruch Zukerman <baruchzu@post.tau.ac.il>
//             Ron Wein <wein@post.tau.ac.il>
//             Efi Fogel <efifogel@gmail.com>

#ifndef CGAL_SURFACE_SWEEP_SUBCURVE_H
#define CGAL_SURFACE_SWEEP_SUBCURVE_H

#include <CGAL/license/Surface_sweep_2.h>

/*! \file
 *
 * Defintion of the Surface_sweep_subcurve class, which is an extended curve
 * type, referred to as Subcurve, used by the surface-sweep framework.
 *
 * The surface-sweep framework is implemented as a template that is
 * parameterized, among the other, by the Subcurve and Event types. That is,
 * instance types of Subcurve and Event must be available when the
 * surface-sweep template is instantiated.
 *
 * Surface_sweep_subcurve derives from an instance of the
 * No_overlap_surface_sweep_subcurve class template. The user is allowed to
 * introduce new types that derive from an instance of the
 * Surface_sweep_subcurve class template. However, some of the fields of this
 * template depends on the Subcurve type.  We use the curiously recurring
 * template pattern (CRTP) idiom to force the correct matching of these types.
 */

#include <CGAL/Surface_sweep_2/Surface_sweep_functors.h>
#include <CGAL/Surface_sweep_2/No_overlap_surface_sweep_subcurve.h>
#include <CGAL/Surface_sweep_2/Surface_sweep_event.h>
#include <CGAL/Multiset.h>
#include <CGAL/assertions.h>
#include <CGAL/Default.h>

namespace CGAL {

/*! \class Surface_sweep_subcurve_base
 *
 * This is the base class of the Surface_sweep_subcurve class template used by
 * the (CRTP) idiom.
 * \tparam GeometryTraits_2 the geometry traits.
 * \tparam Subcurve_ the subcurve actual type.
 *
 * The information contained in this class is:
 * - two pointers to subcurves that are the originating subcurves in case of
 *   an overlap, otherwise thay are both NULL.
 */
template <typename GeometryTraits_2, typename Subcurve_>
class Surface_sweep_subcurve_base :
    public No_overlap_surface_sweep_subcurve<GeometryTraits_2, Subcurve_>
{
public:
  typedef GeometryTraits_2                              Traits_2;
  typedef Subcurve_                                     Subcurve;

  typedef typename Traits_2::X_monotone_curve_2         X_monotone_curve_2;
  typedef No_overlap_surface_sweep_subcurve<Traits_2, Subcurve> Base;

  /*! Construct default.
   */
  Surface_sweep_subcurve_base() :
    m_orig_subcurve1(NULL),
    m_orig_subcurve2(NULL)
  {}

  /*! Construct from a curve.
   */
  Surface_sweep_subcurve_base(const X_monotone_curve_2& curve) :
    Base(curve),
    m_orig_subcurve1(NULL),
    m_orig_subcurve2(NULL)
  {}

protected:
  Subcurve* m_orig_subcurve1;           // The overlapping hierarchy
  Subcurve* m_orig_subcurve2;           // (relevant only in case of overlaps).
};

/*! \class Surface_sweep_subcurve
 *
 * This is a class template that wraps a traits curve of type
 * X_monotone_curve_2.  It contains data that is used when applying the sweep
 * algorithm on a set of x-monotone curves. This class derives from the
 * No_overlap_surface_sweep_subcurve class template.
 * \tparam GeometryTraits_2 the geometry traits.
 * \tparam Subcurve_ the type of the subcurve or Default. If the default is not
 *         overriden it implies that the type is
 *         No_overlap_surface_sweep_subcurve
 */
template <typename GeometryTraits_2, typename Subcurve_ = Default>
class Surface_sweep_subcurve :
    public Surface_sweep_subcurve_base
           <GeometryTraits_2,
            typename Default::Get<Subcurve_,
                                  Surface_sweep_subcurve
                                  <GeometryTraits_2, Subcurve_> >::type>
{
public:
  typedef GeometryTraits_2                                Traits_2;
  typedef typename Traits_2::X_monotone_curve_2           X_monotone_curve_2;
  typedef Surface_sweep_subcurve<Traits_2, Subcurve_>     Self;
  typedef typename Default::Get<Subcurve_, Self>::type    Subcurve;
  typedef Surface_sweep_subcurve_base<Traits_2, Subcurve> Base;

public:
  /*! Construct default.
   */
  Surface_sweep_subcurve() {}

  /*! Construct from a curve.
   */
  Surface_sweep_subcurve(const X_monotone_curve_2& curve) : Base(curve) {}

  /*! Destruct.
   */
  ~Surface_sweep_subcurve() {}

  /*! Get the subcurves that originate an overlap. */
  Subcurve* originating_subcurve1() { return Base::m_orig_subcurve1; }

  Subcurve* originating_subcurve2() { return Base::m_orig_subcurve2; }

  /*! Set the subcurves that originate an overlap. */
  void set_originating_subcurve1(Subcurve* orig_subcurve1)
  { Base::m_orig_subcurve1 = orig_subcurve1; }

  void set_originating_subcurve2(Subcurve* orig_subcurve2)
  { Base::m_orig_subcurve2 = orig_subcurve2; }

  /*! Get all the leaf-nodes in the hierarchy of overlapping subcurves. */
  template <typename OutputIterator>
  OutputIterator all_leaves(OutputIterator oi)
  {
    if (Base::m_orig_subcurve1 == NULL) {
      *oi++ = reinterpret_cast<Subcurve*>(this);
      return oi;
    }

    oi = Base::m_orig_subcurve1->all_leaves(oi);
    oi = Base::m_orig_subcurve2->all_leaves(oi);
    return oi;
  }

  /*! Check if the given subcurve is a node in the overlapping hierarchy. */
  bool is_inner_node(Subcurve* s)
  {
    if (this == s) return true;
    if (Base::m_orig_subcurve1 == NULL) return false;
    return (Base::m_orig_subcurve1->is_inner_node(s) ||
            Base::m_orig_subcurve2->is_inner_node(s));
  }

  /*! Check if the given subcurve is a leaf in the overlapping hierarchy. */
  bool is_leaf(Subcurve* s)
  {
    if (Base::m_orig_subcurve1 == NULL) return (this == s);
    return (Base::m_orig_subcurve1->is_leaf(s) ||
            Base::m_orig_subcurve2->is_leaf(s));
  }

  /*! Check if the two hierarchies contain the same leaf nodes. */
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

  /*! Check if the two hierarchies contain a common leaf node. */
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
    if (Base::m_orig_subcurve1 == NULL) {
      Subcurve* subcurve = reinterpret_cast<Subcurve*>(this);
      if (s->is_leaf(subcurve)) *oi++ = subcurve;
      return oi;
    }

    if (! s->is_inner_node(Base::m_orig_subcurve1))
      *oi++ = Base::m_orig_subcurve1;
    else oi++ = Base::m_orig_subcurve1->distinct_nodes(s, oi);

    if (! s->is_inner_node(Base::m_orig_subcurve2))
      *oi++ = Base::m_orig_subcurve2;
    else oi++ = Base::m_orig_subcurve2->distinct_nodes(s, oi);

    return oi;
  }

  /*! Get the depth of the overlap hierarchy. */
  unsigned int overlap_depth()
  {
    if (Base::m_orig_subcurve1 == NULL) return (1);

    unsigned int depth1 = Base::m_orig_subcurve1->overlap_depth();
    unsigned int depth2 = Base::m_orig_subcurve2->overlap_depth();
    if (depth1 > depth2) return (depth1 + 1);
    else return (depth2 + 1);
  }

#ifdef CGAL_SL_VERBOSE
  void Print() const;
#endif
};

#ifdef CGAL_SL_VERBOSE
  template<class Traits>
  void Surface_sweep_subcurve<Traits>::Print() const
  {
    std::cout << "Curve " << this
              << "  (" << this->last_curve() << ") "
              << " [sc1: " << this->originating_subcurve1()
              << ", sc2: " << this->originating_subcurve2() << "]";
  }
#endif

} //namespace CGAL

#endif
