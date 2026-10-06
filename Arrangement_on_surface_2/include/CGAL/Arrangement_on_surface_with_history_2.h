// Copyright (c) 2007,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s): Ron Wein         <wein@post.tau.ac.il>
//            Efi Fogel        <efif@post.tau.ac.il>
//            Baruch Zukerman  <baruchzu@post.tau.ac.il>

#ifndef CGAL_ARRANGEMENT_ON_SURFACE_WITH_HISTORY_2_H
#define CGAL_ARRANGEMENT_ON_SURFACE_WITH_HISTORY_2_H

#include <CGAL/license/Arrangement_on_surface_2.h>

#include <CGAL/disable_warnings.h>

/*! \file
 * The header file for the Arrangement_on_surface_with_history_2 class.
 */

#include <functional>
#include <list>
#include <map>
#include <memory>
#include <set>
#include <vector>

#include <CGAL/Arrangement_on_surface_2.h>
#include <CGAL/Arr_overlay_2.h>
#include <CGAL/Arr_consolidated_curve_data_traits_2.h>
#include <CGAL/In_place_list.h>
#include <CGAL/Arrangement_2/Arr_with_history_accessor.h>

namespace CGAL {

/*! \class Arrangement_on_surface_with_history_2
 * A class representing planar subdivisions induced by a set of arbitrary
 * input planar curves. These curves are split and form a set of sweepable
 * (`x`-monotone) and pairwise interior-disjoint curves that are associated
 * with the arrangement edges.
 * The `Arrangement_on_surface_with_history_2` class enables tracking the input
 * curve(s) that originated each such subcurve. It also enables keeping track
 * of the edges that resulted from each input curve.
 * The `GeomTraits` parameter corresponds to a traits class that defines the
 * `Point_2`, `X_monotone_curve_2` and `Curve_2` types and implements the geometric
 * predicates and constructions for the family of curves it defines.
 * The `TopolTraits` parameter corresponds to a topology-traits class that defines
 * the topological structure of the surface. Note that the geometry traits
 * class should also be aware of the kind of surface on which its curves and
 * points are defined.
 */

namespace internal {

/*! computes the base class of `Arrangement_on_surface_with_history_2`: an arrangement on surface instantiated with
 * the consolidated-curve-data traits and with the topology traits and the \dcel rebound to these traits.
 */
template <typename GeomTraits, typename TopolTraits>
struct Aos_with_history_base_2 {
  using Data_traits = Arr_consolidated_curve_data_traits_2<GeomTraits, typename GeomTraits::Curve_2*>;
  using Data_dcel = typename TopolTraits::Dcel::template rebind<Data_traits>::other;
  using Data_topol_traits = typename TopolTraits::template rebind<Data_traits, Data_dcel>::other;
  using type = Arrangement_on_surface_2<Data_traits, Data_topol_traits>;
};

} // namespace internal

template <typename GeomTraits_, typename TopolTraits_>
class Arrangement_on_surface_with_history_2 :
    public internal::Aos_with_history_base_2<GeomTraits_, TopolTraits_>::type {
  // The computation of the base class and of the types it is instantiated with.
  using Base_helper = internal::Aos_with_history_base_2<GeomTraits_, TopolTraits_>;

public:
  using Geometry_traits_2 = GeomTraits_;
  using Base_topology_traits = TopolTraits_;

  /*! a shared pointer to the (immutable) geometry traits. */
  using Shared_geometry_traits = std::shared_ptr<const Geometry_traits_2>;

private:
  using Self = Arrangement_on_surface_with_history_2<Geometry_traits_2, Base_topology_traits>;

public:
  using Point_2 = typename Geometry_traits_2::Point_2;
  using Curve_2 = typename Geometry_traits_2::Curve_2;
  using X_monotone_curve_2 = typename Geometry_traits_2::X_monotone_curve_2;

protected:
  friend class Arr_accessor<Self>;
  friend class Arr_with_history_accessor<Self>;

  // The data-traits class, based on Geometry_traits_2.
  using Data_traits_2 = typename Base_helper::Data_traits;
  using Data_curve_2 = typename Data_traits_2::Curve_2;
  using Data_x_curve_2 = typename Data_traits_2::X_monotone_curve_2;
  using Data_iterator = typename Data_traits_2::Data_iterator;

  // The \dcel and the topology traits rebound to the data-traits class.
  using Data_dcel = typename Base_helper::Data_dcel;
  using Data_top_traits = typename Base_helper::Data_topol_traits;

  // The arrangement with history is based on the representation of an
  // arrangement, templated by the data-traits class and the rebound \dcel.
  using Base_arr_2 = typename Base_helper::type;

  // A shared pointer to the data traits, as stored by the base arrangement.
  using Shared_data_traits = typename Base_arr_2::Shared_geometry_traits;

  /*! obtains a pointer to the data-traits view of the given geometry traits
   * that shares ownership with it (using the aliasing constructor). The
   * downcast is the same one performed by the raw-pointer constructor.
   */
  static Shared_data_traits _data_traits(const Shared_geometry_traits& tr) {
    CGAL_precondition(tr != nullptr);
    return Shared_data_traits(tr, static_cast<const Data_traits_2*>(tr.get()));
  }

public:
  using Traits_adaptor_2 = Arr_traits_adaptor_2<Data_traits_2>;
  using Topology_traits = Data_top_traits;
  using Base_arrangement_2 = Base_arr_2;

  // Types inherited from the base arrangement class:
  using Size = typename Base_arr_2::Size;

  using Vertex = typename Base_arr_2::Vertex;
  using Halfedge = typename Base_arr_2::Halfedge;
  using Face = typename Base_arr_2::Face;

  using Vertex_iterator = typename Base_arr_2::Vertex_iterator;
  using Vertex_const_iterator = typename Base_arr_2::Vertex_const_iterator;
  using Halfedge_iterator = typename Base_arr_2::Halfedge_iterator;
  using Halfedge_const_iterator = typename Base_arr_2::Halfedge_const_iterator;
  using Edge_iterator = typename Base_arr_2::Edge_iterator;
  using Edge_const_iterator = typename Base_arr_2::Edge_const_iterator;
  using Face_iterator = typename Base_arr_2::Face_iterator;
  using Face_const_iterator = typename Base_arr_2::Face_const_iterator;
  using Halfedge_around_vertex_circulator = typename Base_arr_2::Halfedge_around_vertex_circulator;
  using Halfedge_around_vertex_const_circulator = typename Base_arr_2::Halfedge_around_vertex_const_circulator;
  using Ccb_halfedge_circulator = typename Base_arr_2::Ccb_halfedge_circulator;
  using Ccb_halfedge_const_circulator = typename Base_arr_2::Ccb_halfedge_const_circulator;
  using Outer_ccb_iterator = typename Base_arr_2::Outer_ccb_iterator;
  using Outer_ccb_const_iterator = typename Base_arr_2::Outer_ccb_const_iterator;
  using Inner_ccb_iterator = typename Base_arr_2::Inner_ccb_iterator;
  using Inner_ccb_const_iterator = typename Base_arr_2::Inner_ccb_const_iterator;
  using Isolated_vertex_iterator = typename Base_arr_2::Isolated_vertex_iterator;
  using Isolated_vertex_const_iterator = typename Base_arr_2::Isolated_vertex_const_iterator;

  using Vertex_handle = typename Base_arr_2::Vertex_handle;
  using Vertex_const_handle = typename Base_arr_2::Vertex_const_handle;
  using Halfedge_handle = typename Base_arr_2::Halfedge_handle;
  using Halfedge_const_handle = typename Base_arr_2::Halfedge_const_handle;
  using Face_handle = typename Base_arr_2::Face_handle;
  using Face_const_handle = typename Base_arr_2::Face_const_handle;

protected:
  /*! compares two halfedge handles by address.
   */
  struct Less_halfedge_handle {
    bool operator()(Halfedge_handle h1, Halfedge_handle h2) const { return std::less<const Halfedge*>()(&*h1, &*h2); }
  };

  // Forward declaration:
  class Curve_halfedges_observer;

public:
  /*! \class Curve_halfedges
   * Extension of a curve with the set of edges that it induces.
   * Each edge is represented by one of the halfedges.
   */
  class Curve_halfedges :
    public Curve_2,
    public In_place_list_base<Curve_halfedges> {
    using Gt = Geometry_traits_2;
    using Btt = Base_topology_traits;
    using Aos_wh = Arrangement_on_surface_with_history_2<Gt, Btt>;

    friend class Curve_halfedges_observer;
    friend class Arrangement_on_surface_with_history_2<Gt, Btt>;
    friend class Arr_with_history_accessor<Aos_wh>;

  private:
    using Halfedges_set = std::set<Halfedge_handle, Less_halfedge_handle>;

    // Data members:
    Halfedges_set m_halfedges;

  public:
    /*! constructs default. */
    Curve_halfedges() = default;

    /*! constructs from a given curve, with an empty set of edges. */
    explicit Curve_halfedges(const Curve_2& curve) : Curve_2(curve) {}

    using iterator = typename Halfedges_set::iterator;
    using const_iterator = typename Halfedges_set::const_iterator;

  private:
    /*! obtains the number of edges induced by the curve. */
    Size size() const { return m_halfedges.size(); }

    /*! obtains an iterator for the first edge in the set (const version). */
    const_iterator begin() const { return m_halfedges.begin(); }

    /*! obtains an iterator for the first edge in the set (non-const version). */
    iterator begin() { return m_halfedges.begin(); }

    /*! obtains a past-the-end iterator for the set edges (const version). */
    const_iterator end() const { return m_halfedges.end(); }

    /*! obtains a past-the-end iterator for the set edges (non-const version). */
    iterator end() { return m_halfedges.end(); }

    /*! inserts an edge to the set. */
    iterator _insert(Halfedge_handle he) {
      auto res = m_halfedges.insert(he);
      CGAL_assertion(res.second);
      return res.first;
    }

    /*! erases an edge, given by its position, from the set. */
    void erase(iterator pos) { m_halfedges.erase(pos); }

    /*! erases an edge from the set. */
    void _erase(Halfedge_handle he) {
      std::size_t res = m_halfedges.erase(he);
      if (res == 0) res = m_halfedges.erase(he->twin());
      CGAL_assertion(res != 0);
    }

    /*! clears the edges set. */
    void clear() { m_halfedges.clear(); }
  };

protected:
  using Curves_alloc = CGAL_ALLOCATOR(Curve_halfedges);
  using Curve_halfedges_list = In_place_list<Curve_halfedges, false>;

  /*! \class Curve_halfedges_observer
   * Observer for the base arrangement. It keeps track of all local changes
   * involving edges and updates the list of halfedges associated with the
   * input curves accordingly.
   */
  class Curve_halfedges_observer : public Base_arr_2::Observer {
  public:
    using Base_aos = typename Base_arr_2::Base_aos;

    using Vertex_handle = typename Base_aos::Vertex_handle;
    using Halfedge_handle = typename Base_aos::Halfedge_handle;
    using X_monotone_curve_2 = typename Base_aos::X_monotone_curve_2;

    /*! notifies after the creation of a new edge.
     * \param e A handle to one of the twin halfedges that were created.
     */
    void after_create_edge(Halfedge_handle e) override { _register_edge(e); }

    /*! notifies before the modification of an existing edge.
     * \param e A handle to one of the twin halfedges to be updated.
     * \param c The `x`-monotone curve to be associated with the edge.
     */
    void before_modify_edge(Halfedge_handle e, const X_monotone_curve_2& /* c */) override { _unregister_edge(e); }

    /*! notifies after an edge was modified.
     * \param e A handle to one of the twin halfedges that were updated.
     */
    void after_modify_edge(Halfedge_handle e) override { _register_edge(e); }

    /*! notifies before the splitting of an edge into two.
     * \param e A handle to one of the existing halfedges.
     * \param c1 The `x`-monotone curve to be associated with the first edge.
     * \param c2 The `x`-monotone curve to be associated with the second edge.
     */
    void before_split_edge(Halfedge_handle e, Vertex_handle /* v */,
                           const X_monotone_curve_2& /* c1 */, const X_monotone_curve_2& /* c2 */) override
    { _unregister_edge(e); }

    /*! notifies after an edge was split.
     * \param e1 A handle to one of the twin halfedges forming the first edge.
     * \param e2 A handle to one of the twin halfedges forming the second edge.
     */
    void after_split_edge(Halfedge_handle e1, Halfedge_handle e2) override {
      _register_edge(e1);
      _register_edge(e2);
    }

    /*! notifies before the merging of two edges.
     * \param e1 A handle to one of the halfedges forming the first edge.
     * \param e2 A handle to one of the halfedges forming the second edge.
     * \param c The `x`-monotone curve to be associated with the merged edge.
     */
    void before_merge_edge(Halfedge_handle e1, Halfedge_handle e2, const X_monotone_curve_2& /* c */) override {
      _unregister_edge(e1);
      _unregister_edge(e2);
    }

    /*! notifies after an edge was merged.
     * \param e A handle to one of the twin halfedges forming the merged edge.
     */
    void after_merge_edge(Halfedge_handle e) override { _register_edge(e); }

    /*! notifies before the removal of an edge.
     * \param e A handle to one of the twin halfedges to be deleted.
     */
    void before_remove_edge(Halfedge_handle e) override { _unregister_edge(e); }

  private:
    /*! registers the given halfedge in the set(s) associated with its curve. */
    void _register_edge(Halfedge_handle e)
    { for (auto* p_cv : e->curve().data()) static_cast<Curve_halfedges*>(p_cv)->_insert(e); }

    /*! unregisters the given halfedge from the set(s) associated with its curve. */
    void _unregister_edge(Halfedge_handle e)
    { for (auto* p_cv : e->curve().data()) static_cast<Curve_halfedges*>(p_cv)->_erase(e); }
  };

  // Data members:
  Curves_alloc m_curves_alloc;
  Curve_halfedges_list m_curves;
  Curve_halfedges_observer m_observer;

public:
  using Curve_iterator = typename Curve_halfedges_list::iterator;
  using Curve_const_iterator = typename Curve_halfedges_list::const_iterator;

  using Curve_handle = Curve_iterator;
  using Curve_const_handle = Curve_const_iterator;

  /// \name Constructors.
  //@{

  /*! constructs default. */
  Arrangement_on_surface_with_history_2();

  /*! constructs copy. */
  Arrangement_on_surface_with_history_2(const Self& arr);

  /*! constructs given a shared geometry-traits object. The arrangement (co-)owns the traits. */
  explicit Arrangement_on_surface_with_history_2(Shared_geometry_traits tr);

  /*! constructs from a traits object. The caller retains ownership of the traits
   * and must keep it alive as long as the arrangement (or any copy of it) exists.
   */
  Arrangement_on_surface_with_history_2(const Geometry_traits_2* tr);
  //@}

  /// \name Assignment functions.
  //@{

  /*! assigns. */
  Self& operator=(const Self& arr);

  /*! assigns an arrangement with history. */
  void assign(const Self& arr);
  //@}

  /// \name Destruction functions.
  //@{

  /*! destructs. */
  virtual ~Arrangement_on_surface_with_history_2();

  /*! clears the arrangement. */
  virtual void clear();
  //@}

  /*! obtains the geometry-traits object. */
  const Geometry_traits_2* geometry_traits() const { return this->m_geom_traits.get(); }

  /*! obtains a shared pointer to the geometry-traits object. If the
   * arrangement was constructed from a raw pointer, the returned pointer is
   * non-owning.
   */
  Shared_geometry_traits shared_geometry_traits() const { return Base_arr_2::shared_geometry_traits(); }

  /// \name Traversal of the arrangement curves.
  //@{
  Size number_of_curves() const { return m_curves.size(); }

  Curve_iterator curves_begin() { return m_curves.begin(); }

  Curve_iterator curves_end() { return m_curves.end(); }

  Curve_const_iterator curves_begin() const { return m_curves.begin(); }

  Curve_const_iterator curves_end() const { return m_curves.end(); }
  //@}

  /*! \class Originating_curve_iterator
   * An iterator over the curves that originated an edge, defined as a derived class to make it convertible to the
   * curve iterator type.
   */
  class Originating_curve_iterator :
    public I_Dereference_iterator<Data_iterator, Curve_2, typename Data_iterator::difference_type,
                                  typename Data_iterator::iterator_category> {
    using Base = I_Dereference_iterator<Data_iterator, Curve_2, typename Data_iterator::difference_type,
                                        typename Data_iterator::iterator_category>;

  public:
    Originating_curve_iterator() = default;

    Originating_curve_iterator(Data_iterator iter) : Base(iter) {}

    // Casting to a curve iterator.
    operator Curve_iterator() const {
      Curve_halfedges* p_cv = static_cast<Curve_halfedges*>(this->ptr());
      return Curve_iterator(p_cv);
    }

    operator Curve_const_iterator() const {
      const Curve_halfedges* p_cv = static_cast<Curve_halfedges*>(this->ptr());
      return Curve_const_iterator(p_cv);
    }
  };

  /// \name Traversal of the origin curves of an edge.
  //@{
  Size number_of_originating_curves(Halfedge_const_handle e) const { return e->curve().data().size(); }

  Originating_curve_iterator originating_curves_begin(Halfedge_const_handle e) const
  { return Originating_curve_iterator(e->curve().data().begin()); }

  Originating_curve_iterator originating_curves_end(Halfedge_const_handle e) const
  { return Originating_curve_iterator(e->curve().data().end()); }
  //@}

  using Induced_edge_iterator = typename Curve_halfedges::const_iterator;

  /// \name Traversal of the edges induced by a curve.
  //@{
  Size number_of_induced_edges(Curve_const_handle c) const { return c->size(); }

  Induced_edge_iterator induced_edges_begin(Curve_const_handle c) const { return c->begin(); }

  Induced_edge_iterator induced_edges_end(Curve_const_handle c) const { return c->end(); }
  //@}

  /// \name Manipulating edges.
  //@{

  /*! splits a given edge into two at the given split point.
   * \param e The edge to split (one of the pair of twin halfedges).
   * \param p The split point.
   * \pre p lies in the interior of the curve associated with e.
   * \return A handle for the halfedge whose source is the source of the
   *         original halfedge e, and whose target is the split point.
   */
  Halfedge_handle split_edge(Halfedge_handle e, const Point_2& p);

  /*! merges two edges to form a single edge.
   * \param e1 The first edge to merge (one of the pair of twin halfedges).
   * \param e2 The second edge to merge (one of the pair of twin halfedges).
   * \pre e1 and e2 must have a common end-vertex of degree 2 and must
   *      be mergeable.
   * \return A handle for the merged halfedge.
   */
  Halfedge_handle merge_edge(Halfedge_handle e1, Halfedge_handle e2);

  /*! checks if two edges can be merged to a single edge.
   * \param e1 The first edge (one of the pair of twin halfedges).
   * \param e2 The second edge (one of the pair of twin halfedges).
   * \return true iff e1 and e2 are mergeable.
   */
  bool are_mergeable(Halfedge_const_handle e1, Halfedge_const_handle e2) const;
  //@}

protected:
  /// \name Curve management.
  //@{

  /*! maps the curves of a source arrangement with history to their duplicates in this arrangement. The keys are the
   * addresses of the `Curve_2` subobjects, which are also the values stored with the edges. They are independent of
   * the type of the source arrangement, which may have a different topology traits (e.g., a different \dcel).
   */
  using Curve_map = std::map<const Curve_2*, Curve_halfedges*>;

  /*! allocates an extended curve with an empty set of edges, and appends it to the list of curves.
   * \param cv The curve.
   * \return A pointer to the new extended curve.
   */
  Curve_halfedges* _new_curve_halfedges(const Curve_2& cv) {
    Curve_halfedges* p_cv = m_curves_alloc.allocate(1);
    std::allocator_traits<Curves_alloc>::construct(m_curves_alloc, p_cv, cv);
    m_curves.push_back(*p_cv);
    return p_cv;
  }

  /*! removes an extended curve from the list of curves and deallocates it.
   * \param p_cv A pointer to the extended curve.
   */
  void _delete_curve_halfedges(Curve_halfedges* p_cv) {
    m_curves.erase(p_cv);
    std::allocator_traits<Curves_alloc>::destroy(m_curves_alloc, p_cv);
    m_curves_alloc.deallocate(p_cv, 1);
  }

  /*! duplicates a range of curves of a source arrangement with history, and records each duplicate in a map.
   * \param begin An iterator to the first curve in the range.
   * \param end A past-the-end iterator for the range.
   * \param cv_map Output: The map from the source curves to their duplicates.
   */
  template <typename CurveIterator>
  void _duplicate_curves(CurveIterator begin, CurveIterator end, Curve_map& cv_map) {
    for (auto it = begin; it != end; ++it) {
      const Curve_2* p_cv = &*it;
      cv_map.emplace(p_cv, _new_curve_halfedges(*p_cv));
    }
  }

  /*! redirects the curve pointers stored with the edges of this arrangement from the curves of a source arrangement
   * to their duplicates, and registers each edge with the duplicates that induce it.
   * \param cv_map The map from the source curves to their duplicates.
   */
  void _relink_edges(const Curve_map& cv_map) {
    std::vector<Curve_halfedges*> dup_curves;
    for (auto eit = this->edges_begin(); eit != this->edges_end(); ++eit) {
      Halfedge_handle e = eit;
      auto& data = e->curve().data();
      dup_curves.clear();
      for (auto* p_cv : data) {
        Curve_halfedges* dup_c = cv_map.find(p_cv)->second;
        dup_curves.push_back(dup_c);
        dup_c->_insert(e);
      }

      // Replace the curve pointers associated with the edge.
      data.clear();
      for (auto* dup_c : dup_curves) data.insert(dup_c);
    }
  }
  //@}

  /// \name Curve insertion and deletion.
  //@{

  /*! inserts a curve into the arrangement.
   * \param cv The curve to be inserted.
   * \param pl a point-location object.
   * \return A handle to the inserted curve.
   */
  template <typename PointLocation>
  Curve_handle _insert_curve(const Curve_2& cv, const PointLocation& pl) {
    // Insert a data-traits curve, which comprises cv and a pointer to a new extended curve, into the base
    // arrangement. Note that the attached observer takes care of updating the edges' set.
    Base_arr_2& base_arr = *this;
    CGAL::insert(base_arr, Data_curve_2(cv, _new_curve_halfedges(cv)), pl);
    return std::prev(m_curves.end());   // the last curve in the list
  }

  /*! inserts a curve into the arrangement, using the default point-location strategy.
   * \param cv The curve to be inserted.
   * \return A handle to the inserted curve.
   */
  Curve_handle _insert_curve(const Curve_2& cv) {
    Base_arr_2& base_arr = *this;
    CGAL::insert(base_arr, Data_curve_2(cv, _new_curve_halfedges(cv)));
    return std::prev(m_curves.end());   // the last curve in the list
  }

  /*! inserts a range of curves into the arrangement.
   * \param begin An iterator pointing to the first curve in the range.
   * \param end A past-the-end iterator for the last curve in the range.
   */
  template <typename InputIterator>
  void _insert_curves(InputIterator begin, InputIterator end) {
    // Create the data-traits curves, each with a pointer to a new extended curve, and perform an aggregated
    // insertion into the base arrangement.
    std::vector<Data_curve_2> data_curves;
    for (auto it = begin; it != end; ++it) data_curves.emplace_back(*it, _new_curve_halfedges(*it));
    Base_arr_2& base_arr = *this;
    CGAL::insert(base_arr, data_curves.begin(), data_curves.end());
  }

  /*! removes a curve from the arrangement (remove all the edges it induces).
   * \param ch A handle to the curve to be removed.
   * \return The number of removed edges.
   */
  Size _remove_curve(Curve_handle ch) {
    // Go over all edges the given curve induces.
    Curve_halfedges* p_cv = &(*ch);
    Size n_removed = 0;
    for (auto it = ch->begin(); it != ch->end();) {
      // Note that we increment the iterator now, as the edge may be removed.
      Halfedge_handle he = *it++;
      if (he->curve().data().size() == 1) {
        // The edge is induced only by our curve; remove it.
        CGAL_assertion(he->curve().data().front() == p_cv);
        Base_arr_2::remove_edge(he);
        ++n_removed;
      }
      // The edge is induced by other curves as well, so we just remove the pointer to our curve from its data.
      else he->curve().data().erase(p_cv);
    }

    _delete_curve_halfedges(p_cv);
    return n_removed;
  }

public:
  /*! sets our arrangement to be the overlay of the two given arrangements.
   * \param arr1 The first arrangement.
   * \param arr2 The second arrangement.
   * \param overlay_tr An overlay-traits class.
   */
  template <typename TopolTraits1, typename TopolTraits2, typename OverlayTraits>
  void _overlay(const Arrangement_on_surface_with_history_2<Geometry_traits_2, TopolTraits1>& arr1,
                const Arrangement_on_surface_with_history_2<Geometry_traits_2, TopolTraits2>& arr2,
                OverlayTraits& overlay_tr) {
    // Clear the current contents of the arrangement.
    clear();

    // Detach the observer from the arrangement, as we do not want to update cross-pointers between the halfedges and
    // the curves during the construction of overlay.
    m_observer.detach();

    // Perform overlay of the base arrangements.
    // Note that the base arrangement types of the inputs are accessed through their public alias.
    using Arr_with_hist1 = Arrangement_on_surface_with_history_2<Geometry_traits_2, TopolTraits1>;
    using Arr_with_hist2 = Arrangement_on_surface_with_history_2<Geometry_traits_2, TopolTraits2>;
    const typename Arr_with_hist1::Base_arrangement_2& base_arr1 = arr1;
    const typename Arr_with_hist2::Base_arrangement_2& base_arr2 = arr2;
    Base_arr_2& base_res = *this;
    CGAL::overlay(base_arr1, base_arr2, base_res, overlay_tr);

    // Duplicate the curves of both input arrangements, and redirect the curve pointers stored with the edges of the
    // result to the duplicates.
    Curve_map cv_map;
    _duplicate_curves(arr1.curves_begin(), arr1.curves_end(), cv_map);
    _duplicate_curves(arr2.curves_begin(), arr2.curves_end(), cv_map);
    _relink_edges(cv_map);

    // Re-attach the observer to the arrangement.
    m_observer.attach(*this);
  }
  //@}
};

//-----------------------------------------------------------------------------
// Global insertion, removal and overlay functions.
//-----------------------------------------------------------------------------

/*! inserts a curve into the arrangement (incremental insertion).
 * The inserted curve may not necessarily be `x`-monotone and may intersect the
 * existing arrangement.
 * \param arr The arrangement-with-history object.
 * \param c The curve to be inserted.
 * \param pl A point-location object associated with the arrangement.
 */
template <typename GeomTraits, typename TopolTraits, typename PointLocation>
typename Arrangement_on_surface_with_history_2<GeomTraits, TopolTraits>::Curve_handle
insert(Arrangement_on_surface_with_history_2<GeomTraits, TopolTraits>& arr,
       const typename GeomTraits::Curve_2& c, const PointLocation& pl) {
  // Obtain an arrangement accessor and perform the insertion.
  using Arr_with_hist_2 = Arrangement_on_surface_with_history_2<GeomTraits, TopolTraits>;
  Arr_with_history_accessor<Arr_with_hist_2> arr_access(arr);
  return arr_access.insert_curve(c, pl);
}

/*! inserts a curve into the arrangement (incremental insertion).
 * The inserted curve may not necessarily be `x`-monotone and may intersect the
 * existing arrangement. The default "walk" point-location strategy is used
 * for inserting the curve.
 * \param arr The arrangement-with-history object.
 * \param c The curve to be inserted.
 */
template <typename GeomTraits, typename TopolTraits>
typename Arrangement_on_surface_with_history_2<GeomTraits, TopolTraits>::Curve_handle
insert(Arrangement_on_surface_with_history_2<GeomTraits, TopolTraits>& arr, const typename GeomTraits::Curve_2& c) {
  // Obtain an arrangement accessor and perform the insertion.
  using Arr_with_hist_2 = Arrangement_on_surface_with_history_2<GeomTraits, TopolTraits>;
  Arr_with_history_accessor<Arr_with_hist_2> arr_access(arr);
  return arr_access.insert_curve(c);
}


/*! inserts a range of curves into the arrangement (aggregated insertion).
 * The inserted curves may intersect one another and may also intersect the
 * existing arrangement.
 * \param arr The arrangement-with-history object.
 * \param begin An iterator for the first curve in the range.
 * \param end A past-the-end iterator for the curve range.
 * \pre The value type of the iterators must be `Curve_2`.
 */
template <typename GeomTraits, typename TopolTraits, typename InputIterator>
void insert(Arrangement_on_surface_with_history_2<GeomTraits, TopolTraits>& arr,
            InputIterator begin, InputIterator end) {
  // Obtain an arrangement accessor and perform the insertion.
  using Arr_with_hist_2 = Arrangement_on_surface_with_history_2<GeomTraits, TopolTraits>;
  Arr_with_history_accessor<Arr_with_hist_2> arr_access(arr);
  arr_access.insert_curves(begin, end);
}

/*! removes a curve from the arrangement (remove all the edges it induces).
 * \param arr The arrangement-with-history object.
 * \param ch A handle to the curve to be removed.
 * \return The number of removed edges.
 */
template <typename GeomTraits, typename TopolTraits>
typename Arrangement_on_surface_with_history_2<GeomTraits, TopolTraits>::Size
remove_curve(Arrangement_on_surface_with_history_2<GeomTraits, TopolTraits>& arr,
             typename Arrangement_on_surface_with_history_2<GeomTraits, TopolTraits>::Curve_handle ch) {
  // Obtain an arrangement accessor and perform the removal.
  using Arr_with_hist_2 = Arrangement_on_surface_with_history_2<GeomTraits, TopolTraits>;
  Arr_with_history_accessor<Arr_with_hist_2> arr_access(arr);
  return arr_access.remove_curve(ch);
}

/*! computes the overlay of two input arrangements.
 * \param arr1 The first arrangement.
 * \param arr2 The second arrangement.
 * \param res Output: The resulting arrangement.
 * \param ovl_traits An overlay-traits class. As `arr1`, `arr2` and `res` are all
 *               templated with the same arrangement-traits class but with
 *               different \dcel types, the overlay-traits class defines the
 *               various overlay operations of pairs of \dcel features from
 *               `TopolTraits1` and `TopolTraits2` to the resulting `ResTopolTraits`.
 */
template <typename GeomTraits, typename TopolTraits1, typename TopolTraits2, typename ResTopolTraits,
          typename OverlayTraits>
void overlay(const Arrangement_on_surface_with_history_2<GeomTraits, TopolTraits1>& arr1,
             const Arrangement_on_surface_with_history_2<GeomTraits, TopolTraits2>& arr2,
             Arrangement_on_surface_with_history_2<GeomTraits, ResTopolTraits>& res, OverlayTraits& ovl_traits)
{ res._overlay(arr1, arr2, ovl_traits); }

/*! computes the overlay of two input arrangements.
 * \param arr1 The first arrangement.
 * \param arr2 The second arrangement.
 * \param res Output The resulting arrangement.
 */
template <typename GeomTraits, typename TopolTraits1, typename TopolTraits2, typename ResTopolTraits>
void overlay(const Arrangement_on_surface_with_history_2<GeomTraits, TopolTraits1>& arr1,
             const Arrangement_on_surface_with_history_2<GeomTraits, TopolTraits2>& arr2,
             Arrangement_on_surface_with_history_2<GeomTraits, ResTopolTraits>& res) {
  using ArrA = Arrangement_on_surface_with_history_2<GeomTraits, TopolTraits1>;
  using ArrB = Arrangement_on_surface_with_history_2<GeomTraits, TopolTraits2>;
  using ArrRes = Arrangement_on_surface_with_history_2<GeomTraits, ResTopolTraits>;
  _Arr_default_overlay_traits_base<ArrA, ArrB, ArrRes> ovl_traits;
  res._overlay(arr1, arr2, ovl_traits);
}

} // namespace CGAL

// The function definitions can be found under:
#include <CGAL/Arrangement_2/Arrangement_on_surface_with_history_2_impl.h>

#include <CGAL/enable_warnings.h>

#endif
