// Copyright (c) 2007,2009,2010,2011 Tel-Aviv University (Israel).
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
// SPDX-License-Identifier: GPL-3.0+
// 
//
// Author(s)     : Ron Wein         <wein@post.tau.ac.il>
//                 Efi Fogel        <efif@post.tau.ac.il>
//                 Baruch Zukerman  <baruchzu@post.tau.ac.il>

#ifndef CGAL_ARRANGEMENT_ON_SURFACE_WITH_HISTORY_2_H
#define CGAL_ARRANGEMENT_ON_SURFACE_WITH_HISTORY_2_H

#include <CGAL/license/Arrangement_on_surface_2.h>

#include <CGAL/disable_warnings.h>

/*! \file
 * The header file for the Arrangement_on_surface_with_history_2 class.
 */

#include <CGAL/Arrangement_on_surface_2.h>
#include <CGAL/Arr_overlay_2.h>
#include <CGAL/Arr_consolidated_curve_data_traits_2.h>
#include <CGAL/Arr_observer.h>
#include <CGAL/In_place_list.h>
#include <CGAL/Arrangement_2/Arr_with_history_accessor.h>

#include <set>

namespace CGAL {

/*! \class
 * A class representing planar subdivisions induced by a set of arbitrary
 * input planar curves. These curves are split and form a set of sweepable
 * (x-monotone) and pairwise interior-disjoint curves that are associated
 * with the arrangement edges.
 * The Arrangement_on_surface_with_history_2 class enables tracking the input
 * curve(s) that originated each such subcurve. It also enables keeping track
 * of the edges that resulted from each input curve. 
 * The Traits parameter corresponds to a traits class that defines the
 * Point_2, X_monotone_curve_2 and Curve_2 types and implements the geometric
 * predicates and constructions for the family of curves it defines.
 * The Dcel parameter should be a model of the ArrDcel concept and support
 * the basic topological operations on a doubly-connected edge-list.
 */
template <class GeomTraits_, class TopTraits_> 
class Arrangement_on_surface_with_history_2 :
  public Arrangement_on_surface_2 
  <Arr_consolidated_curve_data_traits_2
     <GeomTraits_, typename GeomTraits_::Curve_2 *>,
   typename TopTraits_::template rebind
     <Arr_consolidated_curve_data_traits_2
       <GeomTraits_,
        typename GeomTraits_::Curve_2 *>,
      typename TopTraits_::Dcel::template rebind
      <Arr_consolidated_curve_data_traits_2
       <GeomTraits_,
        typename GeomTraits_::Curve_2 *> >::other>::other>
{
public:
  typedef GeomTraits_                                     Geometry_traits_2;
  typedef TopTraits_                                      Base_topology_traits;

private:
  typedef Arrangement_on_surface_with_history_2<Geometry_traits_2,
                                                Base_topology_traits>  Self;

public:
  typedef typename Geometry_traits_2::Point_2             Point_2;
  typedef typename Geometry_traits_2::Curve_2             Curve_2;
  typedef typename Geometry_traits_2::X_monotone_curve_2  X_monotone_curve_2;

  typedef Arr_observer<Self>                              Observer;

protected:
  friend class Arr_observer<Self>;
  friend class Arr_accessor<Self>;
  friend class Arr_with_history_accessor<Self>;

  // Define the data-traits class based on Geometry_traits_2.
  typedef Arr_consolidated_curve_data_traits_2<Geometry_traits_2,Curve_2*>
                                                     Data_traits_2;
  typedef typename Data_traits_2::Curve_2            Data_curve_2;
  typedef typename Data_traits_2::X_monotone_curve_2 Data_x_curve_2;
  typedef typename Data_traits_2::Data_iterator      Data_iterator;

  // Rebind the DCEL and the topology traits to the data-traits class.
  typedef typename Base_topology_traits::Dcel        Base_dcel;
  typedef typename Base_dcel::template 
                 rebind<Data_traits_2>               Dcel_rebind;
  typedef typename Dcel_rebind::other                Data_dcel;

  typedef typename Base_topology_traits::template 
                 rebind<Data_traits_2, Data_dcel>    Top_traits_rebind;
  typedef typename Top_traits_rebind::other          Data_top_traits;

  // The arrangement with history is based on the representation of an
  // arrangement, templated by the data-traits class and the rebound DCEL.
  typedef Arrangement_on_surface_2<Data_traits_2,
                                   Data_top_traits>    Base_arr_2;

public:
  typedef Arr_traits_adaptor_2<Data_traits_2>          Traits_adaptor_2;
 
public:

  typedef Data_top_traits                              Topology_traits;
  typedef Base_arr_2                                   Base_arrangement_2;

  // Types inherited from the base arrangement class:
  typedef typename Base_arr_2::Size                    Size;

  typedef typename Base_arr_2::Vertex                  Vertex;
  typedef typename Base_arr_2::Halfedge                Halfedge;
  typedef typename Base_arr_2::Face                    Face;

  typedef typename Base_arr_2::Vertex_iterator         Vertex_iterator;
  typedef typename Base_arr_2::Vertex_const_iterator   Vertex_const_iterator;
  typedef typename Base_arr_2::Halfedge_iterator       Halfedge_iterator;
  typedef typename Base_arr_2::Halfedge_const_iterator Halfedge_const_iterator;
  typedef typename Base_arr_2::Edge_iterator           Edge_iterator;
  typedef typename Base_arr_2::Edge_const_iterator     Edge_const_iterator;
  typedef typename Base_arr_2::Face_iterator           Face_iterator;
  typedef typename Base_arr_2::Face_const_iterator     Face_const_iterator;
  typedef typename Base_arr_2::Halfedge_around_vertex_circulator
                                      Halfedge_around_vertex_circulator;
  typedef typename Base_arr_2::Halfedge_around_vertex_const_circulator 
                                      Halfedge_around_vertex_const_circulator;
  typedef typename Base_arr_2::Ccb_halfedge_circulator
                                      Ccb_halfedge_circulator;
  typedef typename Base_arr_2::Ccb_halfedge_const_circulator
                                      Ccb_halfedge_const_circulator;
  typedef typename Base_arr_2::Outer_ccb_iterator       
                                      Outer_ccb_iterator;
  typedef typename Base_arr_2::Outer_ccb_const_iterator
                                      Outer_ccb_const_iterator;
  typedef typename Base_arr_2::Inner_ccb_iterator       
                                      Inner_ccb_iterator;
  typedef typename Base_arr_2::Inner_ccb_const_iterator
                                      Inner_ccb_const_iterator;
  typedef typename Base_arr_2::Isolated_vertex_iterator
                                      Isolated_vertex_iterator;
  typedef typename Base_arr_2::Isolated_vertex_const_iterator
                                      Isolated_vertex_const_iterator;

  typedef typename Base_arr_2::Vertex_handle           Vertex_handle;
  typedef typename Base_arr_2::Vertex_const_handle     Vertex_const_handle;
  typedef typename Base_arr_2::Halfedge_handle         Halfedge_handle;
  typedef typename Base_arr_2::Halfedge_const_handle   Halfedge_const_handle;
  typedef typename Base_arr_2::Face_handle             Face_handle;
  typedef typename Base_arr_2::Face_const_handle       Face_const_handle;

protected:
 
  /*! \struct
   * Less functor for comparing two halfedge handles.
   */
  struct Less_halfedge_handle {
    bool operator() (Halfedge_handle h1, Halfedge_handle h2) const
    {
      return (&(*h1) < &(*h2));
    }
  };

  // Forward declaration:
  class Curve_halfedges_observer;

  /*! \class
   * Extension of a curve with the set of edges that it induces.
   * Each edge is represented by one of the halfedges.
   */
  class Curve_halfedges : public Curve_2,
                          public In_place_list_base<Curve_halfedges>
  {
    friend class Curve_halfedges_observer;
    friend class Arrangement_on_surface_with_history_2<Geometry_traits_2,
                                                       Base_topology_traits>;
    friend class Arr_with_history_accessor<
      Arrangement_on_surface_with_history_2<Geometry_traits_2,
                                            Base_topology_traits> >;

  private:
    typedef std::set<Halfedge_handle, Less_halfedge_handle>  Halfedges_set;

    // Data members:
    Halfedges_set m_halfedges;

  public:
    /*! Default constructor. */
    Curve_halfedges ()
    {}

    /*! Constructor from a given curve. */
    Curve_halfedges (const Curve_2& curve) :
      Curve_2(curve)
    {}

    typedef typename Halfedges_set::iterator             iterator;
    typedef typename Halfedges_set::const_iterator       const_iterator;

  private:

    /*! Get the number of edges induced by the curve. */
    Size _size () const
    {
      return (m_halfedges.size());
    }

    /*! Get an iterator for the first edge in the set (const version). */
    const_iterator _begin () const
    {
      return m_halfedges.begin();
    }

    /*! Get an iterator for the first edge in the set (non-const version). */
    iterator _begin ()
    {
      return m_halfedges.begin();
    }
    
    /*! Get a past-the-end iterator for the set edges (const version). */
    const_iterator _end () const
    {
      return m_halfedges.end();
    }

    /*! Get a past-the-end iterator for the set edges (non-const version). */
    iterator _end ()
    {
      return m_halfedges.end();
    }

    /*! Insert an edge to the set. */
    iterator _insert (Halfedge_handle he)
    {
      std::pair<iterator, bool> res = m_halfedges.insert(he);
      CGAL_assertion(res.second);
      return (res.first);
    }
    
    /*! Erase an edge, given by its position, from the set. */
    void _erase(iterator pos)
    {
      m_halfedges.erase(pos);
      return;
    }
    
    /*! Erase an edge from the set. */
    void _erase (Halfedge_handle he)
    {
      size_t res = m_halfedges.erase(he);
      if (res == 0)
        res = m_halfedges.erase(he->twin());
      CGAL_assertion(res != 0);
      return;
    }

    /*! Cleat the edges set. */
    void _clear ()
    {
      m_halfedges.clear();
      return;
    }
  };

  typedef CGAL_ALLOCATOR(Curve_halfedges)               Curves_alloc;
  typedef In_place_list<Curve_halfedges, false>         Curve_halfedges_list;

  /*! \class
   * Observer for the base arrangement. It keeps track of all local changes
   * involving edges and updates the list of halfedges associated with the
   * input curves accordingly.
   */
  class Curve_halfedges_observer : public Arr_observer<Base_arr_2> {
  public:

    typedef typename Base_arr_2::Halfedge_handle     Halfedge_handle;
    typedef typename Base_arr_2::Vertex_handle       Vertex_handle;
    typedef typename Base_arr_2::X_monotone_curve_2  X_monotone_curve_2;

    /*!
     * Notification after the creation of a new edge.
     * \param e A handle to one of the twin halfedges that were created.
     */
    virtual void after_create_edge (Halfedge_handle e)
    {
      _register_edge(e);
    }

    /*!
     * Notification before the modification of an existing edge.
     * \param e A handle to one of the twin halfedges to be updated.
     * \param c The x-monotone curve to be associated with the edge.
     */
    virtual void before_modify_edge (Halfedge_handle e,
                                     const X_monotone_curve_2& /* c */)
    {
      _unregister_edge(e);
    }

    /*!
     * Notification after an edge was modified.
     * \param e A handle to one of the twin halfedges that were updated.
     */
    virtual void after_modify_edge (Halfedge_handle e)
    {
      _register_edge(e);
    }
    
    /*!
     * Notification before the splitting of an edge into two.
     * \param e A handle to one of the existing halfedges.
     * \param c1 The x-monotone curve to be associated with the first edge.
     * \param c2 The x-monotone curve to be associated with the second edge.
     */
    virtual void before_split_edge (Halfedge_handle e,
				    Vertex_handle /* v */,
                                    const X_monotone_curve_2& /* c1 */,
                                    const X_monotone_curve_2& /* c2 */)
    {
      _unregister_edge(e);
    }

    /*!
     * Notification after an edge was split.
     * \param e1 A handle to one of the twin halfedges forming the first edge.
     * \param e2 A handle to one of the twin halfedges forming the second edge.
     */
    virtual void after_split_edge (Halfedge_handle e1, Halfedge_handle e2)
    {
      _register_edge(e1);
      _register_edge(e2);
    }

    /*!
     * Notification before the merging of two edges.
     * \param e1 A handle to one of the halfedges forming the first edge.
     * \param e2 A handle to one of the halfedges forming the second edge.
     * \param c The x-monotone curve to be associated with the merged edge.
     */
    virtual void before_merge_edge (Halfedge_handle e1, Halfedge_handle e2,
                                    const X_monotone_curve_2& /* c */)
    {
      _unregister_edge(e1);
      _unregister_edge(e2);
    }

    /*!
     * Notification after an edge was merged.
     * \param e A handle to one of the twin halfedges forming the merged edge.
     */
    virtual void after_merge_edge (Halfedge_handle e)
    {
      _register_edge(e);
    }

    /*!
     * Notification before the removal of an edge.
     * \param e A handle to one of the twin halfedges to be deleted.
     */
    virtual void before_remove_edge (Halfedge_handle e)
    {
      _unregister_edge(e);
    }

  private:
    
    /*!
     * Register the given halfedge in the set(s) associated with its curve.
     */
    void _register_edge (Halfedge_handle e)
    {
      Curve_halfedges  *curve_halfedges;
      Data_iterator     di;

      for (di = e->curve().data().begin(); di != e->curve().data().end(); ++di)
      {
        curve_halfedges = static_cast<Curve_halfedges*>(*di);
        curve_halfedges->_insert(e);
      }
    }

    /*!
     * Unregister the given halfedge from the set(s) associated with its curve.
     */
    void _unregister_edge (Halfedge_handle e)
    {
      Curve_halfedges  *curve_halfedges;
      Data_iterator     di;

      for (di = e->curve().data().begin(); di != e->curve().data().end(); ++di)
      {
        curve_halfedges = static_cast<Curve_halfedges*>(*di);
        curve_halfedges->_erase(e);
      }
    }
  };

  // Data members:
  Curves_alloc              m_curves_alloc;
  Curve_halfedges_list      m_curves;
  Curve_halfedges_observer  m_observer;
  
public:

  typedef typename Curve_halfedges_list::iterator        Curve_iterator;
  typedef typename Curve_halfedges_list::const_iterator  Curve_const_iterator;

  typedef Curve_iterator                                 Curve_handle;
  typedef Curve_const_iterator                           Curve_const_handle;

  /// \name Constructors.
  //@{

  /*! Default constructor. */
  Arrangement_on_surface_with_history_2 ();

  /*! Copy constructor. */
  Arrangement_on_surface_with_history_2 (const Self& arr);

  /*! Constructor given a traits object. */
  Arrangement_on_surface_with_history_2 (const Geometry_traits_2 *tr);
  //@}

  /// \name Assignment functions.
  //@{

  /*! Assignment operator. */
  Self& operator= (const Self& arr);

  /*! Assign an arrangement with history. */
  void assign (const Self& arr);
  //@}

  /// \name Destruction functions.
  //@{

  /*! Destructor. */
  virtual ~Arrangement_on_surface_with_history_2 ();

  /*! Clear the arrangement. */
  virtual void clear ();
  //@}

  /*! Access the geometry-traits object (const version). */
  inline const Geometry_traits_2 * geometry_traits () const
  {
    return (this->m_geom_traits);
  }

  /*! Access the topology-traits object (non-const version). */
  inline Topology_traits * topology_traits ()
  {
    return (&(this->m_topol_traits));
  }
  
  /*! Access the topology-traits object (const version). */
  inline const Topology_traits* topology_traits () const
  {
    return (&(this->m_topol_traits));
  }

  /// \name Traversal of the arrangement curves.
  //@{
  Size number_of_curves () const { return (m_curves.size()); }

  Curve_iterator curves_begin () { return (m_curves.begin()); }

  Curve_iterator curves_end () { return (m_curves.end()); }

  Curve_const_iterator curves_begin () const { return (m_curves.begin()); }

  Curve_const_iterator curves_end () const { return (m_curves.end()); }
  //@}

  /*! \class
   * Edges iterator - defined as a derived class to make it assignable
   * to the halfedge iterator type.
   */
  class Originating_curve_iterator :
    public I_Dereference_iterator<Data_iterator, Curve_2,
                                  typename Data_iterator::difference_type,
                                  typename Data_iterator::iterator_category>
  {
    typedef I_Dereference_iterator<Data_iterator, Curve_2,
                                   typename Data_iterator::difference_type,
                                   typename Data_iterator::iterator_category>
                                                                          Base;

  public:
 
    Originating_curve_iterator () {}

    Originating_curve_iterator (Data_iterator iter) : Base (iter) {}

    // Casting to a curve iterator.
    operator Curve_iterator () const
    {
      Curve_halfedges   *p_cv = static_cast<Curve_halfedges*>(this->ptr());
      return (Curve_iterator (p_cv));
    }

    operator Curve_const_iterator () const
    {
      const Curve_halfedges   
                        *p_cv = static_cast<Curve_halfedges*>(this->ptr());
      return (Curve_const_iterator  (p_cv));
    }    
  };

  /// \name Traversal of the origin curves of an edge.
  //@{
  Size number_of_originating_curves (Halfedge_const_handle e) const
  {
    return (e->curve().data().size());
  }

  Originating_curve_iterator
  originating_curves_begin (Halfedge_const_handle e) const
  {
    return Originating_curve_iterator (e->curve().data().begin());
  }

  Originating_curve_iterator
  originating_curves_end (Halfedge_const_handle e) const
  {
    return Originating_curve_iterator (e->curve().data().end());
  }
  //@}

  typedef typename Curve_halfedges::const_iterator   Induced_edge_iterator;

  /// \name Traversal of the edges induced by a curve.
  //@{
  Size number_of_induced_edges (Curve_const_handle c) const
  {
    return (c->_size());
  }

  Induced_edge_iterator
  induced_edges_begin (Curve_const_handle c) const { return (c->_begin()); }

  Induced_edge_iterator
  induced_edges_end (Curve_const_handle c) const { return (c->_end()); }
  //@}

  /// \name Manipulating edges.
  //@{

  /*!
   * Split a given edge into two at the given split point.
   * \param e The edge to split (one of the pair of twin halfedges).
   * \param p The split point.
   * \pre p lies in the interior of the curve associated with e.
   * \return A handle for the halfedge whose source is the source of the the
   *         original halfedge e, and whose target is the split point.
   */
  Halfedge_handle split_edge (Halfedge_handle e, const Point_2& p);

  /*!
   * Merge two edges to form a single edge.
   * \param e1 The first edge to merge (one of the pair of twin halfedges).
   * \param e2 The second edge to merge (one of the pair of twin halfedges).
   * \pre e1 and e2 must have a common end-vertex of degree 2 and must
   *      be mergeable.
   * \return A handle for the merged halfedge.
   */
  Halfedge_handle merge_edge (Halfedge_handle e1, Halfedge_handle e2);

  /*!
   * Check if two edges can be merged to a single edge.
   * \param e1 The first edge (one of the pair of twin halfedges).
   * \param e2 The second edge (one of the pair of twin halfedges).
   * \return true iff e1 and e2 are mergeable.
   */
  bool are_mergeable (Halfedge_const_handle e1, Halfedge_const_handle e2) const;

protected:

  /// \name Managing and notifying the arrangement observers.
  //@{

  /*!
   * Register a new observer (so it starts receiving notifications).
   * \param p_obs A pointer to the observer object.
   */
  void _register_observer (Observer *p_obs);

  /*!
   * Unregister an observer (so it stops receiving notifications).
   * \param p_obs A pointer to the observer object.
   * \return Whether the observer was successfully unregistered.
   */
  bool _unregister_observer (Observer *p_obs);
  //@}

  /// \name Curve insertion and deletion.
  //@{

  /*!
   * Insert a curve into the arrangement.
   * \param cv The curve to be inserted.
   * \param pl a point-location object.
   * \return A handle to the inserted curve.
   */
  template <class PointLocation>
  Curve_handle _insert_curve (const Curve_2& cv, const PointLocation& pl)
  {
    // Allocate an extended curve (with an initially empty set of edges)
    // and store it in the curves' list.
    Curve_halfedges   *p_cv = m_curves_alloc.allocate (1);
 #ifdef CGAL_CXX11
    std::allocator_traits<Curves_alloc>::construct(m_curves_alloc, p_cv, cv);
#else
    m_curves_alloc.construct (p_cv, cv);
#endif
    m_curves.push_back (*p_cv);

    // Create a data-traits Curve_2 object, which is comprised of cv and
    // a pointer to the extended curve we have just created.
    // Insert this curve into the base arrangement. Note that the attached
    // observer will take care of updating the edges' set.
    Data_curve_2       data_curve (cv, p_cv);
    Base_arr_2&        base_arr = *this;

    CGAL::insert (base_arr, data_curve, pl);
    
    // Return a handle to the inserted curve (the last in the list).
    Curve_handle       ch = m_curves.end();
    return (--ch);
  }

  /*!
   * Insert a curve into the arrangement, using the default point-location
   * strategy.
   * \param cv The curve to be inserted.
   * \return A handle to the inserted curve.
   */
  Curve_handle _insert_curve (const Curve_2& cv)
  {
    // Allocate an extended curve (with an initially empty set of edges)
    // and store it in the curves' list.
    Curve_halfedges   *p_cv = m_curves_alloc.allocate (1);
    
#ifdef CGAL_CXX11
    std::allocator_traits<Curves_alloc>::construct(m_curves_alloc, p_cv, cv);
#else
    m_curves_alloc.construct (p_cv, cv);
#endif
    m_curves.push_back (*p_cv);

    // Create a data-traits Curve_2 object, which is comprised of cv and
    // a pointer to the extended curve we have just created.
    // Insert this curve into the base arrangement. Note that the attached
    // observer will take care of updating the edges' set.
    Data_curve_2       data_curve (cv, p_cv);
    Base_arr_2&        base_arr = *this;

    CGAL::insert (base_arr, data_curve);
    
    // Return a handle to the inserted curve (the last in the list).
    Curve_handle       ch = m_curves.end();
    return (--ch);
  }
  
  /*!
   * Insert a range of curves into the arrangement.
   * \param begin An iterator pointing to the first curve in the range.
   * \param end A past-the-end iterator for the last curve in the range.
   */
  template <class InputIterator>
  void _insert_curves (InputIterator begin, InputIterator end)
  {
    // Create a list of extended curves (with an initially empty sets of edges)
    std::list<Data_curve_2>   data_curves;

    while (begin != end) {
      Curve_halfedges   *p_cv = m_curves_alloc.allocate (1);
    
#ifdef CGAL_CXX11
      std::allocator_traits<Curves_alloc>::construct(m_curves_alloc, p_cv, *begin);
#else
      m_curves_alloc.construct (p_cv, *begin);
#endif
      m_curves.push_back (*p_cv);

      data_curves.push_back (Data_curve_2 (*begin, p_cv));

      ++begin;
    }

    // Perform an aggregated insertion operation into the base arrangement.
    Base_arr_2&        base_arr = *this;
    CGAL::insert (base_arr, data_curves.begin(), data_curves.end());
  }

  /*!
   * Remove a curve from the arrangement (remove all the edges it induces).
   * \param ch A handle to the curve to be removed.
   * \return The number of removed edges.
   */
  Size _remove_curve (Curve_handle ch)
  {
    // Go over all edges the given curve induces.
    Curve_halfedges                           *p_cv = &(*ch); 
    typename Curve_halfedges::const_iterator   it = ch->_begin();
    Halfedge_handle                            he;
    Size                                       n_removed = 0;

    while (it != ch->_end()) {
      // Check how many curves have originated the current edge.
      // Note we increment the iterator now, as the edge may be removed.
      he = *it;
      ++it;

      if (he->curve().data().size() == 1) {
        // The edge is induced only by out curve - remove it.
        CGAL_assertion (he->curve().data().front() == p_cv);

        Base_arr_2::remove_edge (he);
        n_removed++;
      }
      else {
        // The edge is induced by other curves as well, so we just remove
        // the pointer to out curve from its data container.
        he->curve().data().erase (p_cv);
      }
    }

    // Remove the extended curve object from the list and de-allocate it.
    m_curves.erase (p_cv);

#ifdef CGAL_CXX11
    std::allocator_traits<Curves_alloc>::destroy(m_curves_alloc, p_cv);
#else
    m_curves_alloc.destroy (p_cv);
#endif
    m_curves_alloc.deallocate (p_cv, 1);

    return (n_removed);
  }

public:

  /*!
   * Set our arrangement to be the overlay the two given arrangements.
   * \param arr1 The first arrangement.
   * \param arr2 The second arrangement.
   * \param overlay_tr An overlay-traits class.
   */
  template <class TopTraits1, class TopTraits2, class OverlayTraits>
  void _overlay(const Arrangement_on_surface_with_history_2
                <Geometry_traits_2, TopTraits1>& arr1,
                const Arrangement_on_surface_with_history_2
                <Geometry_traits_2, TopTraits2>& arr2,
                OverlayTraits& overlay_tr)
  {
    // Clear the current contents of the arrangement.
    clear();

    // Detach the observer from the arrangement, as we do not want to update
    // cross-pointers between the halfedges and the curves during the
    // construction of overlay.
    m_observer.detach();

    // Perform overlay of the base arrnagements.
    typedef Arrangement_on_surface_with_history_2<Geometry_traits_2,
                                                  TopTraits1>  Arr_with_hist1;
    typedef typename Arr_with_hist1::Base_arr_2                T_base_arr1;
    typedef Arrangement_on_surface_with_history_2<Geometry_traits_2,
                                                  TopTraits2>  Arr_with_hist2;
    typedef typename Arr_with_hist2::Base_arr_2                T_base_arr2;

    const T_base_arr1&   base_arr1 = arr1;
    const T_base_arr2&   base_arr2 = arr2;
    Base_arr_2&          base_res = *this;

    CGAL::overlay (base_arr1, base_arr2, base_res, overlay_tr);

    // Create duplicates of the curves stored in both input arrangements
    // and map the curves of the original arrangement to their
    // corresponding duplicates.
    typedef std::map<const Curve_halfedges*,
                     Curve_halfedges*>             Curve_map;
    typedef typename Curve_map::value_type         Curve_map_entry;

    Curve_map                 cv_map;
    const Curve_2            *p_cv;
    Curve_halfedges          *dup_c;

    // Duplicate the curves from the first arrangement.
    typename Arr_with_hist1::Curve_const_iterator      ocit1;

    for (ocit1 = arr1.curves_begin(); ocit1 != arr1.curves_end(); ++ocit1) {
      // Create a duplicate of the current curve.
      dup_c = m_curves_alloc.allocate (1);
    
      p_cv = &(*ocit1);
      
#ifdef CGAL_CXX11
      std::allocator_traits<Curves_alloc>::construct(m_curves_alloc, dup_c, *p_cv);
#else
      m_curves_alloc.construct (dup_c, *p_cv);
#endif
      m_curves.push_back (*dup_c);

      // Assign a map entry.
      cv_map.insert (Curve_map_entry (&(*ocit1), dup_c));
    }

    // Duplicate the curves from the second arrangement.
    typename Arr_with_hist2::Curve_const_iterator      ocit2;

    for (ocit2 = arr2.curves_begin(); ocit2 != arr2.curves_end(); ++ocit2) {
      // Create a duplicate of the current curve.
      dup_c = m_curves_alloc.allocate (1);
    
      p_cv = &(*ocit2);
#ifdef CGAL_CXX11
      std::allocator_traits<Curves_alloc>::construct(m_curves_alloc, dup_c, *p_cv);
#else
      m_curves_alloc.construct (dup_c, *p_cv);
#endif
      m_curves.push_back (*dup_c);

      // Assign a map entry.
      cv_map.insert (Curve_map_entry (&(*ocit2), dup_c));
    }

    // Go over the list of halfedges in our arrangement. The curves associated
    // with these edges store pointers to the curves in the original
    // arrangement, so we now have to modify these pointers, according to the
    // mapping we have just created. While doing so, we also construct the set
    // of edges associated with each (duplicated) curve in our arrangement.
    Data_iterator                           dit;
    std::list<Curve_2*>                     dup_curves;
    typename std::list<Curve_2*>::iterator  iter;
    Edge_iterator                           eit;
    Halfedge_handle                         e;
    const Curve_halfedges                  *org_c;

    for (eit = this->edges_begin(); eit != this->edges_end(); ++eit) {
      e = eit;
      dup_curves.clear();
      for (dit = e->curve().data().begin();
           dit != e->curve().data().end(); ++dit)
      {
        org_c = static_cast<Curve_halfedges*>(*dit);
        dup_c = (cv_map.find (org_c))->second;

        dup_curves.push_back (dup_c);
        dup_c->_insert (e);
      }

      // Replace the curve pointers associated with the edge.
      e->curve().data().clear();
      for (iter = dup_curves.begin(); iter != dup_curves.end(); ++iter)
        e->curve().data().insert (*iter);
    }

    // Re-attach the observer to the arrangement.
    m_observer.attach (*this);

    return;
  }
  //@}
};

//-----------------------------------------------------------------------------
// Global insertion, removal and overlay functions.
//-----------------------------------------------------------------------------

/*!
 * Insert a curve into the arrangement (incremental insertion).
 * The inserted curve may not necessarily be x-monotone and may intersect the
 * existing arrangement.
 * \param arr The arrangement-with-history object.
 * \param cv The curve to be inserted.
 * \param pl A point-location object associated with the arrangement.
 */
template <class GeomTraits, class TopTraits, class PointLocation>
typename Arrangement_on_surface_with_history_2<GeomTraits,
                                               TopTraits>::Curve_handle
insert (Arrangement_on_surface_with_history_2<GeomTraits,TopTraits>& arr,
        const typename GeomTraits::Curve_2& c,
        const PointLocation& pl)
{
  // Obtain an arrangement accessor and perform the insertion.
  typedef Arrangement_on_surface_with_history_2<GeomTraits,
                                                TopTraits>     Arr_with_hist_2;
  Arr_with_history_accessor<Arr_with_hist_2>   arr_access (arr);

  return (arr_access.insert_curve (c, pl));
}

/*!
 * Insert a curve into the arrangement (incremental insertion).
 * The inserted curve may not necessarily be x-monotone and may intersect the
 * existing arrangement. The default "walk" point-location strategy is used
 * for inserting the curve.
 * \param arr The arrangement-with-history object.
 * \param cv The curve to be inserted.
 */
template <class GeomTraits, class TopTraits>
typename Arrangement_on_surface_with_history_2<GeomTraits,
                                               TopTraits>::Curve_handle
insert (Arrangement_on_surface_with_history_2<GeomTraits,TopTraits>& arr,
        const typename GeomTraits::Curve_2& c)
{
  // Obtain an arrangement accessor and perform the insertion.
  typedef Arrangement_on_surface_with_history_2<GeomTraits,
                                                TopTraits>     Arr_with_hist_2;
  Arr_with_history_accessor<Arr_with_hist_2>   arr_access (arr);

  return (arr_access.insert_curve (c));
}


/*!
 * Insert a range of curves into the arrangement (aggregated insertion). 
 * The inserted curves may intersect one another and may also intersect the 
 * existing arrangement.
 * \param arr The arrangement-with-history object.
 * \param begin An iterator for the first curve in the range.
 * \param end A past-the-end iterator for the curve range.
 * \pre The value type of the iterators must be Curve_2.
 */
template <class GeomTraits, class TopTraits, class InputIterator>
void insert (Arrangement_on_surface_with_history_2<GeomTraits, TopTraits>& arr,
             InputIterator begin, InputIterator end)
{
  // Obtain an arrangement accessor and perform the insertion.
  typedef Arrangement_on_surface_with_history_2<GeomTraits,
                                                TopTraits>     Arr_with_hist_2;
  Arr_with_history_accessor<Arr_with_hist_2>   arr_access (arr);
  arr_access.insert_curves (begin, end);
}

/*!
 * Remove a curve from the arrangement (remove all the edges it induces).
 * \param ch A handle to the curve to be removed.
 * \return The number of removed edges.
 */
template <class GeomTraits, class TopTraits>
typename Arrangement_on_surface_with_history_2<GeomTraits,
                                               TopTraits>::Size
remove_curve(Arrangement_on_surface_with_history_2<GeomTraits, TopTraits>& arr,
             typename Arrangement_on_surface_with_history_2
               <GeomTraits, TopTraits>::Curve_handle ch)
{
  // Obtain an arrangement accessor and perform the removal.
  typedef Arrangement_on_surface_with_history_2<GeomTraits,
                                                TopTraits>     Arr_with_hist_2;
  Arr_with_history_accessor<Arr_with_hist_2>   arr_access (arr);
  return (arr_access.remove_curve (ch));
}

/*!
 * Compute the overlay of two input arrangement.
 * \param arr1 The first arrangement.
 * \param arr2 The second arrangement.
 * \param res Output: The resulting arrangement.
 * \param traits An overlay-traits class. As arr1, arr2 and res are all
 *               templated with the same arrangement-traits class but with
 *               different DCELs, the overlay-traits class defines the
 *               various overlay operations of pairs of DCEL features from
 *               TopTraits1 and TopTraits2 to the resulting ResTopTraits.
 */
template <class GeomTraits,
          class TopTraits1, class TopTraits2,
          class ResTopTraits, class OverlayTraits>
void
overlay (const Arrangement_on_surface_with_history_2<GeomTraits, TopTraits1>&
         arr1,
         const Arrangement_on_surface_with_history_2<GeomTraits, TopTraits2>&
         arr2,
         Arrangement_on_surface_with_history_2<GeomTraits, ResTopTraits>& res,
         OverlayTraits& ovl_traits)
{
  res._overlay (arr1, arr2, ovl_traits);
}

/*!
 * Compute the overlay of two input arrangement.
 * \param arr1 The first arrangement.
 * \param arr2 The second arrangement.
 * \param res Output: The resulting arrangement.
 */
template <class GeomTraits,
          class TopTraits1, class TopTraits2, class ResTopTraits>
void
overlay (const Arrangement_on_surface_with_history_2<GeomTraits, TopTraits1>&
         arr1,
         const Arrangement_on_surface_with_history_2<GeomTraits, TopTraits2>&
         arr2,
         Arrangement_on_surface_with_history_2<GeomTraits, ResTopTraits>& res)
{
  typedef Arrangement_on_surface_with_history_2<GeomTraits, TopTraits1> ArrA;
  typedef Arrangement_on_surface_with_history_2<GeomTraits, TopTraits2> ArrB;
  typedef Arrangement_on_surface_with_history_2<GeomTraits, ResTopTraits> ArrRes;
  _Arr_default_overlay_traits_base<ArrA, ArrB, ArrRes> ovl_traits;
  res._overlay (arr1, arr2, ovl_traits);
}

} //namespace CGAL

// The function definitions can be found under:
#include <CGAL/Arrangement_2/Arr_on_surface_with_history_2_impl.h>

#include <CGAL/enable_warnings.h>

#endif
