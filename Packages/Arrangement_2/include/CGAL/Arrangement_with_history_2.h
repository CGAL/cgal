// Copyright (c) 2005  Tel-Aviv University (Israel).
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
// Author(s)     : Ron Wein          <wein@post.tau.ac.il>
//                 Efi Fogel         <efif@post.tau.ac.il>

#ifndef CGAL_ARRANGEMENT_WITH_HISTORY_2_H
#define CGAL_ARRANGEMENT_WITH_HISTORY_2_H

/*! \file
 * The header file is for the Arrangement_with_history_2<Traits,Dcel> class.
 */

#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_consolidated_curve_data_traits_2.h>
#include <CGAL/Arr_observer.h>
#include <CGAL/In_place_list.h>

#include <set>

CGAL_BEGIN_NAMESPACE

/*! \class
 * A class representing planar subdivisions induced by a set of arbitrary
 * input planar curves. These curves are split and form a set of x-monotone
 * and pairwise interior-disjoint curves that are associated with the
 * arrangement edges.
 * The Arrangement_with_history_2 class enables tracking the input curve(s)
 * that originated each such subcurve. It also enables keeping track of the
 * edges that resulted from each input curve. 
 * The Traits parameter corresponds to a traits class that defines the
 * Point_2, X_monotone_curve_2 and Curve_2 types and implements the geometric
 * predicates and constructions for the family of curves it defines.
 * The Dcel parameter should be a model of the ArrDcel concept and support
 * the basic topological operations on a doubly-connected edge-list.
 */
template <class Traits_, 
          class Dcel_ = Arr_default_dcel<Traits_> > 
class Arrangement_with_history_2
{
public:

  typedef Traits_                                     Traits_2;
  typedef Dcel_                                       Dcel;
  typedef Arrangement_with_history_2<Traits_2, Dcel>  Self;

  typedef typename Traits_2::Point_2                  Point_2;
  typedef typename Traits_2::Curve_2                  Curve_2;
  typedef typename Traits_2::X_monotone_curve_2       X_monotone_curve_2;

protected:

  typedef Arr_consolidated_curve_data_traits_2<Traits_2,Curve_2*>
                                                      Data_traits;
  typedef typename Data_traits::Curve_2               Data_curve_2;
  typedef typename Data_traits::X_monotone_curve_2    Data_x_curve_2;
  typedef typename Data_x_curve_2::Data_iterator      Data_iterator;
  typedef Arrangement_2<Data_traits,Dcel>             Base_arr_2;
 
public:

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
  typedef typename Base_arr_2::Holes_iterator          Holes_iterator;
  typedef typename Base_arr_2::Holes_const_iterator    Holes_const_iterator;
  typedef typename Base_arr_2::Isolated_vertices_iterator
                                      Isolated_vertices_iterator;
  typedef typename Base_arr_2::Isolated_vertices_const_iterator
                                      Isolated_vertices_const_iterator;

  typedef typename Base_arr_2::Vertex_handle           Vertex_handle;
  typedef typename Base_arr_2::Vertex_const_handle     Vertex_const_handle;
  typedef typename Base_arr_2::Halfedge_handle         Halfedge_handle;
  typedef typename Base_arr_2::Halfedge_const_handle   Halfedge_const_handle;
  typedef typename Base_arr_2::Face_handle             Face_handle;
  typedef typename Base_arr_2::Face_const_handle       Face_const_handle;

protected:
 
  /*! \struct
   * Less functor for comparing tow halfedge handles.
   */
  struct Less_halfedge_handle
  {
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

    /*! Get an iterator for the first edge in the set (non-const version). */
    iterator edges_begin ()
    {
      return (m_halfedges.begin());
    }

    /*! Get an iterator for the first edge in the set (const version). */
    const_iterator edges_begin () const
    {
      return m_halfedges.begin();
    }
    
    /*! Get a past-the-end iterator for the set edges (non-const version). */
    iterator edges_end ()
    {
      return m_halfedges.end();
    }

    /*! Get a past-the-end iterator for the set edges (const version). */
    const_iterator edges_end () const
    {
      return m_halfedges.end();
    }

  private:

    /*! Insert an edge to the set. */
    iterator _insert (Halfedge_handle he)
    {
      std::pair<iterator, bool> res = m_halfedges.insert(he);
      CGAL_assertion(res.second);
      return res.first;
    }
    
    /*! Erase an edge, given by its position, from the set. */
    void _erase(iterator pos)
    {
      m_halfedges.erase(pos);
      return;
    }
    
    /*! Erase an edge from the set. */
    void _erase(Halfedge_handle he)
    {
      size_t res = m_halfedges.erase(he);
      if (res == 0) res = m_halfedges.erase(he->twin());
      CGAL_assertion(res != 0);
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
  class Curve_halfedges_observer : public Arr_observer<Base_arr_2>
  {
  public:

    typedef typename Base_arr_2::Halfedge_handle     Halfedge_handle;
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

      for (di = e->curve().data_begin(); di != e->curve().data_end(); ++di)
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

      for (di = e->curve().data_begin(); di != e->curve().data_end(); ++di)
      {
        curve_halfedges = static_cast<Curve_halfedges*>(*di);
        curve_halfedges->_erase(e);
      }
    }
  };

  // Data members:
  Base_arr_2                m_arr;
  Curves_alloc              m_curves_alloc;
  Curve_halfedges_list      m_curves;
  Curve_halfedges_observer  m_observer;
  
public:

  typedef typename Curve_halfedges_list::iterator        Curve_iterator;
  typedef typename Curve_halfedges_list::const_iterator  Curve_const_iterator;

  typedef Curve_iterator                                 Curve_handle;
  typedef Curve_const_iterator                           Curve_const_handle;

  typedef typename Curve_halfedges::iterator         Curve_edge_iterator;
  typedef typename Curve_halfedges::const_iterator   Curve_edge_const_iterator;

  typedef typename Data_x_curve_2::Data_iterator         Origin_iterator;
  typedef typename Data_x_curve_2::Data_const_iterator   Origin_const_iterator;

  /// \name Constructors.
  //@{

  /*! Default constructor. */
  Arrangement_with_history_2 () :
    m_arr ()
  {
    m_observer.attach (m_arr);
  }

  /*!
   * Copy constructor.
   * \todo not implemented yet.
   */
  Arrangement_with_history_2 (const Self& arr);

  /*! Constructor given a traits object. */
  Arrangement_with_history_2 (Traits_2 *tr) :
    m_arr (static_cast<Data_traits*> (tr))
  {
    m_observer.attach (m_arr);
  }
  //@}

  /// \name Assignment functions.
  //@{

  /*! 
   * Assignment operator.
   * \todo not implemented yet.
   */
  Self& operator= (const Self& arr);

  /*!
   * Assign an arrangement with history.
   * \todo not implemented yet.
   */
  void assign (const Self& arr);
  //@}

  /// \name Destruction functions.
  //@{

  /*! Destructor. */
  virtual ~Arrangement_with_history_2 ()
  {
    clear();
  }

  /*! Clear the arrangement. */
  void clear ()
  {
    // Free all stored curves.
    Curve_iterator         cit = m_curves.begin();
    Curve_halfedges       *p_cv;

    while (cit != m_curves.end())
    {
      p_cv = &(*cit);
      ++cit;

      m_curves.erase (p_cv);
      m_curves_alloc.destroy (p_cv);
      m_curves_alloc.deallocate (p_cv, 1);
    }
    m_curves.destroy();

    // Clear the arrangement.
    m_arr.clear();

    return;
  }
  //@}

  /// \name Access the arrangement dimensions.
  //@{

  /*! Check whether the arrangement is empty. */
  bool is_empty () const
  {
    return (m_arr.is_empty());
  }
   
  /*! Get the number of arrangement vertices. */
  Size number_of_vertices () const
  {
    return (m_arr.number_of_vertices());
  }

  /*! Get the number of arrangement halfedges (the result is always even). */
  Size number_of_halfedges () const
  {
    return (m_arr.number_of_halfedges());
  }

  /*! Get the number of arrangement edges. */
  Size number_of_edges () const
  {
    return (m_arr.number_of_edges());
  }

  /*! Get the number of arrangement faces. */
  Size number_of_faces () const
  {
    return (m_arr.number_of_faces());
  }
  //@}

  /// \name Traversal functions for the arrangement vertices.
  //@{

  /*! Get an iterator for the first vertex in the arrangement. */
  Vertex_iterator vertices_begin() 
  { 
    return (m_arr.vertices_begin());
  }

  /*! Get a past-the-end iterator for the arrangement vertices. */
  Vertex_iterator vertices_end()
  {
    return (m_arr.vertices_end());
  }

  /*! Get a const iterator for the first vertex in the arrangement. */
  Vertex_const_iterator vertices_begin() const
  { 
    return (m_arr.vertices_begin());
  }
  
  /*! Get a past-the-end const iterator for the arrangement vertices. */
  Vertex_const_iterator vertices_end() const
  {
    return (m_arr.vertices_end());
  }
  //@}

  /// \name Traversal functions for the arrangement halfedges.
  //@{

  /*! Get an iterator for the first halfedge in the arrangement. */
  Halfedge_iterator halfedges_begin() 
  { 
    return (m_arr.halfedges_begin());
  }

  /*! Get a past-the-end iterator for the arrangement halfedges. */
  Halfedge_iterator halfedges_end()
  {
    return (m_arr.halfedges_end());
  }

  /*! Get a const iterator for the first halfedge in the arrangement. */
  Halfedge_const_iterator halfedges_begin() const
  { 
    return (m_arr.halfedges_begin());
  }
  
  /*! Get a past-the-end const iterator for the arrangement halfedges. */
  Halfedge_const_iterator halfedges_end() const
  {
    return (m_arr.halfedges_end());
  }
  //@}

  /// \name Traversal functions for the arrangement edges.
  //@{

  /*! Get an iterator for the first edge in the arrangement. */
  Edge_iterator edges_begin() 
  { 
    return (m_arr.edges_begin());
  }

  /*! Get a past-the-end iterator for the arrangement edges. */
  Edge_iterator edges_end()
  {
    return (m_arr.edges_end());
  }

  /*! Get a const iterator for the first edge in the arrangement. */
  Edge_const_iterator edges_begin() const
  { 
    return (m_arr.edges_begin());
  }
  
  /*! Get a past-the-end const iterator for the arrangement edges. */
  Edge_const_iterator edges_end() const

  {
    return (m_arr.edges_end());
  }
  //@}

  /// \name Traversal functions for the arrangement faces.
  //@{

  /*! Get the unbounded face (non-const version). */
  Face_handle unbounded_face ()
  {
    return (m_arr.unbounded_face());
  }

  /*! Get the unbounded face (const version). */
  Face_const_handle unbounded_face () const
  {
    return (m_arr.unbounded_face());
  }

  /*!
   * Get the face containing the given isolated vertex (non-const version).
   * \param The query vertex.
   * \pre v is an isolated vertex (it has no incident halfedges).
   * \return A handle to the face containing v.
   */
  Face_handle incident_face (Vertex_handle v)
  {
    return (m_arr.incident_face (v));
  }

  /*!
   * Get the face containing the given isolated vertex (const version).
   * \param The query vertex.
   * \pre v is an isolated vertex (it has no incident halfedges).
   * \return A const handle to the face containing v.
   */
  inline Face_const_handle incident_face (Vertex_const_handle v) const
  {
    return (m_arr.incident_face (v));
  }

  /*! Get an iterator for the first face in the arrangement. */
  Face_iterator faces_begin() 
  { 
    return (m_arr.faces_begin());
  }

  /*! Get a past-the-end iterator for the arrangement faces. */
  Face_iterator faces_end()
  {
    return (m_arr.faces_end());
  }

  /*! Get a const iterator for the first face in the arrangement. */
  Face_const_iterator faces_begin() const
  { 
    return (m_arr.faces_begin());
  }
  
  /*! Get a past-the-end const iterator for the arrangement faces. */
  Face_const_iterator faces_end() const
  {
    return (m_arr.faces_end());
  }
  //@}

  /// \name Traversal functions for the arrangement curves.
  //@{

  /*! Get an iterator for the first curve in the arrangement. */
  Curve_iterator curves_begin () 
  { 
    return (m_curves.begin()); 
  }

  /*! Get a past-the-end iterator for the arrangement curves. */
  Curve_iterator curves_end ()
  {
    return (m_curves.end()); 
  }

  /*! Get a const iterator for the first curve in the arrangement. */
  Curve_const_iterator curves_begin () const 
  { 
    return (m_curves.begin()); 
  }

  /*! Get a const past-the-end iterator for the arrangement curves. */
  Curve_const_iterator curves_end () const
  {
    return (m_curves.end()); 
  }
  //@}

  /// \name Casting away constness for handle types.
  //@{
  Vertex_handle non_const_handle (Vertex_const_handle vh)
  {
    return (m_arr.non_const_handle (vh));
  }

  Halfedge_handle non_const_handle (Halfedge_const_handle hh)
  {
    return (m_arr.non_const_handle (hh));
  }

  Face_handle non_const_handle (Face_const_handle fh)
  {
    return (m_arr.non_const_handle (fh));
  }
  //@}
  
  /// \name Curve insertion and deletion.
  //@{

  /*!
   * Insert a curve into the arrangement.
   * \param cv The curve to be inserted.
   * \return A handle to the inserted curve.
   */
  Curve_handle insert (const Curve_2& cv)
  {
    // Allocate an extended curve (with an initially empty set of edges)
    // and store it in the curves' list.
    Curve_halfedges   *p_cv = m_curves_alloc.allocate (1);
    
    m_curves_alloc.construct (p_cv, cv);
    m_curves.push_back (*p_cv);

    // Create a data-traits Curve_2 object, which is comprised of cv and
    // a pointer to the extended curve we have just created.
    // Insert this curve into the base arrangement. Note that the attached
    // observer will take care of updating the edges' set.
    Data_curve_2       data_curve (cv, p_cv);

    CGAL::insert (m_arr, data_curve);
    
    // Return a handle to the inserted curve (the last in the list).
    Curve_handle       ch = m_curves.end();
    return (--ch);
  }
  
  /*!
   * Insert a range of curve into the arrangement.
   * \param begin An iterator pointing to the first curve in the range.
   * \param end A past-the-end iterator for the last curve in the range.
   */
  template <class InputIterator>
  void insert (InputIterator begin, InputIterator end)
  {
    // Create a list of extended curves (with an initially empty sets of edges)
    std::list<Data_curve_2>   data_curves;

    while (begin != end)
    {
      Curve_halfedges   *p_cv = m_curves_alloc.allocate (1);
    
      m_curves_alloc.construct (p_cv, *begin);
      m_curves.push_back (*p_cv);

      data_curves.push_back (Data_curve_2 (*begin, p_cv));

      ++begin;
    }

    // Perform an aggregated insertion operation into the base arrangement.
    CGAL::insert (m_arr, data_curves.begin(), data_curves.end());
    return;
  }

  /*!
   * Remove a curve from the arrangement (remove all the edges it induces).
   * \return The number of removed edges.
   */
  size_t remove (Curve_handle& ch)
  {
    // Go over all edges the given curve induces.
    Curve_halfedges                     *p_cv = &(*ch); 
    typename Curve_halfedges::iterator   it;
    Halfedge_handle                      he;
    size_t                               n_removed = 0;

    for (it = ch->begin(); it != ch->end(); ++it)
    {
      // Check how many curves have originated the current edge.
      he = *it;
      if (he->curve().number_of_data_objects() == 1)
      {
	// The edge is induced only by out curve - remove it.
	CGAL_assertion (he->curve().get_data() == p_cv);

	m_arr.remove_edge (he);
	n_removed++;
      }
      else
      {
	// The edge is induced by other curves as well, so we just remove
	// the pointer to out curve from its data container.
	he->curve().remove_data (p_cv);
      }
    }

    // Remove the extended curve object from the list and de-allocate it.
    m_curves.erase (p_cv);
    m_curves_alloc.destroy (p_cv);
    m_curves_alloc.deallocate (p_cv, 1);

    return;
  }
  //@}

};

CGAL_END_NAMESPACE

#endif
