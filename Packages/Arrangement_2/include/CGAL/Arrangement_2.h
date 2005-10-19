// Copyright (c) 2005 Tel-Aviv University (Israel).
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
//                 (based on old version by: Iddo Hanniel,
//                                           Eyal Flato,
//                                           Oren Nechushtan,
//                                           Ester Ezra,
//                                           Shai Hirsch,
//                                           and Eugene Lipovetsky)
#ifndef CGAL_ARRANGEMENT_2_H
#define CGAL_ARRANGEMENT_2_H

/*! \file
 * The header file for the Arrangement_2<Traits,Dcel> class.
 */

#include <CGAL/HalfedgeDS_iterator.h>
#include <CGAL/In_place_list.h>
#include <CGAL/Arr_default_dcel.h>
#include <CGAL/Arr_observer.h>
#include <CGAL/Arrangement_2/Arr_accessor.h>
#include <CGAL/Arrangement_2/Arr_traits_wrapper_2.h>
#include <CGAL/IO/Arrangement_2_reader.h>
#include <map>
#include <vector>
#include <algorithm>

CGAL_BEGIN_NAMESPACE


/*! \class
 * The arrangement class, representing planar subdivisions induced by
 * a set of arbitrary planar curves. 
 * The Traits parameter corresponds to a traits class that defines the
 * Point_2 and X_monotone_curve_2 types and implements the geometric
 * predicates and constructions for the family of curves it defines.
 * The Dcel parameter should be a model of the ArrDcel concept and support
 * the basic topological operations on a doubly-connected edge-list.
 */
template <class Traits_, 
          class Dcel_ = Arr_default_dcel<Traits_> > 
class Arrangement_2
{
public:

  typedef Traits_                               Traits_2;
  typedef Dcel_                                 Dcel;
  typedef Arrangement_2<Traits_2,Dcel>          Self;

  typedef typename Traits_2::Point_2            Point_2;
  typedef typename Traits_2::X_monotone_curve_2 X_monotone_curve_2;

  typedef typename Dcel::Size                   Size;

protected:

  friend class Arr_observer<Self>;
  friend class Arr_accessor<Self>;
  friend class Arrangement_2_reader<Self>;
  
  typedef Arr_traits_basic_wrapper_2<Traits_2>  Traits_wrapper_2;

  // Internal DCEL types:
  typedef typename Dcel::Vertex                      DVertex;
  typedef typename Dcel::Halfedge                    DHalfedge;
  typedef typename Dcel::Face                        DFace;

  typedef typename Dcel::difference_type             DDifference;  

  typedef typename Dcel::iterator_category           DIterator_category;

  typedef typename Dcel::Vertex_iterator             DVertex_iter;
  typedef typename Dcel::Vertex_const_iterator       DVertex_const_iter;

  typedef typename Dcel::Halfedge_iterator           DHalfedge_iter;
  typedef typename Dcel::Halfedge_const_iterator     DHalfedge_const_iter;

  typedef typename Dcel::Edge_iterator               DEdge_iter;
  typedef typename Dcel::Edge_const_iterator         DEdge_const_iter;

  typedef typename Dcel::Face_iterator               DFace_iter;
  typedef typename Dcel::Face_const_iterator         DFace_const_iter;

  typedef typename DFace::Holes_iterator             DHoles_iter;
  typedef typename DFace::Holes_const_iterator       DHoles_const_iter;

  typedef typename DFace::Isolated_vertices_iterator
                                               DIsolated_vertices_iter;
  typedef typename DFace::Isolated_vertices_const_iterator
                                               DIsolated_vertices_const_iter;

public:


  // Forward declerations:
  class Vertex;
  class Halfedge;
  class Face;

  // Definition of the halfedge data-structure itereators and circulators:
  typedef I_HalfedgeDS_iterator
    <DVertex_iter, Vertex, DDifference,
     DIterator_category>                          Vertex_iterator;
  
  typedef I_HalfedgeDS_const_iterator
    <DVertex_const_iter, DVertex_iter, Vertex,
     DDifference, DIterator_category>             Vertex_const_iterator;
  
  typedef I_HalfedgeDS_iterator
    <DHalfedge_iter, Halfedge, DDifference,
     DIterator_category>                          Halfedge_iterator;
  
  typedef I_HalfedgeDS_const_iterator
    <DHalfedge_const_iter, DHalfedge_iter,
     Halfedge, DDifference, DIterator_category>   Halfedge_const_iterator;

  /*! \class
   * Edges iterator - defined as a derived class to make it assignable
   * to the halfedge iterator type.
   */
  class Edge_iterator :
    public I_HalfedgeDS_iterator<DEdge_iter,
                                 Halfedge, DDifference,
                                 DIterator_category>
  {
    typedef I_HalfedgeDS_iterator<DEdge_iter,
                                  Halfedge, DDifference,
                                  DIterator_category>        Base;

  public:
 
    Edge_iterator ()
    {}

    Edge_iterator (DEdge_iter iter) :
      Base (iter)
    {}

    // Casting to a halfedge iterator.
    operator Halfedge_iterator () const
    {
      return (Halfedge_iterator (DHalfedge_iter (this->current_iterator())));
    }

    operator Halfedge_const_iterator () const
    {
      return (Halfedge_const_iterator 
              (DHalfedge_const_iter (this->current_iterator())));
    }    
  };

  class Edge_const_iterator :
    public I_HalfedgeDS_const_iterator<DEdge_const_iter, DEdge_iter,
                                       Halfedge, DDifference,
                                       DIterator_category>
  {
    typedef I_HalfedgeDS_const_iterator<DEdge_const_iter, DEdge_iter,
                                        Halfedge, DDifference,
                                        DIterator_category>           Base;

  public:
 
    Edge_const_iterator ()
    {}

    Edge_const_iterator (Edge_iterator iter) :
      Base (iter)
    {}

    Edge_const_iterator (DEdge_const_iter iter) :
      Base (iter)
    {}

    // Casting to a halfedge iterator.
    operator Halfedge_const_iterator () const
    {
      return (Halfedge_const_iterator 
              (DHalfedge_const_iter (this->current_iterator())));
    }
  };
  
  typedef I_HalfedgeDS_iterator
    <DFace_iter, Face, DDifference,
     DIterator_category>                          Face_iterator;
  
  typedef I_HalfedgeDS_const_iterator
    <DFace_const_iter, DFace_iter, Face,
     DDifference, DIterator_category>             Face_const_iterator;

  typedef _HalfedgeDS_vertex_circ
    <Halfedge, Halfedge_iterator,
     Bidirectional_circulator_tag>    Halfedge_around_vertex_circulator;

  typedef _HalfedgeDS_vertex_const_circ
    <Halfedge,
     Halfedge_const_iterator,
     Bidirectional_circulator_tag>    Halfedge_around_vertex_const_circulator;

  typedef _HalfedgeDS_facet_circ
    <Halfedge,
     Halfedge_iterator,
     Bidirectional_circulator_tag>    Ccb_halfedge_circulator;

  typedef _HalfedgeDS_facet_const_circ
    <Halfedge,
     Halfedge_const_iterator,
     Bidirectional_circulator_tag>    Ccb_halfedge_const_circulator;
  
  typedef I_HalfedgeDS_iterator
    <DHoles_iter, Ccb_halfedge_circulator,
     DDifference,
     DIterator_category>              Holes_iterator;
  
  typedef I_HalfedgeDS_const_iterator
    <DHoles_const_iter, DHoles_iter,
     Ccb_halfedge_const_circulator,
     DDifference,
     DIterator_category>              Holes_const_iterator;

  /*! \class
   * Isolated vertices iterator - defined as a class to make it assignable
   * to the vertex iterator type.
   */
  class Isolated_vertices_iterator :
    public I_HalfedgeDS_iterator<DIsolated_vertices_iter,
                                 Vertex, DDifference,
                                 DIterator_category>
  {
    typedef I_HalfedgeDS_iterator<DIsolated_vertices_iter,
                                  Vertex, DDifference,
                                  DIterator_category>         Base;

  public:

    Isolated_vertices_iterator ()
    {}

    Isolated_vertices_iterator (DIsolated_vertices_iter iter) :
      Base (iter)
    {}

    // Casting to a vertex iterator.
    operator Vertex_iterator () const
    {
      return (Vertex_iterator (DVertex_iter (this->ptr())));
    }

    operator Vertex_const_iterator () const
    {
      return (Vertex_const_iterator (DVertex_const_iter (this->ptr())));
    }
  };
  
  class Isolated_vertices_const_iterator :
    public I_HalfedgeDS_const_iterator <DIsolated_vertices_const_iter,
                                        DIsolated_vertices_iter,
                                        Vertex, DDifference,
                                        DIterator_category>
  {
    typedef I_HalfedgeDS_const_iterator <DIsolated_vertices_const_iter,
                                         DIsolated_vertices_iter,
                                         Vertex, DDifference,
                                         DIterator_category>            Base;

  public:

    Isolated_vertices_const_iterator ()
    {}

    Isolated_vertices_const_iterator (Isolated_vertices_iterator iter) :
      Base (iter)
    {}

    Isolated_vertices_const_iterator (DIsolated_vertices_const_iter iter) :
      Base (iter)
    {}

    // Casting to a vertex iterator.
    operator Vertex_const_iterator () const
    {
      return (Vertex_const_iterator (DVertex_const_iter (this->ptr())));
    }
  };

  // Definition of handles (equivalent to iterators):
  typedef Vertex_iterator              Vertex_handle;
  typedef Halfedge_iterator            Halfedge_handle;
  typedef Face_iterator                Face_handle;

  typedef Vertex_const_iterator        Vertex_const_handle;
  typedef Halfedge_const_iterator      Halfedge_const_handle;
  typedef Face_const_iterator          Face_const_handle;

  /*!
   * \class The arrangement vertex class.
   */
  class Vertex : public DVertex
  {
    typedef DVertex               Base;

  public:

    /*! Default constrcutor. */
    Vertex()
    {}

    /*! Get the vertex degree (number of incident edges). */
    Size degree () const
    {
      if (this->is_isolated())
        return (0);

      // Go around the vertex and count the incident halfedges.
      const DHalfedge  *he_first = Base::halfedge();
      const DHalfedge  *he_curr = he_first;
      Size              n = 0;

      if (he_curr != NULL)
      {
        do
        {
          n++;
          he_curr = he_curr->next()->opposite();
        } while (he_curr != he_first);
      }
      return (n);
    }

    /*!
     * Get the incident halfedges (non-const version).
     * \pre The vertex is not isolated.
     */
    Halfedge_around_vertex_circulator incident_halfedges() 
    {
      CGAL_precondition (! this->is_isolated());

      return Halfedge_around_vertex_circulator
        (DHalfedge_iter (Base::halfedge()));
    }

    /*!
     * Get the incident halfedges (const version).
     * \pre The vertex is not isolated.
     */
    Halfedge_around_vertex_const_circulator incident_halfedges() const 
    {
      CGAL_precondition (! this->is_isolated());

      return Halfedge_around_vertex_const_circulator
        (DHalfedge_const_iter (Base::halfedge())); 
    }

    /*! 
     * Get the face that contains the vertex (non-const version).
     * \pre The vertex is isolated.
     */
    Face_handle face()
    {
      CGAL_precondition (this->is_isolated());

      return (DFace_iter (Base::face()));
    }

    /*! 
     * Get the face that contains the vertex (const version).
     * \pre The vertex is isolated.
     */
    Face_const_handle face() const
    {
      CGAL_precondition (this->is_isolated());

      return (DFace_const_iter (Base::face()));
    }

  private:

    // Blocking access to inherited functions from the Dcel::Vertex.
    void set_point (Point_2* );
    const DHalfedge* halfedge () const;
    DHalfedge* halfedge ();
    void set_halfedge (DHalfedge* );
    void set_face (DFace* );

  };

  /*!
   * \class The arrangement halfedge class.
   */
  class Halfedge : public DHalfedge
  {
    typedef DHalfedge             Base;

  public:

    /*! Default constrcutor. */
    Halfedge ()
    {}

    /*! Get the source vertex (non-const version). */
    Vertex_handle source ()
    {
      return (DVertex_iter (Base::opposite()->vertex()));
    }

    /*! Get the source vertex (const version). */
    Vertex_const_handle source () const
    {
      return (DVertex_const_iter (Base::opposite()->vertex()));
    }
    
    /*! Get the target vertex (non-const version). */
    Vertex_handle target ()
    {
      return (DVertex_iter (Base::vertex()));
    }

    /*! Get the target vertex (const version). */
    Vertex_const_handle target () const
    {
      return (DVertex_const_iter (Base::vertex()));
    }
    
    /*! Get the incident face (non-const version). */
    Face_handle face()
    {
      return (DFace_iter (Base::face()));
    }

    /*! Get the incident face (const version). */
    Face_const_handle face() const
    {
      return (DFace_const_iter (Base::face()));
    }

    /*! Get the twin halfedge (non-const version). */
    Halfedge_handle twin()
    {
      return (DHalfedge_iter (Base::opposite()));
    }

    /*! Get the twin halfedge (const version). */
    Halfedge_const_handle twin() const 
    {
      return (DHalfedge_const_iter (Base::opposite()));
    }

   /*! Get the previous halfegde in the chain (non-const version). */
    Halfedge_handle prev () 
    { 
      return (DHalfedge_iter (Base::prev()));
    }

    /*! Get the previous halfegde in the chain (const version). */
    Halfedge_const_handle prev () const 
    { 
      return (DHalfedge_const_iter (Base::prev()));
    }

    /*! Get the next halfegde in the chain (non-const version). */
    Halfedge_handle next () 
    { 
      return (DHalfedge_iter (Base::next()));
    }

    /*! Get the next halfegde in the chain (const version). */
    Halfedge_const_handle next () const 
    { 
      return (DHalfedge_const_iter (Base::next()));
    }

    
    /*! Get the connected component of the halfedge (non-const version). */
    Ccb_halfedge_circulator ccb ()
    { 
      return Ccb_halfedge_circulator (DHalfedge_iter (this));
    }

    /*! Get the connected component of the halfedge (const version). */
    Ccb_halfedge_const_circulator ccb () const
    { 
      return Ccb_halfedge_const_circulator (DHalfedge_const_iter (this));
    }

  private:

    // Blocking access to inherited functions from the Dcel::Halfedge.
    void set_curve (X_monotone_curve_2* );
    const DHalfedge* opposite () const;
    DHalfedge* opposite ();
    void set_opposite (DHalfedge* );
    void set_direction (Comparison_result );
    void set_prev (DHalfedge* he);
    void set_next (DHalfedge* he);
    const DVertex* vertex () const ;
    DVertex* vertex ();
    void set_vertex (DVertex* v);
    void set_face (DFace* f);
  };

  /*!
   * \class The arrangement face class.
   */
  class Face : public DFace
  {
    typedef DFace                 Base;

  public:

    /*! Default constrcutor. */
    Face()
    {}

    /*! Check whether the face is unbounded. */
    bool is_unbounded () const
    {
      // Check whether the outer-boundary edge exists or not.
      return (Base::halfedge() == NULL);
    }

    /*! Get an iterator for the holes inside the face (non-const version). */
    Holes_iterator holes_begin() 
    {
      return (DHoles_iter (Base::holes_begin()));
    }

    /*! Get an iterator for the holes inside the face (const version). */
    Holes_const_iterator holes_begin() const
    {
      return (DHoles_const_iter (Base::holes_begin()));
    }
    
    /*! Get a past-the-end iterator for the holes (non-const version). */
    Holes_iterator holes_end() 
    {
      return (DHoles_iter (Base::holes_end()));
    }

    /*! Get a past-the-end iterator for the holes (const version). */
    Holes_const_iterator holes_end() const 
    {
      return (DHoles_const_iter (Base::holes_end()));
    }

    /*! Get an iterator for the isolated_vertices inside the face
     * (non-const version). */
    Isolated_vertices_iterator isolated_vertices_begin ()
    {
      return (DIsolated_vertices_iter (Base::isolated_vertices_begin()));
    }

    /*! Get an iterator for the isolated_vertices inside the face
     * (const version). */
    Isolated_vertices_const_iterator isolated_vertices_begin () const
    {
      return (DIsolated_vertices_const_iter (Base::isolated_vertices_begin()));
    }
    
    /*! Get a past-the-end iterator for the isolated_vertices 
     * (non-const version). */
    Isolated_vertices_iterator isolated_vertices_end () 
    {
      return (DIsolated_vertices_iter (Base::isolated_vertices_end()));
    }

    /*! Get a past-the-end iterator for the isolated_vertices
     * (const version). */
    Isolated_vertices_const_iterator isolated_vertices_end () const 
    {
      return (DIsolated_vertices_const_iter (Base::isolated_vertices_end()));
    }

    /*! Get a circulator for the outer boundary (non-const version). */
    Ccb_halfedge_circulator outer_ccb () 
    {
      CGAL_precondition(Base::halfedge() != NULL); 
      return Ccb_halfedge_circulator (DHalfedge_iter (Base::halfedge()));
    }

    /*! Get a circulator for the outer boundary (const version). */
    Ccb_halfedge_const_circulator outer_ccb () const
    {
      CGAL_precondition(Base::halfedge() != NULL); 
      return Ccb_halfedge_const_circulator
        (DHalfedge_const_iter (Base::halfedge()));
    }

  private:

    // Blocking access to inherited functions from the Dcel::Face.
    const DHalfedge* halfedge () const;
    DHalfedge* halfedge ();
    void set_halfedge (DHalfedge* );
    void add_hole (DHalfedge* );
    void erase_hole (DHoles_iter );
    void add_isolated_vertex (DVertex* );
    void erase_isolated_vertex (DIsolated_vertices_iter );
  };

protected:

  /*!
   * \class Representation of a point object stored in the points' container.
   */
  typedef typename Traits_2::Point_2                  Base_point_2;

  class Stored_point_2 : public Base_point_2,
                         public In_place_list_base<Stored_point_2>
  {
  public:

    /*! Default constructor. */
    Stored_point_2 ()
    {}

    /*! Constructor from a point. */
    Stored_point_2 (const Base_point_2& p) :
      Base_point_2 (p)
    {}
  };

  /*!
   * \class Representation of an x-monotone curve object stored in the 
   *        curves' container.
   */  
  class Stored_curve_2 : public X_monotone_curve_2,
                         public In_place_list_base<Stored_curve_2>
  {
  public:

    /*! Default constructor. */
    Stored_curve_2 ()
    {}

    /*! Constructor from an x-monotone curve. */
    Stored_curve_2 (const X_monotone_curve_2& cv) :
      X_monotone_curve_2 (cv)
    {}
  };

  typedef CGAL_ALLOCATOR(Stored_point_2)          Points_alloc;
  typedef CGAL_ALLOCATOR(Stored_curve_2)          Curves_alloc;

  typedef In_place_list<Stored_point_2, false>    Points_container;
  typedef In_place_list<Stored_curve_2, false>    X_curves_container;

  typedef Arr_observer<Self>                      Observer;
  typedef std::list<Observer*>                    Observers_container;
  typedef typename Observers_container::iterator  Observers_iterator;

  typedef typename Observers_container::reverse_iterator  
                                                  Observers_rev_iterator;

  // Data members:
  Dcel                dcel;         // The DCEL representing the arrangement.
  DFace              *un_face;      // The unbounded face of the DCEL.
  Size                n_iso_verts;  // Number of isolated vertices.
  Points_container    points;       // Container for the points that
                                    // correspond to the vertices.
  Points_alloc        points_alloc; // Allocator for the points.

  X_curves_container  curves;       // Container for the x-monotone curves
                                    // that correspond to the edges.
  Curves_alloc        curves_alloc; // Allocator for the curves.
  Observers_container observers;    // Storing pointers to existing observers.
  Traits_wrapper_2   *traits;       // The traits wrapper.
  bool                own_traits;   // Inidicate whether we should evetually
                                    // free the traits object.

public:

  /// \name Constructors.
  //@{

  /*! Default constructor. */
  Arrangement_2 ();

  /*! Copy constructor. */
  Arrangement_2 (const Self& arr);

  /*! Constructor given a traits object. */
  Arrangement_2 (Traits_2 *tr);
  //@}

  /// \name Assignment functions.
  //@{

  /*! Assignment operator. */
  Self& operator= (const Self& arr);

  /*! Assign an arrangement. */
  void assign (const Self& arr);
  //@}

  /// \name Destruction functions.
  //@{

  /*! Destructor. */
  virtual ~Arrangement_2 ();

  /*! Clear the arrangement. */
  virtual void clear();
  //@}

  /*! Access the traits object (non-const version). */
  Traits_2* get_traits ()
  {
    return (traits);
  }
  
  /*! Access the traits object (const version). */
  const Traits_2* get_traits () const
  {
    return (traits);
  }

  /// \name Access the arrangement dimensions.
  //@{

  /*! Check whether the arrangement is empty. */
  bool is_empty () const
  {
    return (vertices_begin() == vertices_end() &&
            halfedges_begin() == halfedges_end());
  }

  /*!
   * Check whether the arrangement is valid. In particular, check the
   * validity of each vertex, halfedge and face, their incidence relations
   * and the geometric properties of the arrangement.
   */
  bool is_valid() const;
  
  /*! Get the number of arrangement vertices. */
  Size number_of_vertices () const
  {
    return (dcel.size_of_vertices());
  }

  /*! Get the number of isolated arrangement vertices. */
  Size number_of_isolated_vertices () const
  {
    return (n_iso_verts);
  }

  /*! Get the number of arrangement halfedges (the result is always even). */
  Size number_of_halfedges () const
  {
    return (dcel.size_of_halfedges());
  }

  /*! Get the number of arrangement edges. */
  Size number_of_edges () const
  {
    return (dcel.size_of_halfedges() / 2);
  }

  /*! Get the number of arrangement faces. */
  Size number_of_faces () const
  {
    return (dcel.size_of_faces());
  }
  //@}

  /// \name Traversal functions for the arrangement vertices.
  //@{

  /*! Get an iterator for the first vertex in the arrangement. */
  Vertex_iterator vertices_begin() 
  { 
    return (Vertex_iterator (dcel.vertices_begin())); 
  }

  /*! Get a past-the-end iterator for the arrangement vertices. */
  Vertex_iterator vertices_end()
  {
    return (Vertex_iterator(dcel.vertices_end())); 
  }

  /*! Get a const iterator for the first vertex in the arrangement. */
  Vertex_const_iterator vertices_begin() const
  { 
    return (Vertex_const_iterator (dcel.vertices_begin())); 
  }
  
  /*! Get a past-the-end const iterator for the arrangement vertices. */
  Vertex_const_iterator vertices_end() const
  {
    return (Vertex_const_iterator (dcel.vertices_end())); 
  }
  //@}

  /// \name Traversal functions for the arrangement halfedges.
  //@{

  /*! Get an iterator for the first halfedge in the arrangement. */
  Halfedge_iterator halfedges_begin() 
  { 
    return (Halfedge_iterator (dcel.halfedges_begin())); 
  }

  /*! Get a past-the-end iterator for the arrangement halfedges. */
  Halfedge_iterator halfedges_end()
  {
    return (Halfedge_iterator(dcel.halfedges_end())); 
  }

  /*! Get a const iterator for the first halfedge in the arrangement. */
  Halfedge_const_iterator halfedges_begin() const
  { 
    return (Halfedge_const_iterator (dcel.halfedges_begin())); 
  }
  
  /*! Get a past-the-end const iterator for the arrangement halfedges. */
  Halfedge_const_iterator halfedges_end() const
  {
    return (Halfedge_const_iterator (dcel.halfedges_end())); 
  }
  //@}

  /// \name Traversal functions for the arrangement edges.
  //@{

  /*! Get an iterator for the first edge in the arrangement. */
  Edge_iterator edges_begin() 
  { 
    return (Edge_iterator (dcel.edges_begin())); 
  }

  /*! Get a past-the-end iterator for the arrangement edges. */
  Edge_iterator edges_end()
  {
    return (Edge_iterator(dcel.edges_end())); 
  }

  /*! Get a const iterator for the first edge in the arrangement. */
  Edge_const_iterator edges_begin() const
  { 
    return (Edge_const_iterator (dcel.edges_begin())); 
  }
  
  /*! Get a past-the-end const iterator for the arrangement edges. */
  Edge_const_iterator edges_end() const

  {
    return (Edge_const_iterator (dcel.edges_end())); 
  }
  //@}

  /// \name Traversal functions for the arrangement faces.
  //@{

  /*! Get the unbounded face (non-const version). */
  Face_handle unbounded_face ()
  {
    return (Face_handle (un_face));
  }

  /*! Get the unbounded face (const version). */
  Face_const_handle unbounded_face () const
  {

    return (Face_const_handle (un_face));
  }

  /*! Get an iterator for the first face in the arrangement. */
  Face_iterator faces_begin() 
  { 
    return (Face_iterator (dcel.faces_begin())); 
  }

  /*! Get a past-the-end iterator for the arrangement faces. */
  Face_iterator faces_end()
  {
    return (Face_iterator(dcel.faces_end())); 
  }

  /*! Get a const iterator for the first face in the arrangement. */
  Face_const_iterator faces_begin() const
  { 
    return (Face_const_iterator (dcel.faces_begin())); 
  }
  
  /*! Get a past-the-end const iterator for the arrangement faces. */
  Face_const_iterator faces_end() const
  {
    return (Face_const_iterator (dcel.faces_end())); 
  }
  //@}

  /// \name Casting away constness for handle types.
  //@{
  Vertex_handle non_const_handle (Vertex_const_handle vh)
  {
    DVertex    *p_v = (DVertex*) &(*vh);
    return (Vertex_handle (p_v));
  }

  Halfedge_handle non_const_handle (Halfedge_const_handle hh)
  {
    DHalfedge  *p_he = (DHalfedge*) &(*hh);
    return (Halfedge_handle (p_he));
  }

  Face_handle non_const_handle (Face_const_handle fh)
  {
    DFace      *p_f = (DFace*) &(*fh);
    return (Face_handle (p_f));
  }
  //@}

  /// \name Vertex manipulation functions.
  //@{

  /*!
   * Replace the point associated with the given vertex.
   * \param v The vertex to modify.
   * \param p The point that should be associated with the edge.
   * \pre p is geometrically equivalent to the current point
   *      associated with v.
   * \return A handle for a the modified vertex (same as v).
   */
  Vertex_handle modify_vertex (Vertex_handle v,
                               const Point_2& p);

  /*!
   * Insert an isolated vertex in the interior of a given face.
   * \param p The isolated point.
   * \param f The face into which we insert the new isolated vertex.
   * \return A handle for the isolated vertex that has been created.
   */
  Vertex_handle insert_isolated_vertex (const Point_2& p,
                                        Face_handle f);

  /*!
   * Remove an isolated vertex from the interior of a given face.
   * \param v The vertex to remove.
   * \pre v is an isolated vertex (it has no incident halfedges).
   * \return A handle for the face containing v.
   */
  Face_handle remove_isolated_vertex (Vertex_handle v);
  
  ///@}

  /// \name Specilaized insertion functions for x-monotone curves.
  //@{

  /*!
   * Insert an x-monotone curve into the arrangement as a new hole (inner
   * component) inside the given face.
   * \param cv The given x-monotone curve.
   * \param f The face into which we insert the new hole.
   * \return A handle for one of the halfedges corresponding to the inserted
   *         curve, directed (lexicographically) from left to right.
   */
  Halfedge_handle insert_in_face_interior (const X_monotone_curve_2& cv, 
                                           Face_handle f);

  /*!
   * Insert an x-monotone curve into the arrangement, such that its left
   * endpoint corresponds to a given arrangement vertex.
   * \param cv The given x-monotone curve.
   * \param v The given vertex.
   * \pre The left endpoint of cv is incident to the vertex v.
   * \return A handle for one of the halfedges corresponding to the inserted
   *         curve, whose target is the new vertex.
   */
  Halfedge_handle insert_from_left_vertex (const X_monotone_curve_2& cv, 
                                           Vertex_handle v);

  /*!
   * Insert an x-monotone curve into the arrangement, such that its left
   * endpoints corresponds to a given arrangement vertex, given the exact
   * place for the curve in the circular list around this vertex.
   * \param cv The given x-monotone curve.
   * \param prev The reference halfedge. We should represent cv as a pair
   *             of edges, one of them should become prev's successor.
   * \pre The target vertex of prev is cv's left endpoint.
   * \return A handle for one of the halfedges corresponding to the inserted
   *         curve, whose target is the new vertex that was created.
   */
  Halfedge_handle insert_from_left_vertex (const X_monotone_curve_2& cv,
                                           Halfedge_handle prev);

  /*!
   * Insert an x-monotone curve into the arrangement, such that its right
   * endpoint corresponds to a given arrangement vertex.
   * \param cv The given x-monotone curve.
   * \param v The given vertex.
   * \pre The right endpoint of cv is incident to the vertex v.
   * \return A handle for one of the halfedges corresponding to the inserted
   *         curve, whose target is the new vertex.
   */
  Halfedge_handle insert_from_right_vertex (const X_monotone_curve_2& cv, 
                                            Vertex_handle v);

  /*! 
   * Insert an x-monotone curve into the arrangement, such that its right
   * endpoints corresponds to a given arrangement vertex, given the exact
   * place for the curve in the circular list around this vertex.

   * \param cv The given x-monotone curve.
   * \param prev The reference halfedge. We should represent cv as a pair
   *             of edges, one of them should become prev's successor.
   * \pre The target vertex of prev is cv's right endpoint.
   * \return A handle for one of the halfedges corresponding to the inserted
   *         curve, whose target is the new vertex that was created.
   */
  Halfedge_handle insert_from_right_vertex (const X_monotone_curve_2& cv,
                                            Halfedge_handle prev);
  
  /*! 
   * Insert an x-monotone curve into the arrangement, such that both its
   * endpoints correspond to given arrangement vertices.
   * \param cv The given x-monotone curve.
   * \param v1 The first vertex.
   * \param v2 The second vertex.
   * \pre v1 and v2 corresponds to cv's endpoints.
   * \return A handle for one of the halfedges corresponding to the inserted
   *         curve directed from v1 to v2.
   */
  Halfedge_handle insert_at_vertices (const X_monotone_curve_2& cv, 
                                      Vertex_handle v1, 
                                      Vertex_handle v2);

  /*! 
   * Insert an x-monotone curve into the arrangement, such that both its
   * endpoints correspond to given arrangement vertices, given the exact
   * place for the curve in one of the circular lists around a vertex.
   * \param cv The given x-monotone curve.
   * \param prev1 The reference halfedge for the first vertex.
   * \param v2 The second vertex.
   * \pre The target vertex of prev1 and v2 corresponds to cv's endpoints.
   * \return A handle for one of the halfedges corresponding to the inserted
   *         curve directed from prev1 to v2.
   */
  Halfedge_handle insert_at_vertices (const X_monotone_curve_2& cv, 
                                      Halfedge_handle h1, 
                                      Vertex_handle v2);

  /*!
   * Insert an x-monotone curve into the arrangement, such that both its
   * endpoints correspond to given arrangement vertices, given the exact
   * place for the curve in both circular lists around these two vertices.
   * \param cv the given curve.
   * \param prev1 The reference halfedge for the first vertex.
   * \param prev2 The reference halfedge for the second vertex.
   * \pre The target vertices of prev1 and prev2 are cv's endpoints.
   * \return A handle for one of the halfedges corresponding to the inserted
   *         curve directed from prev1's target to prev2's target.
   */
  Halfedge_handle insert_at_vertices (const X_monotone_curve_2 & cv,
                                      Halfedge_handle prev1, 
                                      Halfedge_handle prev2);

  //@}

  /// \name Halfedge manipulation functions.
  //@{

  /*!
   * Replace the x-monotone curve associated with the given edge.
   * \param e The edge to modify.
   * \param cv The curve that should be associated with the edge.
   * \pre cv is geometrically equivalent to the current curve
   *      associated with e.
   * \return A handle for a the modified halfedge (same as e).
   */
  Halfedge_handle modify_edge (Halfedge_handle e, 
                               const X_monotone_curve_2& cv);

  /*!
   * Split a given edge into two, and associate the given x-monotone
   * curves with the split edges.
   * \param e The edge to split (one of the pair of twin halfegdes).
   * \param cv1 The curve that should be associated with the first split edge.
   * \param cv2 The curve that should be associated with the second split edge.

   * \pre cv1's source and cv2's target equal the endpoints of the curve
   *      currently assoicated with e (respectively), and cv1's target equals
   *      cv2's target, and this is the split point (ot vice versa).
   * \return A handle for the halfedge whose source is the source of the the
   *         original halfedge e, and whose target is the split point.
   */
  Halfedge_handle split_edge (Halfedge_handle e,
                              const X_monotone_curve_2& cv1, 
                              const X_monotone_curve_2& cv2);

  /*!
   * Merge two edges to form a single edge, and associate the given x-monotone
   * curve with the merged edge.
   * \param e1 The first edge to merge (one of the pair of twin halfegdes).
   * \param e2 The second edge to merge (one of the pair of twin halfegdes).
   * \param cv The curve that should be associated with merged edge.
   * \return A handle for the merged halfedge.
   */
  Halfedge_handle merge_edge (Halfedge_handle e1, 
                              Halfedge_handle e2, 
                              const X_monotone_curve_2& cv);              

  /*!
   * Remove an edge from the arrangement.
   * \param e The edge to remove (one of the pair of twin halfegdes).
   * \param remove_source Should the source vertex of e be removed if it
   *                      becomes isolated (true by default).
   * \param remove_target Should the target vertex of e be removed if it
   *                      becomes isolated (true by default).

   * \return A handle for the remaining face.
   */
  Face_handle remove_edge (Halfedge_handle e,
                           bool remove_source = true,
			   bool remove_target = true);
  
  //@}

protected:

  /// \name Allocating and de-allocating points and curves.
  //@{

  /*! Allocate a new point. */
  Stored_point_2 *_new_point (const Point_2& p)
  {
    Stored_point_2   *p_sp = points_alloc.allocate (1);

    points_alloc.construct (p_sp, p);
    points.push_back (*p_sp);
    return (p_sp);
  }

  /*! De-allocate a point. */
  void _delete_point (Point_2& p)
  {
    Stored_point_2   *p_sp = static_cast<Stored_point_2*> (&p);

    points.erase (p_sp);
    points_alloc.destroy (p_sp);
    points_alloc.deallocate (p_sp, 1);
  }

  /*! Allocate a new curve. */
  Stored_curve_2 *_new_curve (const X_monotone_curve_2& cv)
  {
    Stored_curve_2   *p_scv = curves_alloc.allocate (1);

    curves_alloc.construct (p_scv, cv);
    curves.push_back (*p_scv);
    return (p_scv);
  }

  /*! De-allocate a curve. */
  void _delete_curve (X_monotone_curve_2& cv)
  {
    Stored_curve_2   *p_scv = static_cast<Stored_curve_2*> (&cv);

    curves.erase (p_scv);
    curves_alloc.destroy (p_scv);
    curves_alloc.deallocate (p_scv, 1);
  }
  //@}

  /// \name Converting handles to pointers (for the arrangement accessor).
  //@{

  /*! Convert a vertex handle to a pointer to a DCEL vertex. */
  /*! Convert a vertex handle to a pointer to a DCEL vertex. */
  inline DVertex* _vertex (Vertex_handle vh) const
  {
    return (&(*vh));
  }

  /*! Convert a constant vertex handle to a pointer to a DCEL vertex. */
  inline const DVertex* _vertex (Vertex_const_handle vh) const
  {
    return (&(*vh));
  }

  /*! Convert a halfedge handle to a pointer to a DCEL halfedge. */
  inline DHalfedge* _halfedge (Halfedge_handle hh) const
  {
    return (&(*hh));
  }

  /*! Convert a constant halfedge handle to a pointer to a DCEL halfedge. */
  inline const DHalfedge* _halfedge (Halfedge_const_handle hh) const
  {
    return (&(*hh));
  }

  /*! Convert a face handle to a pointer to a DCEL face. */
  inline DFace* _face (Face_handle fh) const
  {
    return (&(*fh));
  }

  /*! Convert a constant face handle to a pointer to a DCEL face. */
  inline const DFace* _face (Face_const_handle fh) const
  {
    return (&(*fh));
  }
  //@}

  /// \name Converting pointers to handles (for the arrangement accessor).
  //@{

  /*! Convert a pointer to a DCEL vertex to a vertex handle. */
  Vertex_handle _handle_for (DVertex *v)
  {
    return (Vertex_handle (v));
  }

  /*! Convert a pointer to a DCEL vertex to a constant vertex handle. */
  Vertex_const_handle _const_handle_for (const DVertex *v) const
  {
    return (Vertex_const_handle (v));
  }

  /*! Convert a pointer to a DCEL halfedge to a halfedge handle. */
  Halfedge_handle _handle_for (DHalfedge *he)
  {
    return (Halfedge_handle (he));
  }


  /*! Convert a pointer to a DCEL halfedge to a constant halfedge handle. */
  Halfedge_const_handle _const_handle_for (const DHalfedge *he) const
  {
    return (Halfedge_const_handle (he));
  }

  /*! Convert a pointer to a DCEL face to a face handle. */
  Face_handle _handle_for (DFace *f)
  {
    return (Face_handle (f));
  }

  /*! Convert a pointer to a DCEL face to a constant face handle. */
  Face_const_handle _const_handle_for (const DFace *f) const
  {
    return (Face_const_handle (f));
  }
  //@}

  /// \name Auxiliary (protected) functions.
  //@{
   
  /*!
   * Locate the place for the given curve around the given vertex.
   * \param v The given arrangement vertex.
   * \param cv The given x-monotone curve.
   * \pre v is one of cv's endpoints.
   * \return A pointer to a halfedge whose target is v, where cv should be
   *         inserted between this halfedge and the next halfedge around this
   *         vertex (in a clockwise order).
   *         A NULL return value indicates a precondition violation.
   */
  DHalfedge* _locate_around_vertex (DVertex *v,
                                    const X_monotone_curve_2& cv) const;

  /*!
   * Compute the distance (in halfedges) between two halfedges.
   * \param e1 The source halfedge.
   * \param e2 The destination halfedge.
   * \return In case e1 and e2 belong to the same connected component, the 
   *         function returns number of boundary halfedges between the two 
   *         halfedges. Otherwise, it returns (-1).
   */
  int _halfedge_distance (const DHalfedge *e1, const DHalfedge *e2) const;

  /*!
   * Determine whether a given query halfedge lies in the interior of a new
   * face we are about to create, by connecting it with another halfedge
   * using a given x-monotone curve.
   * \param prev1 The query halfedge.
   * \param prev2 The other halfedge we are about to connect with prev1.
   * \param cv The x-monotone curve we use to connect prev1 and prev2.
   * \pre prev1 and prev2 belong to the same connected component, and by
   *      connecting them using cv we form a new face.
   * \return (true) if prev1 lies in the interior of the face we are about
   *         to create, (false) otherwise - in which case prev2 must lie
   *         inside this new face.
   */
  bool _is_inside_new_face (const DHalfedge *prev1,
                            const DHalfedge *prev2,
                            const X_monotone_curve_2& cv) const;

  /*!
   * Determine whether a given point lies within the region bounded by
   * a boundary of a connected component.
   * \param p The query point.
   * \param v A vertex associated with the point p (or NULL if no such vertex
   *          exists or a vertex may exist but it is not known).
   * \param he A halfedge on the boundary of the connected component.
   * \return (true) if the point lies within region, (false) otherwise.
   */
  bool _point_is_in (const Point_2& p,
                     const DVertex* v,
                     const DHalfedge* he) const;

  /*!
   * Move a given hole from one face to another.
   * \param from_face The face currently containing the hole.
   * \param to_face The face into which we should move the hole.
   * \param hole A DCEL holes iterator pointing at the hole.
   */
  void _move_hole (DFace *from_face,
                   DFace *to_face,
                   DHoles_iter hole);

  /*!
   * Insert the given vertex as an isolated vertex inside the given face.
   * \param f The face that should contain the isolated vertex.
   * \param v The isolated vertex.
   */
  void _insert_isolated_vertex (DFace *f,
				DVertex *v);

  /*!
   * Move a given isolated vertex from one face to another.
   * \param from_face The face currently containing the isolated vertex.
   * \param to_face The face into which we should move the isolated vertex.
   * \param vit A DCEL isolated vertices iterator pointing at the vertex.
   */
  void _move_isolated_vertex (DFace *from_face,
			      DFace *to_face,
			      DIsolated_vertices_iter vit);

  /*!
   * Check whether the given halfedge lies on the outer boundary of the given
   * face.
   * \param f The given face.
   * \param e The given halfedge.
   * \return A pointer to a halfedge on the outer boundary of f in case e lies
   *         on this outer boundary, or NULL if it does not.
   */
  DHalfedge* _is_on_outer_boundary (DFace *f, DHalfedge *e) const;

  /*!
   * Check whether the given halfedge lies on the inner boundary of the given
   * face.
   * \param f The given face.
   * \param e The given halfedge.
   * \return A pointer to a halfedge on the inner boundary of f in case e lies
   *         on this outer boundary, or NULL if it does not.
   */
  DHalfedge* _is_on_inner_boundary (DFace *f, DHalfedge *e) const;

  /*!
   * Find the hole represented by a given halfedge from the holes container
   * of a given face and earse this hole once it is found.
   * \param f The given face.
   * \param e The given halfedge.
   * \return Whether the hole was found and erased or not.
   */
  bool _find_and_erase_hole (DFace *f, DHalfedge* e);

  /*!
   * Find the vertex in the isolated vertices container of a given face and
   * earse this vertex once it is found.
   * \param f The given face.
   * \param v The isolated vertex.
   * \return Whether the vertex was found and erased or not.
   */
  bool _find_and_erase_isolated_vertex (DFace *f, DVertex* v);

  /*!
   * Create a new vertex and associate it with the given point.
   * \param p The point.
   * \return A pointer to the newly created vertex.
   */
  DVertex* _create_vertex (const Point_2& p);
  
  /*!
   * Insert an x-monotone curve into the arrangement, such that both its
   * endpoints correspond to free arrangement vertices (newly created vertices
   * or existing isolated vertices), so a new hole is formed in the face
   * that contains the two vertices.
   * \param cv The given x-monotone curve.
   * \param f The face containing the two end vertices.
   * \param v1 The free vertex that corresponds to the left endpoint of cv.
   * \param v2 The free vertex that corresponds to the right endpoint of cv.
   * \param res The comparison result of the points associated with v1 and v2.
   * \return A pointer to one of the halfedges corresponding to the inserted
   *         curve, directed from v1 to v2.
   */
  DHalfedge* _insert_in_face_interior (const X_monotone_curve_2& cv,
                                       DFace *f,
                                       DVertex *v1, DVertex *v2,
                                       Comparison_result res);

  /*! 
   * Insert an x-monotone curve into the arrangement, such that one of its
   * endpoints corresponds to a given arrangement vertex, given the exact
   * place for the curve in the circular list around this vertex. The other
   * endpoint corrsponds to a free vertex (a newly created vertex or an
   * isolated vertex).
   * \param cv The given x-monotone curve.
   * \param prev The reference halfedge. We should represent cv as a pair
   *             of edges, one of them should become prev's successor.
   * \param v The free vertex that corresponds to the other endpoint.
   * \param res The comparison result of the points associated with prev's
   *            target and v.
   * \return A pointer to one of the halfedges corresponding to the inserted
   *         curve, whose target is the vertex v.
   */
  DHalfedge* _insert_from_vertex (const X_monotone_curve_2& cv,
                                  DHalfedge *prev,
                                  DVertex *v,
                                  Comparison_result res);

  /*!
   * Insert an x-monotone curve into the arrangement, where the end vertices
   * are given by the target points of two given halfedges.
   * The two halfedges should be given such that in case a new face is formed,
   * it will be the incident face of the halfedge directed from the first
   * vertex to the second vertex.
   * \param cv the given curve.
   * \param prev1 The reference halfedge for the first vertex.
   * \param prev2 The reference halfedge for the second vertex.
   * \param res The comparison result of the points associated with prev1's
   *            target vertex and prev2's target vertex.
   * \param new_face Output - whether a new face has been created.
   * \return A pointer to one of the halfedges corresponding to the inserted
   *         curve directed from prev1's target to prev2's target.
   *         In case a new face has been created, it is given as the incident
   *         face of this halfedge.
   */
  DHalfedge* _insert_at_vertices (const X_monotone_curve_2& cv,
                                  DHalfedge *prev1, 
                                  DHalfedge *prev2,
                                  Comparison_result res,
                                  bool& new_face);

  /*!
   * Relocate all holes and isolated vertices to their proper position,
   * immediately after a face has split due to the insertion of a new halfedge.
   * \param new_he The new halfedge that caused the split, such that the new
   *               face lies to its left and the old face to its right.
   */
  void _relocate_in_new_face (DHalfedge *new_he);

  /*!
   * Replace the point associated with the given vertex.
   * \param v The vertex to modify.
   * \param p The point that should be associated with the edge.
   */
  void _modify_vertex (DVertex *v, 
                       const Point_2& p);

  /*!
   * Replace the x-monotone curve associated with the given edge.
   * \param e The edge to modify.
   * \param cv The curve that should be associated with the edge.
   */
  void _modify_edge (DHalfedge *he, 
                     const X_monotone_curve_2& cv);

  /*!
   * Split a given edge into two at a given point, and associate the given
   * x-monotone curves with the split edges.
   * \param e The edge to split (one of the pair of twin halfegdes).
   * \param p The split point.
   * \param cv1 The curve that should be associated with the first split edge,
   *            whose source equals e's source and its target is p.
   * \param cv2 The curve that should be associated with the second split edge,
   *            whose source is p and its target equals e's target.
   * \return A pointer to the first split halfedge, whose source equals the
   *         source of e, and whose target is the split point.
   */
  DHalfedge* _split_edge (DHalfedge *e,
                          const Point_2& p,
                          const X_monotone_curve_2& cv1, 
                          const X_monotone_curve_2& cv2);

  /*!
   * Split a given edge into two at a given vertex, and associate the given
   * x-monotone curves with the split edges.
   * \param e The edge to split (one of the pair of twin halfegdes).
   * \param v The split vertex.
   * \param cv1 The curve that should be associated with the first split edge,
   *            whose source equals e's source and its target is v.
   * \param cv2 The curve that should be associated with the second split edge,
   *            whose source is v and its target equals e's target.
   * \return A pointer to the first split halfedge, whose source equals the
   *         source of e, and whose target is v.
   */
  DHalfedge* _split_edge (DHalfedge *e,
                          DVertex *v,
                          const X_monotone_curve_2& cv1,
                          const X_monotone_curve_2& cv2);

  /*!
   * Remove a pair of twin halfedges from the arrangement.
   * \param e One of the halfedges to be removed.
   * \param remove_source Should the source vertex of e be removed if it
   *                      becomes isolated.
   * \param remove_target Should the target vertex of e be removed if it
   *                      becomes isolated.
   * \pre In case the removal causes the creation of a new hole, e should 
   *      point at this hole.
   * \return A pointer to the remaining face.
   */
  DFace *_remove_edge (DHalfedge *e,
		       bool remove_source, bool remove_target);

  //@}

  /// \name Auxiliary (protected) functions for validity checking.
  //@{

  /*! Check the validity of a given vertex. */
  bool _is_valid (Vertex_const_handle v) const;

  /*! Check the validity of a given halfedge. */
  bool _is_valid (Halfedge_const_handle he) const;

  /*! Check the validity of a given face. */
  bool _is_valid (Face_const_handle f) const;

  /*! Check the validity of a CCB of a given face. */
  bool _is_valid (Ccb_halfedge_const_circulator start,
                  Face_const_handle f) const;

  /*!
   * Check that all vertices are unique (no two vertices with the same 
   * geometric point.
   */
  bool _are_vertices_unique() const;
  
  /*! Check that the curves around a given vertex are ordered clockwise. */
  bool _are_curves_ordered_cw_around_vertrex (Vertex_const_handle v) const;
  
  //@}

protected:

  /// \name Managing and notifying the arrangement observers.
  //@{

  /*!
   * Register a new observer (so it starts receiving notifications).
   * \param p_obs A pointer to the observer object.
   */
  void _register_observer (Observer *p_obs)
  {
    observers.push_back (p_obs);
  }

  /*!
   * Unregister a new observer (so it stops receiving notifications).
   * \param p_obs A pointer to the observer object.
   * \return Whether the observer was successfully unregistered.
   */
  bool _unregister_observer (Observer *p_obs)
  {
    Observers_iterator   iter;
    Observers_iterator   end = observers.end();

    for (iter = observers.begin(); iter != end; ++iter)
    {
      if ((*iter) == p_obs)
      {
        // Remove the p_ob pointer from the list of observers.
        observers.erase (iter);
        return (true);
      }
    }

    // If we reached here, the observer was not registered.
    return (false);
  }

private:

  /* Notify the observers on global arrangement operations: */

  void _notify_before_assign (const Self& arr)
  {
    Observers_iterator   iter;
    Observers_iterator   end = observers.end();

    for (iter = observers.begin(); iter != end; ++iter)
      (*iter)->before_assign (arr);
  }

  void _notify_after_assign ()
  {
    Observers_rev_iterator   iter;
    Observers_rev_iterator   end = observers.rend();

    for (iter = observers.rbegin(); iter != end; ++iter)
      (*iter)->after_assign();
  }

  void _notify_before_clear ()
  {
    Observers_iterator   iter;
    Observers_iterator   end = observers.end();

    for (iter = observers.begin(); iter != end; ++iter)
      (*iter)->before_clear();
  }

  void _notify_after_clear (Face_handle u)
  {
    Observers_rev_iterator   iter;
    Observers_rev_iterator   end = observers.rend();

    for (iter = observers.rbegin(); iter != end; ++iter)
      (*iter)->after_clear (u);
  }

  void _notify_before_global_change ()
  {
    Observers_iterator   iter;

    Observers_iterator   end = observers.end();

    for (iter = observers.begin(); iter != end; ++iter)
      (*iter)->before_global_change();

  }

  void _notify_after_global_change ()
  {
    Observers_rev_iterator   iter;
    Observers_rev_iterator   end = observers.rend();

    for (iter = observers.rbegin(); iter != end; ++iter)
      (*iter)->after_global_change();
  }

  /* Notify the observers on local changes in the arrangement: */
  
  void _notify_before_create_vertex (const Point_2& p)
  {
    Observers_iterator   iter;
    Observers_iterator   end = observers.end();

    for (iter = observers.begin(); iter != end; ++iter)
      (*iter)->before_create_vertex (p);
  }

  void _notify_after_create_vertex (Vertex_handle v)
  {
    Observers_rev_iterator   iter;
    Observers_rev_iterator   end = observers.rend();

    for (iter = observers.rbegin(); iter != end; ++iter)
      (*iter)->after_create_vertex (v);
  }

  void _notify_before_create_edge (const X_monotone_curve_2& c,
                                   Vertex_handle v1, Vertex_handle v2)
  {
    Observers_iterator   iter;
    Observers_iterator   end = observers.end();

    for (iter = observers.begin(); iter != end; ++iter)
      (*iter)->before_create_edge (c, v1, v2);
  }

  void _notify_after_create_edge (Halfedge_handle e)
  {
    Observers_rev_iterator   iter;
    Observers_rev_iterator   end = observers.rend();

    for (iter = observers.rbegin(); iter != end; ++iter)
      (*iter)->after_create_edge (e);
  }

  void _notify_before_modify_vertex (Vertex_handle v,
                                     const Point_2& p)
  {
    Observers_iterator   iter;
    Observers_iterator   end = observers.end();

    for (iter = observers.begin(); iter != end; ++iter)
      (*iter)->before_modify_vertex (v, p);
  }

  void _notify_after_modify_vertex (Vertex_handle v)
  {
    Observers_rev_iterator   iter;
    Observers_rev_iterator   end = observers.rend();

    for (iter = observers.rbegin(); iter != end; ++iter)
      (*iter)->after_modify_vertex (v);
  }

  void _notify_before_modify_edge (Halfedge_handle e,
                                   const X_monotone_curve_2& c)
  {
    Observers_iterator   iter;
    Observers_iterator   end = observers.end();

    for (iter = observers.begin(); iter != end; ++iter)
      (*iter)->before_modify_edge (e, c);
  }

  void _notify_after_modify_edge (Halfedge_handle e)
  {
    Observers_rev_iterator   iter;
    Observers_rev_iterator   end = observers.rend();

    for (iter = observers.rbegin(); iter != end; ++iter)
      (*iter)->after_modify_edge (e);
  }

  void _notify_before_split_edge (Halfedge_handle e,
                                  Vertex_handle v,
                                  const X_monotone_curve_2& c1,
                                  const X_monotone_curve_2& c2)
  {
    Observers_iterator   iter;
    Observers_iterator   end = observers.end();

    for (iter = observers.begin(); iter != end; ++iter)
      (*iter)->before_split_edge (e, v, c1, c2);
  }

  void _notify_after_split_edge (Halfedge_handle e1,
                                 Halfedge_handle e2)
  {
    Observers_rev_iterator   iter;
    Observers_rev_iterator   end = observers.rend();

    for (iter = observers.rbegin(); iter != end; ++iter)
      (*iter)->after_split_edge (e1, e2);
  }

  void _notify_before_split_face (Face_handle f,
                                  Halfedge_handle e)
  {
    Observers_iterator   iter;
    Observers_iterator   end = observers.end();

    for (iter = observers.begin(); iter != end; ++iter)
      (*iter)->before_split_face (f, e);
  }

  void _notify_after_split_face (Face_handle f,
                                 Face_handle new_f,
                                 bool is_hole)
  {
    Observers_rev_iterator   iter;
    Observers_rev_iterator   end = observers.rend();

    for (iter = observers.rbegin(); iter != end; ++iter)
      (*iter)->after_split_face (f, new_f, is_hole);
  }

  void _notify_before_split_hole (Face_handle f,
				  Ccb_halfedge_circulator h,
                                  Halfedge_handle e)
  {
    Observers_iterator   iter;
    Observers_iterator   end = observers.end();

    for (iter = observers.begin(); iter != end; ++iter)
      (*iter)->before_split_hole (f, h, e);
  }

  void _notify_after_split_hole (Face_handle f,
                                 Ccb_halfedge_circulator h1,
                                 Ccb_halfedge_circulator h2)
  {
    Observers_rev_iterator   iter;
    Observers_rev_iterator   end = observers.rend();

    for (iter = observers.rbegin(); iter != end; ++iter)
      (*iter)->after_split_hole (f, h1, h2);
  }

  void _notify_before_add_hole (Face_handle f,
                                Halfedge_handle e)
  {
    Observers_iterator   iter;
    Observers_iterator   end = observers.end();

    for (iter = observers.begin(); iter != end; ++iter)
      (*iter)->before_add_hole (f, e);
  }

  void _notify_after_add_hole (Ccb_halfedge_circulator h)
  {
    Observers_rev_iterator   iter;
    Observers_rev_iterator   end = observers.rend();

    for (iter = observers.rbegin(); iter != end; ++iter)
      (*iter)->after_add_hole (h);
  }

  void _notify_before_add_isolated_vertex (Face_handle f,
                                           Vertex_handle v)
  {
    Observers_iterator   iter;
    Observers_iterator   end = observers.end();

    for (iter = observers.begin(); iter != end; ++iter)
      (*iter)->before_add_isolated_vertex (f, v);
  }

  void _notify_after_add_isolated_vertex (Vertex_handle v)
  {
    Observers_rev_iterator   iter;
    Observers_rev_iterator   end = observers.rend();

    for (iter = observers.rbegin(); iter != end; ++iter)
      (*iter)->after_add_isolated_vertex (v);
  }

  void _notify_before_merge_edge (Halfedge_handle e1,
                                  Halfedge_handle e2,
                                  const X_monotone_curve_2& c)
  {
    Observers_iterator   iter;
    Observers_iterator   end = observers.end();

    for (iter = observers.begin(); iter != end; ++iter)
      (*iter)->before_merge_edge (e1, e2, c);
  }

  void _notify_after_merge_edge (Halfedge_handle e)
  {
    Observers_rev_iterator   iter;
    Observers_rev_iterator   end = observers.rend();

    for (iter = observers.rbegin(); iter != end; ++iter)
      (*iter)->after_merge_edge (e);
  }

  void _notify_before_merge_face (Face_handle f1,
                                  Face_handle f2,
                                  Halfedge_handle e)
  {
    Observers_iterator   iter;
    Observers_iterator   end = observers.end();

    for (iter = observers.begin(); iter != end; ++iter)
      (*iter)->before_merge_face (f1, f2, e);
  }

  void _notify_after_merge_face (Face_handle f)
  {
    Observers_rev_iterator   iter;
    Observers_rev_iterator   end = observers.rend();

    for (iter = observers.rbegin(); iter != end; ++iter)
      (*iter)->after_merge_face (f);
  }

  void _notify_before_merge_hole (Face_handle f,
				  Ccb_halfedge_circulator h1,
                                  Ccb_halfedge_circulator h2,
                                  Halfedge_handle e)
  {
    Observers_iterator   iter;
    Observers_iterator   end = observers.end();

    for (iter = observers.begin(); iter != end; ++iter)
      (*iter)->before_merge_hole (f, h1, h2, e);
  }

  void _notify_after_merge_hole (Face_handle f,
				 Ccb_halfedge_circulator h)
  {
    Observers_rev_iterator   iter;
    Observers_rev_iterator   end = observers.rend();

    for (iter = observers.rbegin(); iter != end; ++iter)
      (*iter)->after_merge_hole (f, h);
  }

  void _notify_before_move_hole (Face_handle from_f,
                                 Face_handle to_f,
                                 Ccb_halfedge_circulator h)
  {
    Observers_iterator   iter;
    Observers_iterator   end = observers.end();

    for (iter = observers.begin(); iter != end; ++iter)
      (*iter)->before_move_hole (from_f, to_f, h);
  }

  void _notify_after_move_hole (Ccb_halfedge_circulator h)
  {
    Observers_rev_iterator   iter;
    Observers_rev_iterator   end = observers.rend();

    for (iter = observers.rbegin(); iter != end; ++iter)
      (*iter)->after_move_hole (h);
  }

  void _notify_before_move_isolated_vertex (Face_handle from_f,
                                            Face_handle to_f,
                                            Vertex_handle v)
  {
    Observers_iterator   iter;
    Observers_iterator   end = observers.end();

    for (iter = observers.begin(); iter != end; ++iter)
      (*iter)->before_move_isolated_vertex (from_f, to_f, v);
  }


  void _notify_after_move_isolated_vertex (Vertex_handle v)
  {
    Observers_rev_iterator   iter;
    Observers_rev_iterator   end = observers.rend();

    for (iter = observers.rbegin(); iter != end; ++iter)
      (*iter)->after_move_isolated_vertex (v);
  }

  void _notify_before_remove_vertex (Vertex_handle v)
  {
    Observers_iterator   iter;
    Observers_iterator   end = observers.end();

    for (iter = observers.begin(); iter != end; ++iter)
      (*iter)->before_remove_vertex (v);
  }

  void _notify_after_remove_vertex ()
  {
    Observers_rev_iterator   iter;
    Observers_rev_iterator   end = observers.rend();

    for (iter = observers.rbegin(); iter != end; ++iter)
      (*iter)->after_remove_vertex ();
  }

  void _notify_before_remove_edge (Halfedge_handle e)
  {
    Observers_iterator   iter;
    Observers_iterator   end = observers.end();

    for (iter = observers.begin(); iter != end; ++iter)
      (*iter)->before_remove_edge (e);
  }

  void _notify_after_remove_edge ()
  {
    Observers_rev_iterator   iter;
    Observers_rev_iterator   end = observers.rend();

    for (iter = observers.rbegin(); iter != end; ++iter)
      (*iter)->after_remove_edge ();
  }

  void _notify_before_remove_hole (Face_handle f,
                                   Ccb_halfedge_circulator h)
  {
    Observers_iterator   iter;
    Observers_iterator   end = observers.end();

    for (iter = observers.begin(); iter != end; ++iter)
      (*iter)->before_remove_hole (f, h);
  }

  void _notify_after_remove_hole (Face_handle f)
  {
    Observers_rev_iterator   iter;
    Observers_rev_iterator   end = observers.rend();

    for (iter = observers.rbegin(); iter != end; ++iter)
      (*iter)->after_remove_hole (f);
  }
  //@}

};

//-----------------------------------------------------------------------------
// Declarations of the various global insertion and removal functions.
//-----------------------------------------------------------------------------

/*!
 * Insert a curve into the arrangement (incremental insertion).
 * The inserted curve may not necessarily be x-monotone and may intersect the
 * existing arrangement.
 * \param arr The arrangement.
 * \param cv The curve to be inserted.
 * \param pl A point-location object associated with the arrangement.
 */
template <class Traits, class Dcel, class PointLocation>
void insert_curve (Arrangement_2<Traits,Dcel>& arr,
                   const typename Traits::Curve_2& c,
                   const PointLocation& pl);

/*!
 * Insert a curve into the arrangement (incremental insertion).
 * The inserted curve may not necessarily be x-monotone and may intersect the
 * existing arrangement. The default "walk" point-location strategy is used
 * for the curve insertion.
 * \param arr The arrangement.
 * \param cv The curve to be inserted.
 */
template <class Traits, class Dcel>
void insert_curve (Arrangement_2<Traits,Dcel>& arr,
                   const typename Traits::Curve_2& c);

/*!
 * Insert a range of curves into the arrangement (aggregated insertion). 
 * The inserted curves may intersect one another and may also intersect the 
 * existing arrangement.
 * \param arr The arrangement.
 * \param begin An iterator for the first curve in the range.
 * \param end A past-the-end iterator for the curve range.
 * \pre The value type of the iterators must be Curve_2.
 */
template <class Traits, class Dcel, class InputIterator>
void insert_curves (Arrangement_2<Traits,Dcel>& arr,
                    InputIterator begin, InputIterator end);

/*!
 * Insert an x-monotone curve into the arrangement (incremental insertion).
 * The inserted x-monotone curve may intersect the existing arrangement.
 * \param arr The arrangement.
 * \param cv The x-monotone curve to be inserted.
 * \param pl A point-location object associated with the arrangement.
 */
template <class Traits, class Dcel, class PointLocation>
void insert_x_monotone_curve (Arrangement_2<Traits,Dcel>& arr,
                              const typename Traits::X_monotone_curve_2& c,
                              const PointLocation& pl);

/*!
 * Insert an x-monotone curve into the arrangement (incremental insertion)
 * when the location of the left endpoint of the curve is known and is
 * given as an isertion hint.
 * The inserted x-monotone curve may intersect the existing arrangement.
 * \param arr The arrangement.
 * \param cv The x-monotone curve to be inserted.
 * \param obj An object that represents the location of cv's left endpoint
 *            in the arrangement.
 */

template <class Traits, class Dcel>
void insert_x_monotone_curve (Arrangement_2<Traits,Dcel>& arr,
                              const typename Traits::X_monotone_curve_2& c,
                              const Object& obj);

/*!
 * Insert an x-monotone curve into the arrangement (incremental insertion).
 * The inserted x-monotone curve may intersect the existing arrangement.
 * The default "walk" point-location strategy is used for the curve insertion.
 * \param arr The arrangement.
 * \param cv The x-monotone curve to be inserted.
 */
template <class Traits, class Dcel>
void insert_x_monotone_curve (Arrangement_2<Traits,Dcel>& arr,
                              const typename Traits::X_monotone_curve_2& c);

/*!
 * Insert a range of x-monotone curves into the arrangement (aggregated
 * insertion). The inserted x-monotone curves may intersect one another and
 * may also intersect the existing arrangement.
 * \param arr The arrangement.
 * \param begin An iterator for the first curve in the range.
 * \param end A past-the-end iterator for the curve range.
 * \pre The value type of the iterators must be X_monotone_curve_2.
 */
template <class Traits, class Dcel, class InputIterator>
void insert_x_monotone_curves (Arrangement_2<Traits,Dcel>& arr,
                               InputIterator begin, InputIterator end);

/*!
 * Insert an x-monotone curve into the arrangement, such that the curve
 * interior does not intersect with any existing edge or vertex in the
 * arragement (incremental insertion).
 * \param arr The arrangement.
 * \param c The x-monotone curve to be inserted.
 * \param pl A point-location object associated with the arrangement.
 * \pre The interior of c does not intersect any existing edge or vertex.
 * \return A handle for one of the new halfedges corresponding to the inserted
 *         curve, directed (lexicographically) from left to right.
 */
template <class Traits, class Dcel, class PointLocation>
typename Arrangement_2<Traits,Dcel>::Halfedge_handle
insert_non_intersecting_curve (Arrangement_2<Traits,Dcel>& arr,
                               const typename Traits::X_monotone_curve_2& c,
                               const PointLocation& pl);

/*!
 * Insert an x-monotone curve into the arrangement, such that the curve
 * interior does not intersect with any existing edge or vertex in the
 * arragement (incremental insertion). The default "walk" point-location
 * strategy is used for the curve insertion.
 * \param arr The arrangement.
 * \param c The x-monotone curve to be inserted.
 * \pre The interior of c does not intersect any existing edge or vertex.

 * \return A handle for one of the new halfedges corresponding to the inserted
 *         curve, directed (lexicographically) from left to right.
 */
template <class Traits, class Dcel>
typename Arrangement_2<Traits,Dcel>::Halfedge_handle
insert_non_intersecting_curve (Arrangement_2<Traits,Dcel>& arr,
                               const typename Traits::X_monotone_curve_2& c);

/*!
 * Insert a range of pairwise interior-disjoint x-monotone curves into
 * the arrangement, such that the curve interiors do not intersect with
 * any existing edge or vertex in the arragement (aggregated insertion).
 * \param arr The arrangement.
 * \param begin An iterator for the first x-monotone curve in the range.
 * \param end A past-the-end iterator for the x-monotone curve range.
 * \pre The value type of the iterators must be X_monotone_curve_2.
 *      The curves in the range are pairwise interior-disjoint, and their
 *      interiors do not intersect any existing edge or vertex.
 */
template <class Traits, class Dcel, class InputIterator>
void insert_non_intersecting_curves (Arrangement_2<Traits,Dcel>& arr,
                                     InputIterator begin, InputIterator end);

/*!
 * Remove an edge from the arrangement. In case it is possible to merge
 * the edges incident to the end-vertices of the removed edge after its
 * deletion, the function performs these merges as well.
 * \param arr The arrangement.
 * \param e The edge to remove (one of the pair of twin halfegdes).
 * \return A handle for the remaining face.
 */
template <class Traits, class Dcel>
typename Arrangement_2<Traits,Dcel>::Face_handle
remove_edge (Arrangement_2<Traits,Dcel>& arr,
             typename Arrangement_2<Traits,Dcel>::Halfedge_handle e);

/*!
 * Insert a vertex that corresponds to a given point into the arrangement.
 * The inserted point may lie on any existing arrangement feature.
 * \param arr The arrangement.
 * \param p The point to be inserted.
 * \param pl A point-location object associated with the arrangement.
 * \return A handle to the vertex that corresponds to the given point.
 */
template <class Traits, class Dcel, class PointLocation>
typename Arrangement_2<Traits,Dcel>::Vertex_handle
insert_point (Arrangement_2<Traits,Dcel>& arr,
              const typename Traits::Point_2& p,
              const PointLocation& pl);

/*!
 * Insert a vertex that corresponds to a given point into the arrangement.
 * The inserted point may lie on any existing arrangement feature.
 * \param arr The arrangement.
 * \param p The point to be inserted.
 * \return A handle to the vertex that corresponds to the given point.
 */
template <class Traits, class Dcel>
typename Arrangement_2<Traits,Dcel>::Vertex_handle
insert_point (Arrangement_2<Traits,Dcel>& arr,
              const typename Traits::Point_2& p);

/*!
 * Remove a vertex from the arrangement.
 * \param arr The arrangement.
 * \param v The vertex to remove.
 * \return Whether the vertex has been removed or not.
 */
template <class Traits, class Dcel>
bool remove_vertex (Arrangement_2<Traits,Dcel>& arr,
                    typename Arrangement_2<Traits,Dcel>::Vertex_handle v);

CGAL_END_NAMESPACE

// The function definitions can be found under:
#include <CGAL/Arrangement_2/Arrangement_2_functions.h>
#include <CGAL/Arrangement_2/Arrangement_2_insert.h>

#endif

