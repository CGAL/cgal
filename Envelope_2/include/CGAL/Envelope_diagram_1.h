// Copyright (c) 2011  Tel-Aviv University (Israel).
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
// Author(s)     : Ron Wein   <wein@post.tau.ac.il>

#ifndef CGAL_ENVELOPE_DIAGRAM_1_H
#define CGAL_ENVELOPE_DIAGRAM_1_H

#include <CGAL/license/Envelope_2.h>


#include <list>

#include <CGAL/basic.h>
#include <CGAL/memory.h>

namespace CGAL {

/*! \class
 * A minimization (or a maximization) diagram that represents the lower (or the
 * upper) envelope of a set of curves in the plane.
 */
template <class Traits_, class Allocator = CGAL_ALLOCATOR(int) >
class Envelope_diagram_1 {
public:
  typedef Traits_                                   Traits_2;
  typedef typename Traits_2::Point_2                Point_2;
  typedef typename Traits_2::X_monotone_curve_2     X_monotone_curve_2;
  typedef std::list<X_monotone_curve_2>             Curve_container;
  typedef typename Curve_container::iterator        Curve_iterator;
  typedef typename Curve_container::const_iterator  Curve_const_iterator;
  typedef std::size_t                               Size;

  /*!
   * Representation of a diagram vertex, which stores the point it represents
   * and a list of all curves that pass at that point.
   */
  class Edge;

  class Vertex {
  private:
    Point_2          _p;
    Curve_container  _cvs;
    Edge* _leftP;
    Edge* _rightP;

  public:
    /*! Constructor. */
    Vertex () :
      _leftP(NULL),
      _rightP(NULL)
    {}

    /*! Constructor with a point. */
    Vertex (const Point_2& p) :
      _p(p),
      _leftP(NULL),
      _rightP(NULL)
    {}

    /*! Get the point. */
    const Point_2& point () const
    {
      return (_p);
    }

    /*! Get the number of curves associated with the vertex. */
    Size number_of_curves () const
    {
      return (_cvs.size());
    }

    /*! Get the range of curves associated with the vertex. */
    Curve_const_iterator curves_begin () const
    {
      return (_cvs.begin());
    }

    Curve_const_iterator curves_end () const
    {
      return (_cvs.end());
    }

    /*! Get the left edge (const version). */
    const Edge* left () const
    {
      return (_leftP);
    }

    /*! Get the left edge (non-const version). */
    Edge* left ()
    {
      return (_leftP);
    }

    /*! Get the right edge (const version). */
    const Edge* right () const
    {
      return (_rightP);
    }

    /*! Get the right edge (non-const version). */
    Edge* right ()
    {
      return (_rightP);
    }

    /*! Set the point. */
    void set_point (const Point_2& p)
    {
      _p = p;
    }

    /*! Clear the list of curves. */
    void clear_curves ()
    {
      _cvs.clear();
    }

    /*! Associate a new curve with the vertex. */
    void add_curve (const X_monotone_curve_2& cv)
    {
      _cvs.push_back (cv);
    }

    /*! Associate a range of new curves with the vertex. */
    void add_curves (Curve_const_iterator begin, Curve_const_iterator end)
    {
      Curve_const_iterator  iter;

      for (iter = begin; iter != end; ++iter)
        _cvs.push_back (*iter);
    }

    /*! Set the left edge. */
    void set_left (Edge* e)
    {
      _leftP = e;
    }

    /*! Get the right edge. */
    void set_right (Edge* e)
    {
      _rightP = e;
    }
  };

  /*!
   * Representation of a diagram edge, which represents an interval and
   * stores all curves that form the envelope on it (usually there's just one
   * curve or zero curves, unless the input contains overlapping curves).
   */
  class Edge
  {
  private:
    Curve_container  _cvs;
    Vertex* _leftP;
    Vertex* _rightP;

  public:
    /*! Constructor. */
    Edge () :
      _leftP(NULL),
      _rightP(NULL)
    {}

    /*! Check if the edge represents an empty interval. */
    bool is_empty () const
    {
      return (_cvs.empty());
    }

    /*! Get the number of curves associated with the edge. */
    Size number_of_curves () const
    {
      return (_cvs.size());
    }

    /*!
     * Get a representative x-monotone curve.
     * \pre The edge does not represent an empty interval.
     */
    const X_monotone_curve_2& curve () const
    {
      CGAL_precondition (! _cvs.empty());
      return (*(_cvs.begin()));
    }

    /*! Get the range of curves associated with the edge. */
    Curve_const_iterator curves_begin () const
    {
      return (_cvs.begin());
    }

    Curve_const_iterator curves_end () const
    {
      return (_cvs.end());
    }

    /*! Get the left vertex (const version). */
    const Vertex* left () const
    {
      return (_leftP);
    }

    /*! Get the left vertex (const version). */
    Vertex* left ()
    {
      return (_leftP);
    }

    /*! Get the right vertex (const version). */
    const Vertex* right () const
    {
      return (_rightP);
    }

    /*! Get the right vertex (non-const version). */
    Vertex* right ()
    {
      return (_rightP);
    }

    /*! Clear the list of curves. */
    void clear_curves ()
    {
      _cvs.clear();
    }

    /*! Associate a new curve with the edge. */
    void add_curve (const X_monotone_curve_2& cv)
    {
      _cvs.push_back (cv);
    }

    /*! Associate a range of new curves with the edge. */
    void add_curves (Curve_const_iterator begin, Curve_const_iterator end)
    {
      Curve_const_iterator  iter;

      for (iter = begin; iter != end; ++iter)
        _cvs.push_back (*iter);
    }

    /*! Set the left vertex. */
    void set_left (Vertex* v)
    {
      _leftP = v;
    }

    /*! Get the right vertex. */
    void set_right (Vertex* v)
    {
      _rightP = v;
    }
  };

  typedef Vertex*                        Vertex_handle;
  typedef const Vertex*                  Vertex_const_handle;
  typedef Edge*                          Edge_handle;
  typedef const Edge*                    Edge_const_handle;

private:

  // Vertex allocator.
#ifdef CGAL_CXX11
    typedef std::allocator_traits<Allocator> Allocator_traits;
    typedef typename Allocator_traits::template rebind_alloc<Vertex> Vertex_allocator;
#else
  typedef typename Allocator::template rebind<Vertex>    Vertex_alloc_rebind;
  typedef typename Vertex_alloc_rebind::other            Vertex_allocator;
#endif

  // Halfedge allocator.
#ifdef CGAL_CXX11
    typedef typename Allocator_traits::template rebind_alloc<Edge> Edge_allocator;
#else
  typedef typename Allocator::template rebind<Edge>      Edge_alloc_rebind;
  typedef typename Edge_alloc_rebind::other              Edge_allocator;
#endif

  Edge* _leftmostP;                   // The leftmost edge of the diagram
                                      // (representing the range from -oo).
  Edge* _rightmostP;                  // The rightmost edge of the diagram
                                      // (representing the range to +oo).

  Vertex_allocator  vertex_alloc;     // An allocator for vertices.
  Edge_allocator    edge_alloc;       // An allocator for edges.

public:

  /*!
   * Constructor.
   */
  Envelope_diagram_1 ()
  {
    // An empty diagram contains one (empty) edge, representing (-oo, +oo):
    _leftmostP = _rightmostP = new_edge();
  }

  /*!
   * Destructor.
   */
  ~Envelope_diagram_1 ()
  {
    _free();
  }

  /*!
   * Get the leftmost edge of the diagram (const version).
   */
  Edge_const_handle leftmost () const
  {
    return (_leftmostP);
  }

  /*!
   * Get the leftmost edge of the diagram (non-const version).
   */
  Edge_handle leftmost ()
  {
    return (_leftmostP);
  }

  /*!
   * Get the rightmost edge of the diagram (const version).
   */
  Edge_const_handle rightmost () const
  {
    return (_rightmostP);
  }

  /*!
   * Get the rightmost edge of the diagram (non-const version).
   */
  Edge_handle rightmost ()
  {
    return (_rightmostP);
  }

  /*!
   * Clear the diagram.
   */
  void clear ()
  {
    _free();

    // An empty diagram contains one (empty) edge, representing (-oo, +oo):
    _leftmostP = _rightmostP = new_edge();
  }

  /*!
   * Set the leftmost edge.
   */
  void set_leftmost (Edge_handle e)
  {
    _leftmostP = e;
  }

  /*!
   * Set the rightmost edge.
   */
  void set_rightmost (Edge_handle e)
  {
    _rightmostP = e;
  }

  /*! Create a new vertex. */
  Vertex_handle new_vertex (const Point_2& p)
  {
    Vertex* v = vertex_alloc.allocate (1);
#ifdef CGAL_CXX11
    std::allocator_traits<Vertex_allocator>::construct(vertex_alloc, v, p);
#else
    vertex_alloc.construct (v, Vertex(p));
#endif
    return (v);
  }

  /*! Create a new edge. */
  Edge_handle new_edge ()
  {
    Edge* e = edge_alloc.allocate (1);
#ifdef CGAL_CXX11
    std::allocator_traits<Edge_allocator>::construct(edge_alloc, e);
#else
    edge_alloc.construct (e, Edge());
#endif
    return (e);
  }
   
  /*! Delete an existing vertex. */
  void delete_vertex (Vertex_handle v)
  {
#ifdef CGAL_CXX11
    std::allocator_traits<Vertex_allocator>::destroy(vertex_alloc, v);
#else
    vertex_alloc.destroy (v);
#endif
    vertex_alloc.deallocate (v, 1);
  }
  
  /*! Delete an existing edge. */
  void delete_edge (Edge_handle e)
  {
#ifdef CGAL_CXX11
    std::allocator_traits<Edge_allocator>::destroy(edge_alloc, e);
#else
    edge_alloc.destroy (e);
#endif
    edge_alloc.deallocate (e, 1);
  }
  
private:
  /*!
   * Free all diagram elements.
   */
  void _free ()
  {
    Vertex* v;
    Edge* e = _leftmostP;

    while (e != NULL) {
      // Get a pointer to the next vertex.
      v = e->right();

      // Free the edge and update it to be the next one after v.
      delete_edge (e);

      if (v != NULL) {
        e = v->right();

        // Free the current vertex.
        delete_vertex (v);
      }
      else
      {
        e = NULL;
      }
    }
     
    _leftmostP = NULL;
    _rightmostP = NULL;
  }
};

} //namespace CGAL

#endif
