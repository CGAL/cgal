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
// 
// Author(s)     : Idit Haran   <haranidi@post.tau.ac.il>
//                 Ron Wein     <wein@post.tau.ac.il>
//                 Efi Fogel    <efif@post.tau.ac.il>

#ifndef CGAL_ARR_LANDMARKS_NEAREST_NEIGHBOR_H
#define CGAL_ARR_LANDMARKS_NEAREST_NEIGHBOR_H

#include <CGAL/license/Arrangement_on_surface_2.h>


/*! \file
 * Definition of the Arr_landmarks_nearest_neighbor<Arrangement> template.
 */
#include <CGAL/basic.h>
#include <CGAL/Search_traits.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Arr_point_location_result.h>
#include <CGAL/Arrangement_2/Arr_traits_adaptor_2.h>

#include <boost/variant.hpp>
#include <boost/optional.hpp>

namespace CGAL {

/*! \class
 * A class that answers nearest neighbor queries.
 * It recieves a set of points, and builds a kd-tree for them.
 * Given a query point, it finds the closest point to the query.
 */
template <typename Arrangement_>
class Arr_landmarks_nearest_neighbor {
public:
  typedef Arrangement_                                  Arrangement_2;

  typedef typename Arrangement_2::Vertex_const_handle   Vertex_const_handle;
  typedef typename Arrangement_2::Halfedge_const_handle Halfedge_const_handle;
  typedef typename Arrangement_2::Face_const_handle     Face_const_handle;

  typedef Arr_point_location_result<Arrangement_2>      PL_result;
  typedef typename PL_result::Type                      PL_result_type;

  typedef typename Arrangement_2::Geometry_traits_2     Geometry_traits_2;
  typedef typename Geometry_traits_2::Approximate_number_type
    Approximate_number_type;
  typedef typename Geometry_traits_2::Point_2           Point_2;

  /*! \class NN_Point_2
   * Stores a point along with its approximate coordinates and its location
   * in the arrangement.
   */
  class  NN_Point_2 {
  public:
    Point_2 m_point;            // The point.
    PL_result_type m_object;    // The arrangement feature containing the point.
    Approximate_number_type m_vec[2]; // Approximate point x and y-coordinates.

  public:
    /*! Default constructor. */
    NN_Point_2() { m_vec[0] = m_vec[1] = 0; }

    /*! Constructor from a point. */
    NN_Point_2(const Point_2& p) :
      m_point(p)
    {
      // Obtain the coordinate approximations,
      Geometry_traits_2  m_traits;
      m_vec[0] = m_traits.approximate_2_object()(p, 0);
      m_vec[1] = m_traits.approximate_2_object()(p, 1);
    }

    /*! Constructor from a point and an its location in the arrangement. */
    NN_Point_2(const Point_2& p, const PL_result_type obj) :
      m_point(p),
      m_object(obj)
    {
      // Obtain the coordinate approximations,
      Geometry_traits_2  m_traits;
      m_vec[0] = m_traits.approximate_2_object()(p, 0);
      m_vec[1] = m_traits.approximate_2_object()(p, 1);
    }

    /* Get the point. */
    const Point_2& point() const { return (m_point); }

    /* Get the object representing the location in the arrangement. */
    const PL_result_type& object() const { return (m_object); }

    /*! Get an iterator for the approximate coordinates. */
    const Approximate_number_type* begin() const { return (m_vec); }

    /*! Get a past-the-end iterator for the approximate coordinates. */
    const Approximate_number_type* end() const { return (m_vec + 2); }

    /*! Equality operators. */
    bool operator== (const NN_Point_2& nnp) const
    { return (m_vec[0] == nnp.m_vec[0] && m_vec[1] == nnp.m_vec[1]); }

    bool operator!= (const NN_Point_2& nnp) const 
    { return (m_vec[0] != nnp.m_vec[0] || m_vec[1] != nnp.m_vec[1]); }
  };

  /*! \struct Construct_coord_iterator
   * An auxiliary structure that generates iterators (actually pointers) for
   * traversing the approximated point coordinates.
   */
  struct Construct_coord_iterator
  {
    typedef const Approximate_number_type*      result_type;
    
    /*! Get an iterator for the approximate coordinates. */
    const Approximate_number_type* operator()(const NN_Point_2& nnp) const
    { return (nnp.begin()); }

    /*! Get a past-the-end iterator for the approximate coordinates. */
    const Approximate_number_type* operator()(const NN_Point_2& nnp, int) const
    { return (nnp.end()); }
  };

protected:
  typedef CGAL::Search_traits<Approximate_number_type, NN_Point_2,
                              const Approximate_number_type*,
                              Construct_coord_iterator>     Search_traits;
  typedef CGAL::Orthogonal_k_neighbor_search<Search_traits> Neighbor_search;
  typedef typename Neighbor_search::iterator                Neighbor_iterator;
  typedef typename Neighbor_search::Tree                    Tree;

  // Data members:
  Tree* m_tree;        // The search tree.
  bool  m_is_empty;    // Is the search tree empty.

public: 
  bool is_empty() const { return m_is_empty; }
  
private:
  typedef Arr_landmarks_nearest_neighbor<Arrangement_2>     Self;

  /*! Copy constructor - not supported. */
  Arr_landmarks_nearest_neighbor(const Self&);

  /*! Assignment operator - not supported. */
  Self& operator=(const Self&);

public:
  /*! Default constructor. */
  Arr_landmarks_nearest_neighbor () :
    m_tree(NULL),
    m_is_empty(true)
  {}

  /*! Destructor. */
  ~Arr_landmarks_nearest_neighbor() { clear(); }

  /*!
   * Allocate the search tree and initialize it with landmark points.
   * \param begin An iterator for the first landmark point.
   * \param end A past-the-end iterator for the landmark points.
   * \pre The search tree is not initialized.
   */
  template <class InputIterator>
  void init(InputIterator begin, InputIterator end)
  {
    CGAL_precondition_msg(m_tree == NULL,
                          "The search tree is already initialized.");

    if (begin != end) {
      m_tree = new Tree(begin, end);
      m_is_empty = false;
    }
    else {
      m_tree = new Tree();
      m_is_empty = true;
    }
  }

  /*! Clear the search tree. */
  void clear() 
  {
    if (m_tree != NULL)
      delete m_tree;
    m_tree = NULL;
    m_is_empty = true;
  }

  /*!
   * Find the nearest landmark point to the query point.
   * \param q The query point.
   * \param obj Output: The location of the nearest landmark point in the
   *                    arrangement (a vertex, halfedge, or face handle).
   * \pre The search tree has been initialized and is not empty.
   * \return The nearest landmark point.
   */
  Point_2 find_nearest_neighbor(const Point_2& q, PL_result_type &obj) const
  {
    CGAL_precondition_msg(m_tree != NULL && ! m_is_empty,
                          "The search tree is not initialized.");

    // Create an NN_Point_2 object from the query point and use it to
    // query the search tree to find the nearest landmark point.
    NN_Point_2         nn_query(q);
    Neighbor_search    search(*m_tree, nn_query, 1);

    // For some reason search.begin()->first fails
    const NN_Point_2&  nearest_p = (*(search.begin())).first;
    obj = nearest_p.object();   
    return nearest_p.point();
  }
};

} //namespace CGAL

#endif
