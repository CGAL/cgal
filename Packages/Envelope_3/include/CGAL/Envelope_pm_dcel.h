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
// Author(s)     : Michal Meyerovitch     <gorgymic@post.tau.ac.il>

#ifndef CGAL_ENVELOPE_PM_DCEL_H
#define CGAL_ENVELOPE_PM_DCEL_H

#include <CGAL/basic.h>
#include <CGAL/Arr_default_dcel.h>
#include "CGAL/Envelope_base.h"
#include <CGAL/Unique_hash_map.h>
#include <iostream>
#include <set>

CGAL_BEGIN_NAMESPACE

template <class Data>
class Dcel_data
{
public:
  typedef Dcel_data<Data>                         Self;
  typedef std::list<Data>                         Data_container;
  typedef typename Data_container::iterator       Data_iterator;
  typedef typename Data_container::const_iterator Data_const_iterator;

protected:
  /*! data container */
  Data_container m_data;

  /*! Indicates that the data (surfaces) have been set already */
  bool m_is_set;

  // the decision that was made
  Dac_decision m_decision;
public:
  /*! Constructor */
  Dcel_data() : m_is_set(false), m_decision(NOT_SET)
  {}

  /*! \brief returns true iff data has been set already */
  bool get_is_set() const { return m_is_set; }

  /*! \brief resets the flag  */
  void set_is_set(bool flag) { m_is_set = flag; }

  bool is_decision_set()
  {
    return (m_decision != NOT_SET);
  }
  
  Dac_decision get_decision()
  {
    return m_decision;
  }                                                   

  void set_decision(Comparison_result comp)
  {
  	if (comp == SMALLER)
      m_decision = FIRST;
  	else if (comp == LARGER)
      m_decision = SECOND;
  	else
      m_decision = BOTH;
  }

  void set_decision(Dac_decision dec)
  {
    m_decision = dec;
  }
  
  /*!
   * Get the number of data objects associated with the face.
   */
  int number_of_data_objects() const
  {
    return (m_data.size());
  }

  /*!
   * check if the data is set to be empty
   */
  bool has_no_data() const
  {
    return (m_is_set && number_of_data_objects() == 0);
  }
  
  /*!
   * Get the first data object associated with the face.
   * \pre number_of_data_objects() is not 0.
   */
  const Data& get_data() const
  {
    CGAL_precondition (m_data.size() > 0);

    return (m_data.front());
  }

  /*!
   * Get the data iterators (const version).
   */
  Data_const_iterator begin_data() const
  {
    return (m_data.begin());
  }

  Data_const_iterator end_data() const
  {
    return (m_data.end());
  }

  /*!
   * Get the data iterators (non-const version).
   */
  Data_iterator begin_data()
  {
    return (m_data.begin());
  }

  Data_iterator end_data()
  {
    return (m_data.end());
  }


  /*!
   * Set a data object to the face.
   * \param data The data object to set.
   */
  void set_data (const Data & data)
  {
    clear_data();
    add_data(data);
  }

  /*!
   * Set a range of data objects to the face.
   * \param begin A begin iterator for the data range.
   * \param end A past-the-end iterator for the data range.
   */
  template <class InputIterator>
  void set_data(const InputIterator & begin, const InputIterator & end)
  {
    clear_data();
    add_data(begin, end);
  }

  /*!
   * set the data to be empty.
   */
  void set_no_data()
  {
    clear_data();
    m_is_set = true;
  }

   /*!
   * Add a data object to the face.
   * \param data The additional data object.
   */
  void add_data (const Data & data)
  {
    m_data.push_back(data);
    m_is_set = true;
  }

  /*!
   * Add a range of data objects to the face.
   * \param begin A begin iterator for the data range.
   * \param end A past-the-end iterator for the data range.
   */
  template <class InputIterator>
  void add_data (const InputIterator & begin, const InputIterator & end)
  {
    InputIterator    it;
    for (it = begin; it != end; it++)
      m_data.push_back(*it);
    m_is_set = true;
  }

  /*!
   * Clear the data objects.
   */
  void clear_data ()
  {
    m_data.clear();
    m_is_set = false;
  }

  /*!
   * Check if the set of data objects in the input range is equal to our
   * set of data objects
   */
  template <class InputIterator>
  bool is_equal_data(const InputIterator & begin, const InputIterator & end) const
  {
    if (!get_is_set())
      return false;
//    // insert all input data into a map, then try to find there all our data
//    Unique_hash_map<Data, bool> input_data;
//    InputIterator input = begin;
//    int input_size = 0;
//    for(; input != end; ++input, ++input_size)
//      input_data[*input] = true;
//    if (input_size != number_of_data_objects())
//      return false;
//    Data_const_iterator my_data = begin_data();
//    for(; my_data != end_data(); ++my_data)
//      if (!input_data.is_defined(*my_data))
//        return false;
//
//    return true;

    // insert the input data objects into a set
    std::set<Data> input_data(begin, end);
    std::set<Data> my_data(begin_data(), end_data());
    if (input_data.size() != my_data.size())
      return false;
    return (my_data == input_data);
  }

  template <class InputIterator>
  bool has_equal_data(const InputIterator & begin, const InputIterator & end) const
  {
    if (!get_is_set())
      return false;
    // insert the input data objects into a set
    std::set<Data> input_data(begin, end);
    std::set<Data> my_data(begin_data(), end_data());
    std::list<Data> intersection;
    std::set_intersection(my_data.begin(), my_data.end(),
                          input_data.begin(), input_data.end(),
                          std::back_inserter(intersection));
    return (intersection.size() > 0);
  }

protected:

  /*! Place holder for the source of the overlay data */
  Object m_aux_source[2];

//  /*! Indicates that the overlay data (surfaces) have been set already */
//  bool m_aux_is_set[2];

public:
  template<class HandleType>
  void set_aux_source(unsigned int id, HandleType h)
  {
    CGAL_precondition(id < 2);
    m_aux_source[id] = make_object(h);
  }
  void set_aux_source(unsigned int id, Object o)
  {
    CGAL_precondition(id < 2);
    CGAL_precondition(!o.is_empty());
    m_aux_source[id] = o;
  }

  Object get_aux_source(unsigned int id)
  {
    CGAL_precondition(id < 2);
    CGAL_precondition (!m_aux_source[id].is_empty());
    return m_aux_source[id];
  }

//  /*!
//   * Get the number of aux data objects associated with the face.
//   */
//  int number_of_aux_data_objects(unsigned int id) const
//  {
//    CGAL_precondition(id < 2);
//    if (!m_aux_is_set[id]) return 0;
//    return ((Self*)(m_aux_source[id]))->number_of_data_objects();
//  }
//
//  /*! \bried obtains idth map data */
//  const Data & get_aux_data(unsigned int id) const
//  {
//    CGAL_precondition(id < 2);
//    CGAL_precondition (m_aux_is_set[id]);
//    return ((Self*)(m_aux_source[id]))->get_data();
//  }
//
//  /*! \bried Get the auxiliary data iterators (const version). */
//  Data_const_iterator begin_aux_data(unsigned int id) const
//  {
//    CGAL_precondition(id < 2);
//    return ((Self*)(m_aux_source[id]))->begin_data();
//  }
//
//  Data_const_iterator end_aux_data(unsigned int id) const
//  {
//    CGAL_precondition(id < 2);
//    return ((Self*)(m_aux_source[id]))->end_data();
//  }
//
//  /*! \bried Get the auxiliary data iterators (non-const version). */
//  Data_iterator begin_aux_data(unsigned int id)
//  {
//    CGAL_precondition(id < 2);
//    return ((Self*)(m_aux_source[id]))->begin_data();
//  }
//
//  Data_iterator end_aux_data(unsigned int id)
//  {
//    CGAL_precondition(id < 2);
//    return ((Self*)(m_aux_source[id]))->end_data();
//  }

  /*! \brief returns true iff the point has been set already */
  bool get_aux_is_set(unsigned int id) const
  {
    CGAL_precondition(id < 2);
	return (!m_aux_source[id].is_empty());
//    return m_aux_is_set[id];
  }

//  /*! check if aux data is set to be empty */
//  bool aux_has_no_data(unsigned int id) const
//  {
//    CGAL_precondition(id < 2);
//    return (m_aux_is_set[id] && ((Self*)(m_aux_source[id]))->has_no_data());
//  }
//
//  /*!
//   * Check if the set of data objects in the input range is equal to our
//   * aux set of data objects
//   */
//  template <class InputIterator>
//  bool is_equal_aux_data(unsigned int id, const InputIterator & begin,
//                         const InputIterator & end) const
//  {
//    CGAL_precondition(id < 2);
//    if (!get_aux_is_set(id))
//      return false;
//    // insert the input data objects into a set
//    std::set<Data> input_data(begin, end);
//    std::set<Data> my_data(begin_aux_data(id), end_aux_data(id));
//    if (input_data.size() != my_data.size())
//      return false;
//    return (my_data == input_data);
//  }
//
//  /*!
//   * Check if the set of data objects has an equal data element as
//   * aux set of data objects
//   */
//  template <class InputIterator>
//  bool has_equal_aux_data(unsigned int id, const InputIterator & begin,
//                         const InputIterator & end) const
//  {
//    CGAL_precondition(id < 2);
//    if (!get_aux_is_set(id))
//      return false;
//    // insert the input data objects into a set
//    std::set<Data> input_data(begin, end);
//    std::set<Data> my_data(begin_aux_data(id), end_aux_data(id));
//    std::list<Data> intersection;
//    std::set_intersection(my_data.begin(), my_data.end(),
//                          input_data.begin(), input_data.end(),
//                          std::back_inserter(intersection));
//    return (intersection.size() > 0);
//  }
//

};



/*! Extend the planar-map vertex */
template <class Point_2, class Data>
class Envelope_pm_vertex : public CGAL::Arr_vertex_base<Point_2>,
                           public Dcel_data<Data>
{
private:
//  void* m_copy_from_halfedge;
  // indicate if the edge was added in the decomposition process
  // and is not part of the arrangement
  bool         m_is_fake;

  // is this vertex an intersection vertex?
  // used in the partial vd, to eliminate vertical edges from
  // intersection points
  bool         m_is_intersection;
  
  // indications for the Envelope algorithm (for an isolated vertex only)
  bool         m_is_equal_data_in_face; 
  bool         m_has_equal_data_in_face;
  bool         m_is_equal_aux_data_in_face[2]; 
  bool         m_has_equal_aux_data_in_face[2];

public:
  /*! Constructor */
  Envelope_pm_vertex() : Dcel_data<Data>()
      , m_is_fake(false)
      , m_is_intersection(false)
      , m_is_equal_data_in_face(false)
  	  , m_has_equal_data_in_face(false)
//    , m_copy_from_halfedge(NULL)
  {
    m_is_equal_aux_data_in_face[0] = m_is_equal_aux_data_in_face[1] = false;
    m_has_equal_aux_data_in_face[0] = m_has_equal_aux_data_in_face[1] = false;
  }

  void set_is_fake(bool b)
  {
    m_is_fake = b;
  }
  bool get_is_fake() const
  {
    return m_is_fake;
  }

  void set_is_intersection(bool b)
  {
    m_is_intersection = b;
  }
  bool get_is_intersection() const
  {
    return m_is_intersection;
  }

  void set_is_equal_data_in_face(bool b)
  {
    m_is_equal_data_in_face = b;
  }
  bool get_is_equal_data_in_face() const
  {
    return m_is_equal_data_in_face;
  }

  void set_has_equal_data_in_face(bool b)
  {
    m_has_equal_data_in_face = b;
  }
  bool get_has_equal_data_in_face() const
  {
    return m_has_equal_data_in_face;
  }

  void set_is_equal_aux_data_in_face(unsigned int id, bool b)
  {
    CGAL_assertion(id < 2);
    m_is_equal_aux_data_in_face[id] = b;
  }
  bool get_is_equal_aux_data_in_face(unsigned int id) const
  {
    CGAL_assertion(id < 2);
    return m_is_equal_aux_data_in_face[id];
  }

  void set_has_equal_aux_data_in_face(unsigned int id, bool b)
  {
    CGAL_assertion(id < 2);
    m_has_equal_aux_data_in_face[id] = b;
  }
  bool get_has_equal_aux_data_in_face(unsigned int id) const
  {
    CGAL_assertion(id < 2);
    return m_has_equal_aux_data_in_face[id];
  }

};

/*! Extend the planar-map halfedge */
template <class X_monotone_curve_2, class Data>
class Envelope_pm_halfedge : public CGAL::Arr_halfedge_base<X_monotone_curve_2>,
                             public Dcel_data<Data>
{
private:

  // indicate if the edge was added in the decomposition process
  // and is not part of the arrangement
  bool         m_is_fake;

  // indications for the Envelope algorithm
  bool         m_is_equal_data_in_face; 
  bool         m_has_equal_data_in_face;
  bool         m_is_equal_aux_data_in_face[2]; 
  bool         m_has_equal_aux_data_in_face[2];

  bool         m_is_equal_data_in_target;
  bool         m_has_equal_data_in_target;
  bool         m_is_equal_aux_data_in_target[2];
  bool         m_has_equal_aux_data_in_target[2];

public:
  Envelope_pm_halfedge() : Dcel_data<Data>()
                           , m_is_fake(false)
                           , m_is_equal_data_in_face(false)
							             , m_has_equal_data_in_face(false)
                           , m_is_equal_data_in_target(false)
							             , m_has_equal_data_in_target(false)
  {
    m_is_equal_aux_data_in_face[0] = m_is_equal_aux_data_in_face[1] = false;
    m_has_equal_aux_data_in_face[0] = m_has_equal_aux_data_in_face[1] = false;
    m_is_equal_aux_data_in_target[0] = m_is_equal_aux_data_in_target[1] = false;
    m_has_equal_aux_data_in_target[0] = m_has_equal_aux_data_in_target[1] = false;
  } 

  void set_is_fake(bool b)
  {
    m_is_fake = b;
  }
  bool get_is_fake() const
  {
    return m_is_fake;
  }

  void set_is_equal_data_in_face(bool b)
  {
    m_is_equal_data_in_face = b;
  }
  bool get_is_equal_data_in_face() const
  {
    return m_is_equal_data_in_face;
  }

  void set_has_equal_data_in_face(bool b)
  {
    m_has_equal_data_in_face = b;
  }
  bool get_has_equal_data_in_face() const
  {
    return m_has_equal_data_in_face;
  }

  void set_is_equal_aux_data_in_face(unsigned int id, bool b)
  {
    CGAL_assertion(id < 2);
    m_is_equal_aux_data_in_face[id] = b;
  }
  bool get_is_equal_aux_data_in_face(unsigned int id) const
  {
    CGAL_assertion(id < 2);
    return m_is_equal_aux_data_in_face[id];
  }

  void set_has_equal_aux_data_in_face(unsigned int id, bool b)
  {
    CGAL_assertion(id < 2);
    m_has_equal_aux_data_in_face[id] = b;
  }
  bool get_has_equal_aux_data_in_face(unsigned int id) const
  {
    CGAL_assertion(id < 2);
    return m_has_equal_aux_data_in_face[id];
  }

  void set_is_equal_data_in_target(bool b)
  {
    m_is_equal_data_in_target = b;
  }
  bool get_is_equal_data_in_target() const
  {
    return m_is_equal_data_in_target;
  }

  void set_has_equal_data_in_target(bool b)
  {
    m_has_equal_data_in_target = b;
  }
  bool get_has_equal_data_in_target() const
  {
    return m_has_equal_data_in_target;
  }

  void set_is_equal_aux_data_in_target(unsigned int id, bool b)
  {
    CGAL_assertion(id < 2);
    m_is_equal_aux_data_in_target[id] = b;
  }
  bool get_is_equal_aux_data_in_target(unsigned int id) const
  {
    CGAL_assertion(id < 2);
    return m_is_equal_aux_data_in_target[id];
  }

  void set_has_equal_aux_data_in_target(unsigned int id, bool b)
  {
    CGAL_assertion(id < 2);
    m_has_equal_aux_data_in_target[id] = b;
  }
  bool get_has_equal_aux_data_in_target(unsigned int id) const
  {
    CGAL_assertion(id < 2);
    return m_has_equal_aux_data_in_target[id];
  }

};

/*! Extend the planar-map face */
template <class Data>
class Envelope_pm_face : public CGAL::Arr_face_base,
                         public Dcel_data<Data>
{
public:
  typedef std::list<Data>                         Data_container;
  typedef typename Data_container::iterator       Data_iterator;
  typedef typename Data_container::const_iterator Data_const_iterator;

  /*! Constructor */
  Envelope_pm_face() : Dcel_data<Data>()
  {}  
};

/*! A new dcel builder with full Envelope features */
template <class Traits, class Data>
class Envelope_pm_dcel :
  public CGAL::Arr_dcel_base<Envelope_pm_vertex<typename Traits::Point_2, Data>,
                       Envelope_pm_halfedge<typename Traits::X_monotone_curve_2, Data>,
                       Envelope_pm_face<Data> >
{
public:
  typedef Data                                                    Face_data;
  typedef typename Envelope_pm_face<Data>::Data_iterator		      Face_data_iterator;
  typedef typename Envelope_pm_face<Data>::Data_const_iterator    Face_data_const_iterator;
  typedef Data                                                    Edge_data;
  typedef Face_data_iterator                                      Edge_data_iterator;
  typedef Face_data_const_iterator                                Edge_data_const_iterator;
  typedef Data                                                    Vertex_data;
  typedef Face_data_iterator                                      Vertex_data_iterator;
  typedef Face_data_const_iterator                                Vertex_data_const_iterator;

  typedef Dcel_data<Data>                                         Dcel_elem_with_data;

  typedef Data                                                    Dcel_data;
  typedef Face_data_iterator                                      Dcel_data_iterator;
  typedef Face_data_const_iterator                                Dcel_data_const_iterator;

  /*! Constructor */
  Envelope_pm_dcel() {}
};

CGAL_END_NAMESPACE

class Curve_data
{
private:
  /*! The id of the source planar map in the overlay */
  unsigned int m_pm_id;

  /*! The halfedge the curve is mapped to in the source planar map */
  void * m_halfedge;

  /*! Is the input curve oriented the same as the input halfedge */
  bool m_is_same_direction_in;

  // a temporary id
  int m_id;

public:
  /*! Constructor */
  Curve_data(unsigned int pm_id = (unsigned int) -1,
             void * halfedge = NULL,
             bool is_same_direction_in = true,
             int id = 0) :
    m_pm_id(pm_id),
    m_halfedge(halfedge),
    m_is_same_direction_in(is_same_direction_in),
    m_id(id)
  {}

  /*! obtains the halfedge the curve is mapped to in the source planar map */
  void * get_halfedge() const { return m_halfedge; }

  /*! obtains the flag that indicates whether the input curve is in the
   * same direction as the returned halfedge
   */
  bool get_is_same_direction_in() const { return m_is_same_direction_in; }

  /*! obtains the id of the source planar-map */
  unsigned int get_pm_id() const { return m_pm_id; }

  // return the temporary id
  int get_id() const { return m_id; }

  /*! Equality operator. */
  bool operator== (const Curve_data& cd) const
  {
    return true;
  }

};


#endif  // CGAL_ENVELOPE_PM_DCEL_H
