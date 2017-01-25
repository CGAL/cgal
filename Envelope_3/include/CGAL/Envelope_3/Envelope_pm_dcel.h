// Copyright (c) 2005  Tel-Aviv University (Israel).
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
// Author(s)     : Michal Meyerovitch     <gorgymic@post.tau.ac.il>
//                 Baruch Zukerman        <baruchzu@post.tau.ac.il>
//                 Ron Wein               <wein@post.tau.ac.il>
//                 Efi Fogel              <efif@post.tau.ac.il>

#ifndef CGAL_ENVELOPE_PM_DCEL_H
#define CGAL_ENVELOPE_PM_DCEL_H

#include <CGAL/license/Envelope_3.h>


#include <CGAL/Arr_default_dcel.h>
#include <CGAL/Envelope_3/Envelope_base.h>

namespace CGAL
{
namespace Envelope_3
{

template <class Data>
class Dcel_info
{
public:
  typedef Dcel_info<Data>                         Self;
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
  Dcel_info() : m_is_set(false), m_decision(DAC_DECISION_NOT_SET)
  {}

  /*! \brief returns true iff data has been set already */
  bool get_is_set() const { return m_is_set; }

  /*! \brief resets the flag  */
  void set_is_set(bool flag) { m_is_set = flag; }

  bool is_decision_set()
  {
    return (m_decision != DAC_DECISION_NOT_SET);
  }
  
  Dac_decision get_decision() const
  {
    return m_decision;
  }                                                   

  void set_decision(Comparison_result comp)
  {
    m_decision = enum_cast<Dac_decision>(comp);
  }

  void set_decision(Dac_decision dec)
  {
    m_decision = dec;
  }

  /*! User-friendly interface: */
  size_t number_of_surfaces () const
  {
    return (m_data.size());
  }

  Data_const_iterator surfaces_begin () const
  {
    return (m_data.begin());
  }

  Data_const_iterator surfaces_end () const
  {
    return (m_data.end());
  }

   /*!
   * Get the first Xy-monotone surface associated with the face.
   * \pre number_of_surfaces() is not 0.
   */
  const Data& surface() const
  {
    CGAL_precondition (m_data.size() > 0);
    return (m_data.front());
  }

  /*!
   * Get the number of data objects associated with the face.
   */
  int number_of_data_objects() const
  {
    return static_cast<int>(m_data.size());
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
  bool is_equal_data(const InputIterator & begin, const InputIterator & end)
    const
  {
    if (!get_is_set())
      return false;

    // insert the input data objects into a set
    std::set<Data> input_data(begin, end);
    std::set<Data> my_data(begin_data(), end_data());
    if (input_data.size() != my_data.size())
      return false;
    return (my_data == input_data);
  }

  template <class InputIterator>
  bool has_equal_data(const InputIterator & begin, const InputIterator & end)
    const
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

public:
  template<class HandleType>
  void set_aux_source(unsigned int id, HandleType h)
  {
    CGAL_precondition(id < 2);
    m_aux_source[id] = make_object(h);
  }
  void set_aux_source(unsigned int id, const Object& o)
  {
    CGAL_precondition(id < 2);
    CGAL_precondition(!o.is_empty());
    m_aux_source[id] = o;
  }

  const Object& get_aux_source(unsigned int id)
  {
    CGAL_precondition(id < 2);
    CGAL_precondition (!m_aux_source[id].is_empty());
    return m_aux_source[id];
  }

  /*! \brief returns true iff the point has been set already */
  bool get_aux_is_set(unsigned int id) const
  {
    CGAL_precondition(id < 2);
	return (!m_aux_source[id].is_empty());
  }
};

/*! Extend the planar-map vertex */
template <class Point_2, class Data>
class Envelope_pm_vertex : public CGAL::Arr_vertex_base<Point_2>,
                           public Dcel_info<Data>
{
  typedef CGAL::Arr_vertex_base<Point_2>       Base_vertex;
  typedef Dcel_info<Data>                      Base_info;
  typedef Envelope_pm_vertex<Point_2, Data>    Self;

protected:
  // all flags are bits in this variable:
  unsigned short flags;

  // the flags indications:
  enum Bit_pos
  {
    // for an isolated vertex only
    IS_EQUAL   = 0,
    IS_EQUAL_AUX = 1,
    HAS_EQUAL = 3,
    HAS_EQUAL_AUX = 4,
    // indicate if the edge was added in the decomposition process
    // and is not part of the arrangement
    IS_FAKE = 6,
    // is this vertex an intersection vertex?
    // used in the partial vd, to eliminate vertical edges from
    // intersection points
    IS_INTERSECTION = 7
  };
public:
  /*! Constructor */
  Envelope_pm_vertex() : Dcel_info<Data>(), flags(0)
  {}

  /*void set_is_fake(bool b)
  {
    set_bit(IS_FAKE, b);
  }
  bool get_is_fake() const
  {
    return get_bit(IS_FAKE);
  }*/

 /* void set_is_intersection(bool b)
  {
    set_bit(IS_INTERSECTION, b);
  }*/
  /*bool get_is_intersection() const
  {
    return get_bit(IS_FAKE);
  }*/

  void set_is_equal_data_in_face(bool b)
  {
    set_bit(IS_EQUAL, b);
  }
  bool get_is_equal_data_in_face() const
  {
    return get_bit(IS_EQUAL);
  }

  void set_has_equal_data_in_face(bool b)
  {
    set_bit(HAS_EQUAL, b);
  }
  bool get_has_equal_data_in_face() const
  {
    return get_bit(HAS_EQUAL);
  }

  void set_is_equal_aux_data_in_face(unsigned int id, bool b)
  {
    CGAL_assertion(id < 2);
    set_bit(IS_EQUAL_AUX+id, b);
  }
  bool get_is_equal_aux_data_in_face(unsigned int id) const
  {
    CGAL_assertion(id < 2);
    return get_bit(IS_EQUAL_AUX+id);
  }

  void set_has_equal_aux_data_in_face(unsigned int id, bool b)
  {
    CGAL_assertion(id < 2);
    set_bit(HAS_EQUAL_AUX+id, b);
  }
  bool get_has_equal_aux_data_in_face(unsigned int id) const
  {
    CGAL_assertion(id < 2);
    return get_bit(HAS_EQUAL_AUX+id);
  }

  /*! Assign from another vertex.
   * \param v the other vertex.
   */
  virtual void assign(const Base_vertex & v)
  {
    Base_vertex::assign(v);

    const Self & ex_v = static_cast<const Self&>(v);
    this->Base_info::operator=(ex_v);
    flags = ex_v.flags;
  }
  
protected:
  void set_bit(unsigned int ind, bool b)
  {
    if (b)
      // set bit "ind" to 1:
      flags |= (1 << ind);
    else
      // set bit "ind" to 0:
      flags &= ~(1 << ind);
  }

  bool get_bit(unsigned int ind) const
  {
    // (1 << i) is bit i on, other bits off (start counting from 0)
    unsigned int mask = 1 << ind;
    return ((flags & mask) == mask);
  }
};

/*! Extend the planar-map halfedge */
template <class X_monotone_curve_2, class Data>
class
Envelope_pm_halfedge : public CGAL::Arr_halfedge_base<X_monotone_curve_2>,
                       public Dcel_info<Data>
{
  typedef CGAL::Arr_halfedge_base<X_monotone_curve_2>     Base_halfedge;
  typedef Dcel_info<Data>                                 Base_info;
  typedef Envelope_pm_halfedge<X_monotone_curve_2, Data>  Self;

protected:

  // all flags are bits in this variable:
  unsigned int flags;

  // flags indications
  enum Bit_pos
  {
    // indications for the Envelope algorithm
    // relation between halfedge and incident face
    IS_EQUAL_FACE   = 0,
    IS_EQUAL_AUX_FACE = 1,
    HAS_EQUAL_FACE = 3,
    HAS_EQUAL_AUX_FACE = 4,
    // relation between halfedge and target vertex
    IS_EQUAL_TARGET   = 6,
    IS_EQUAL_AUX_TARGET = 7,
    HAS_EQUAL_TARGET = 9,
    HAS_EQUAL_AUX_TARGET = 10,
    // relation between target vertex and incident face
    HAS_EQUAL_F_T = 12,
    HAS_EQUAL_AUX_F_T = 13,
    // indicate if the edge was added in the decomposition process
    // and is not part of the arrangement
    IS_FAKE = 15
  };
  
public:
  Envelope_pm_halfedge() : Dcel_info<Data>(), flags(0)
  {} 

 /* void set_is_fake(bool b)
  {
    set_bit(IS_FAKE, b);
  }
  bool get_is_fake() const
  {
    return get_bit(IS_FAKE);
  }*/

  void set_is_equal_data_in_face(bool b)
  {
    set_bit(IS_EQUAL_FACE, b);
  }
  bool get_is_equal_data_in_face() const
  {
    return get_bit(IS_EQUAL_FACE);
  }

  void set_has_equal_data_in_face(bool b)
  {
    set_bit(HAS_EQUAL_FACE, b);
  }
  bool get_has_equal_data_in_face() const
  {
    return get_bit(HAS_EQUAL_FACE);
  }

  void set_is_equal_aux_data_in_face(unsigned int id, bool b)
  {
    CGAL_assertion(id < 2);
    set_bit(IS_EQUAL_AUX_FACE+id, b);
  }
  bool get_is_equal_aux_data_in_face(unsigned int id) const
  {
    CGAL_assertion(id < 2);
    return get_bit(IS_EQUAL_AUX_FACE+id);
  }

  void set_has_equal_aux_data_in_face(unsigned int id, bool b)
  {
    CGAL_assertion(id < 2);
    set_bit(HAS_EQUAL_AUX_FACE+id, b);
  }
  bool get_has_equal_aux_data_in_face(unsigned int id) const
  {
    CGAL_assertion(id < 2);
    return get_bit(HAS_EQUAL_AUX_FACE+id);
  }

  void set_is_equal_data_in_target(bool b)
  {
    set_bit(IS_EQUAL_TARGET, b);
  }
  bool get_is_equal_data_in_target() const
  {
    return get_bit(IS_EQUAL_TARGET);
  }

  void set_has_equal_data_in_target(bool b)
  {
    set_bit(HAS_EQUAL_TARGET, b);
  }
  bool get_has_equal_data_in_target() const
  {
    return get_bit(HAS_EQUAL_TARGET);
  }

  void set_is_equal_aux_data_in_target(unsigned int id, bool b)
  {
    CGAL_assertion(id < 2);
    set_bit(IS_EQUAL_AUX_TARGET+id, b);
  }
  bool get_is_equal_aux_data_in_target(unsigned int id) const
  {
    CGAL_assertion(id < 2);
    return get_bit(IS_EQUAL_AUX_TARGET+id);
  }

  void set_has_equal_aux_data_in_target(unsigned int id, bool b)
  {
    CGAL_assertion(id < 2);
    set_bit(HAS_EQUAL_AUX_TARGET+id, b);
  }
  bool get_has_equal_aux_data_in_target(unsigned int id) const
  {
    CGAL_assertion(id < 2);
    return get_bit(HAS_EQUAL_AUX_TARGET+id);
  }
  // access to flags that contain relation between target and face
  void set_has_equal_data_in_target_and_face(bool b)
  {
    set_bit(HAS_EQUAL_F_T, b);
  }
  bool get_has_equal_data_in_target_and_face() const
  {
    return get_bit(HAS_EQUAL_F_T);
  }

  void set_has_equal_aux_data_in_target_and_face(unsigned int id, bool b)
  {
    CGAL_assertion(id < 2);
    set_bit(HAS_EQUAL_AUX_F_T+id, b);
  }
  bool get_has_equal_aux_data_in_target_and_face(unsigned int id) const
  {
    CGAL_assertion(id < 2);
    return get_bit(HAS_EQUAL_AUX_F_T+id);
  }

  /*! Assign from another halfedge.
   * \param h the other halfedge.
   */
  virtual void assign(const Base_halfedge & h)
  {
    Base_halfedge::assign(h);

    const Self & ex_h = static_cast<const Self&>(h);
    this->Base_info::operator=(ex_h);
    flags = ex_h.flags;
  }

protected:
  void set_bit(unsigned int ind, bool b)
  {
    if (b)
      // set bit "ind" to 1:
      flags |= (1 << ind);
    else
      // set bit "ind" to 0:
      flags &= ~(1 << ind);
    CGAL_assertion(get_bit(ind) == b);
  }

  bool get_bit(unsigned int ind) const
  {
    // (1 << i) is bit i on, other bits off (start counting from 0)
    unsigned int mask = 1 << ind;
    return ((flags & mask) == mask);
  }
};

/*! Extend the planar-map face */
template <class Data>
class Envelope_pm_face : public CGAL::Arr_face_base,
                         public Dcel_info<Data>
{
  typedef CGAL::Arr_face_base       Base_face;
  typedef Dcel_info<Data>           Base_info;
  typedef Envelope_pm_face<Data>    Self;

public:
  typedef std::list<Data>                         Data_container;
  typedef typename Data_container::iterator       Data_iterator;
  typedef typename Data_container::const_iterator Data_const_iterator;

  /*! Constructor */
  Envelope_pm_face() : Dcel_info<Data>()
  {}  

  /*! Assign from another face.
   * \param f the other face.
   */
  virtual void assign (const Base_face & f)
  {
    Base_face::assign(f);

    const Self & ex_f = static_cast<const Self&>(f);
    this->Base_info::operator=(ex_f);
  }
};

/*! A new dcel builder with full Envelope features */
template <class Traits, class Data>
class Envelope_pm_dcel : public
CGAL::Arr_dcel_base<Envelope_pm_vertex<typename Traits::Point_2, Data>,
                    Envelope_pm_halfedge<typename Traits::X_monotone_curve_2,
                                         Data>,
                    Envelope_pm_face<Data> >
{
public:

  typedef Data                                    Face_data;
  typedef typename Envelope_pm_face<Data>::Data_iterator
                                                  Face_data_iterator;
  typedef typename Envelope_pm_face<Data>::Data_const_iterator
                                                  Face_data_const_iterator;

  typedef Data                                    Edge_data;
  typedef Face_data_iterator                      Edge_data_iterator;
  typedef Face_data_const_iterator                Edge_data_const_iterator;

  typedef Data                                    Vertex_data;
  typedef Face_data_iterator                      Vertex_data_iterator;
  typedef Face_data_const_iterator                Vertex_data_const_iterator;

  typedef Dcel_info<Data>                         Dcel_elem_with_data;

  typedef Data                                    Dcel_data;
  typedef Face_data_iterator                      Dcel_data_iterator;
  typedef Face_data_const_iterator                Dcel_data_const_iterator;

  /*! \struct
   * An auxiliary structure for rebinding the DCEL with a new traits class.
   */
  template<typename T>
  struct rebind
  {
    typedef Envelope_pm_dcel<T, Face_data> other;
  };

  /*! Constructor */
  Envelope_pm_dcel() {}
};

} // namespace Envelope_3
} // namespace CGAL

#endif
