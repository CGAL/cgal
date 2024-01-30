// Copyright (c) 2005  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Michal Meyerovitch     <gorgymic@post.tau.ac.il>
//             Baruch Zukerman        <baruchzu@post.tau.ac.il>
//             Ron Wein               <wein@post.tau.ac.il>
//             Efi Fogel              <efif@post.tau.ac.il>

#ifndef CGAL_ENVELOPE_PM_DCEL_H
#define CGAL_ENVELOPE_PM_DCEL_H

#include <CGAL/license/Envelope_3.h>


#include <CGAL/Arr_default_dcel.h>
#include <CGAL/Envelope_3/Envelope_base.h>

namespace CGAL {

namespace Envelope_3 {

template <typename Data_>
class Dcel_info {
public:
  using Data = Data_;
  using Self = Dcel_info<Data>;
  using Data_container = std::list<Data>;
  using Data_iterator =  typename Data_container::iterator;
  using Data_const_iterator = typename Data_container::const_iterator;

protected:
  /*! data container */
  Data_container m_env_data;

  /*! Indicates that the data (surfaces) have been set already */
  bool m_is_env_set;

  // the decision that was made
  Dac_decision m_decision;
public:
  /*! Constructor */
  Dcel_info() : m_is_env_set(false), m_decision(DAC_DECISION_NOT_SET) {}

  /*! \brief returns true iff data has been set already */
  bool is_env_set() const { return m_is_env_set; }

  /*! \brief resets the flag  */
  void set_is_env_set(bool flag) { m_is_env_set = flag; }

  bool is_decision_set() { return (m_decision != DAC_DECISION_NOT_SET); }

  Dac_decision decision() const { return m_decision; }

  void set_decision(Comparison_result comp)
  { m_decision = enum_cast<Dac_decision>(comp); }

  void set_decision(Dac_decision dec) { m_decision = dec; }

  /*! User-friendly interface: */
  size_t number_of_surfaces() const { return m_env_data.size(); }

  Data_const_iterator surfaces_begin() const { return m_env_data.begin(); }

  Data_const_iterator surfaces_end() const { return m_env_data.end(); }

   /*! Obtain the first Xy-monotone surface associated with the face.
    * \pre number_of_surfaces() is not 0.
    */
  const Data& surface() const {
    CGAL_precondition(m_env_data.size() > 0);
    return m_env_data.front();
  }

  /*! Obtain the number of data objects associated with the cell.
   */
  int env_data_size() const
  { return static_cast<int>(m_env_data.size()); }

  /*! Check whether the data is set to be empty
   */
  bool has_no_env_data() const
  { return (m_is_env_set && (env_data_size() == 0)); }

  /*! Obtain the first data object associated with the cell.
   * \pre m_env_data.size() is not 0.
   */
  const Data& env_data_front() const {
    CGAL_precondition(m_env_data.size() > 0);
    return m_env_data.front();
  }

  /*! Obtain the data iterators (const version).
   */
  Data_const_iterator begin_env_data() const { return m_env_data.begin(); }

  Data_const_iterator end_env_data() const { return m_env_data.end(); }

  /*! Obtain the data iterators (non-const version).
   */
  Data_iterator begin_env_data() { return m_env_data.begin(); }

  Data_iterator end_env_data() { return m_env_data.end(); }

  /*! Set a data object to the face.
   * \param data The data object to set.
   */
  void set_env_data(const Data& data) {
    clear_env_data();
    add_env_data(data);
  }

  /*! Set a range of data objects to the face.
   * \param begin A begin iterator for the data range.
   * \param end A past-the-end iterator for the data range.
   */
  template <typename InputIterator>
  void set_env_data(InputIterator begin, InputIterator end) {
    clear_env_data();
    add_env_data(begin, end);
  }

  /*! set the data to be empty.
   */
  void set_no_env_data() {
    clear_env_data();
    m_is_env_set = true;
  }

   /*! Add a data object to the face.
   * \param data The additional data object.
   */
  void add_env_data(const Data& data) {
    m_env_data.push_back(data);
    m_is_env_set = true;
  }

  /*! Add a range of data objects to the face.
   * \param begin A begin iterator for the data range.
   * \param end A past-the-end iterator for the data range.
   */
  template <typename InputIterator>
  void add_env_data(InputIterator begin, InputIterator end) {
    for (auto it = begin; it != end; ++it) m_env_data.push_back(*it);
    m_is_env_set = true;
  }

  /*! Clear the data objects.
   */
  void clear_env_data() {
    m_env_data.clear();
    m_is_env_set = false;
  }

  /*! Check whether the set of data objects in the input range is equal to our
   * set of data objects
   */
  template <typename InputIterator>
  bool is_equal_env_data(InputIterator begin, InputIterator end) const {
    if (! is_env_set()) return false;

    // insert the input data objects into a set
    std::set<Data> input_data(begin, end);
    std::set<Data> my_data(begin_env_data(), end_env_data());
    if (input_data.size() != my_data.size()) return false;
    return (my_data == input_data);
  }

  template <typename InputIterator>
  bool has_equal_env_data(InputIterator begin, InputIterator end) const {
    if (! is_env_set()) return false;

    // insert the input data objects into a set
    std::set<Data> input_data(begin, end);
    std::set<Data> my_data(begin_env_data(), end_env_data());
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
  void set_aux_source(unsigned int id, HandleType h) {
    CGAL_precondition(id < 2);
    m_aux_source[id] = make_object(h);
  }

  void set_aux_source(unsigned int id, const Object& o) {
    CGAL_precondition(id < 2);
    CGAL_precondition(!o.is_empty());
    m_aux_source[id] = o;
  }

  const Object& aux_source(unsigned int id) {
    CGAL_precondition(id < 2);
    CGAL_precondition (!m_aux_source[id].is_empty());
    return m_aux_source[id];
  }

  /*! \brief returns true iff the point has been set already */
  bool aux_is_set(unsigned int id) const {
    CGAL_precondition(id < 2);
    return (! m_aux_source[id].is_empty());
  }
};

/*! Extend the planar-map vertex */
template <typename BaseVertex, class VertexData>
class Envelope_pm_vertex : public BaseVertex, public Dcel_info<VertexData> {
public:
  using Base_vertex = BaseVertex;
  using Vertex_data = VertexData;

private:
  using Self = Envelope_pm_vertex<Base_vertex, Vertex_data>;

protected:
  // all flags are bits in this variable:
  unsigned short flags;

  // the flags indications:
  enum Bit_pos {
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
  using Base_info = Dcel_info<Vertex_data>;

  /*! Constructor */
  Envelope_pm_vertex() : Dcel_info<Vertex_data>(), flags(0) {}

  /* void set_is_fake(bool b) { set_bit(IS_FAKE, b); }
   * bool is_fake() const { return get_bit(IS_FAKE); }
   */

  /* void set_is_intersection(bool b) { set_bit(IS_INTERSECTION, b); }
   * bool is_intersection() const { return get_bit(IS_FAKE); }
   */

  void set_is_equal_env_data_in_face(bool b) { set_bit(IS_EQUAL, b); }
  bool is_equal_env_data_in_face() const { return get_bit(IS_EQUAL); }

  void set_has_equal_env_data_in_face(bool b) { set_bit(HAS_EQUAL, b); }
  bool has_equal_env_data_in_face() const { return get_bit(HAS_EQUAL); }

  void set_is_equal_aux_data_in_face(unsigned int id, bool b) {
    CGAL_assertion(id < 2);
    set_bit(IS_EQUAL_AUX+id, b);
  }
  bool is_equal_aux_data_in_face(unsigned int id) const {
    CGAL_assertion(id < 2);
    return get_bit(IS_EQUAL_AUX+id);
  }

  void set_has_equal_aux_data_in_face(unsigned int id, bool b) {
    CGAL_assertion(id < 2);
    set_bit(HAS_EQUAL_AUX+id, b);
  }
  bool has_equal_aux_data_in_face(unsigned int id) const {
    CGAL_assertion(id < 2);
    return get_bit(HAS_EQUAL_AUX+id);
  }

  /*! Assign from another vertex.
   * \param v the other vertex.
   */
  virtual void assign(const Base_vertex& v) {
    Base_vertex::assign(v);

    const Self & ex_v = static_cast<const Self&>(v);
    this->Base_info::operator=(ex_v);
    flags = ex_v.flags;
  }

protected:
  void set_bit(unsigned int ind, bool b) {
    // set bit "ind" to 1 or 0:
    if (b) flags |= (1 << ind);
    else flags &= ~(1 << ind);
  }

  bool get_bit(unsigned int ind) const {
    // (1 << i) is bit i on, other bits off (start counting from 0)
    unsigned int mask = 1 << ind;
    return ((flags & mask) == mask);
  }
};

/*! Extend the planar-map halfedge */
template <typename BaseHalfedge, typename HalfedgeData>
class Envelope_pm_halfedge : public BaseHalfedge,
                             public Dcel_info<HalfedgeData> {
public:
  using Base_halfedge = BaseHalfedge;
  using Halfedge_data = HalfedgeData;

private:
  using Self = Envelope_pm_halfedge<Base_halfedge, Halfedge_data>;

protected:
  // all flags are bits in this variable:
  unsigned int flags;

  // flags indications
  enum Bit_pos {
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
  using Base_info = Dcel_info<Halfedge_data>;

  Envelope_pm_halfedge() : Dcel_info<Halfedge_data>(), flags(0) {}

 /* void set_is_fake(bool b) { set_bit(IS_FAKE, b); }
  * bool is_fake() const { return get_bit(IS_FAKE); }
  */

  void set_is_equal_env_data_in_face(bool b) { set_bit(IS_EQUAL_FACE, b); }
  bool is_equal_env_data_in_face() const { return get_bit(IS_EQUAL_FACE); }

  void set_has_equal_env_data_in_face(bool b) { set_bit(HAS_EQUAL_FACE, b); }
  bool has_equal_env_data_in_face() const { return get_bit(HAS_EQUAL_FACE); }

  void set_is_equal_aux_data_in_face(unsigned int id, bool b) {
    CGAL_assertion(id < 2);
    set_bit(IS_EQUAL_AUX_FACE+id, b);
  }
  bool is_equal_aux_data_in_face(unsigned int id) const {
    CGAL_assertion(id < 2);
    return get_bit(IS_EQUAL_AUX_FACE+id);
  }

  void set_has_equal_aux_data_in_face(unsigned int id, bool b) {
    CGAL_assertion(id < 2);
    set_bit(HAS_EQUAL_AUX_FACE+id, b);
  }
  bool has_equal_aux_data_in_face(unsigned int id) const {
    CGAL_assertion(id < 2);
    return get_bit(HAS_EQUAL_AUX_FACE+id);
  }

  void set_is_equal_env_data_in_target(bool b) { set_bit(IS_EQUAL_TARGET, b); }
  bool is_equal_env_data_in_target() const { return get_bit(IS_EQUAL_TARGET); }

  void set_has_equal_env_data_in_target(bool b) { set_bit(HAS_EQUAL_TARGET, b); }
  bool has_equal_env_data_in_target() const { return get_bit(HAS_EQUAL_TARGET); }

  void set_is_equal_aux_data_in_target(unsigned int id, bool b) {
    CGAL_assertion(id < 2);
    set_bit(IS_EQUAL_AUX_TARGET+id, b);
  }
  bool is_equal_aux_data_in_target(unsigned int id) const {
    CGAL_assertion(id < 2);
    return get_bit(IS_EQUAL_AUX_TARGET+id);
  }

  void set_has_equal_aux_data_in_target(unsigned int id, bool b) {
    CGAL_assertion(id < 2);
    set_bit(HAS_EQUAL_AUX_TARGET+id, b);
  }
  bool has_equal_aux_data_in_target(unsigned int id) const {
    CGAL_assertion(id < 2);
    return get_bit(HAS_EQUAL_AUX_TARGET+id);
  }
  // access to flags that contain relation between target and face
  void set_has_equal_env_data_in_target_and_face(bool b)
  { set_bit(HAS_EQUAL_F_T, b); }
  bool has_equal_env_data_in_target_and_face() const
  { return get_bit(HAS_EQUAL_F_T); }

  void set_has_equal_aux_data_in_target_and_face(unsigned int id, bool b) {
    CGAL_assertion(id < 2);
    set_bit(HAS_EQUAL_AUX_F_T+id, b);
  }
  bool has_equal_aux_data_in_target_and_face(unsigned int id) const {
    CGAL_assertion(id < 2);
    return get_bit(HAS_EQUAL_AUX_F_T+id);
  }

  /*! Assign from another halfedge.
   * \param h the other halfedge.
   */
  virtual void assign(const Base_halfedge& h) {
    Base_halfedge::assign(h);

    const Self & ex_h = static_cast<const Self&>(h);
    this->Base_info::operator=(ex_h);
    flags = ex_h.flags;
  }

protected:
  void set_bit(unsigned int ind, bool b) {
    if (b) flags |= (1 << ind);
    else flags &= ~(1 << ind);
    CGAL_assertion(get_bit(ind) == b);
  }

  bool get_bit(unsigned int ind) const {
    // (1 << i) is bit i on, other bits off (start counting from 0)
    unsigned int mask = 1 << ind;
    return ((flags & mask) == mask);
  }
};

  //! Extend the planar-map face.
template <typename BaseFace, typename FaceData>
class Envelope_pm_face : public BaseFace, public Dcel_info<FaceData> {
public:
  using Base_face = BaseFace;
  using Face_data = FaceData;

private:
  using Self = Envelope_pm_face<Base_face, Face_data>;

public:
  using Base_info = Dcel_info<Face_data>;
  using Data_container = std::list<Face_data>;
  using Data_iterator = typename Data_container::iterator;
  using Data_const_iterator = typename Data_container::const_iterator;

  /*! Construct.
   */
  Envelope_pm_face() : Dcel_info<Face_data>() {}

  /*! Assign from another face.
   * \param f the other face.
   */
  virtual void assign(const Base_face& f) {
    Base_face::assign(f);

    const Self& ex_f = static_cast<const Self&>(f);
    this->Base_info::operator=(ex_f);
  }
};

/*! A new dcel builder with full Envelope features */
template <typename Traits_, typename DcelData,
          typename VertexBase = Arr_vertex_base<typename Traits_::Point_2>,
          typename HalfedgeBase =
            Arr_halfedge_base<typename Traits_::X_monotone_curve_2>,
          typename FaceBase = Arr_face_base>
class Envelope_pm_dcel :
  public CGAL::Arr_dcel_base<Envelope_pm_vertex<VertexBase, DcelData>,
                             Envelope_pm_halfedge<HalfedgeBase, DcelData>,
                             Envelope_pm_face<FaceBase, DcelData>> {
public:
  using Dcel_data = DcelData;
  using Vertex_base = VertexBase;
  using Halfedge_base = HalfedgeBase;
  using Face_base = FaceBase;

  using Face_data = Dcel_data;
  using Env_pm_face = Envelope_pm_face<Face_base, Dcel_data>;
  using Face_data_iterator = typename Env_pm_face::Data_iterator;
  using Face_data_const_iterator = typename Env_pm_face::Data_const_iterator;

  using Edge_data = Dcel_data;
  using Edge_data_iterator = Face_data_iterator;
  using Edge_data_const_iterator = Face_data_const_iterator;

  using Vertex_data = Dcel_data;
  using Vertex_data_iterator = Face_data_iterator;
  using Vertex_data_const_iterator = Face_data_const_iterator;

  using Dcel_elem_with_data = Dcel_info<Dcel_data>;
  using Dcel_data_iterator = Face_data_iterator;
  using Dcel_data_const_iterator = Face_data_const_iterator;

  /*! \struct
   * An auxiliary structure for rebinding the DCEL with a new traits class.
   */
  template <typename T>
  struct rebind {
    using Pnt = typename T::Point_2;
    using Xcv = typename T::X_monotone_curve_2;
    using Rebind_v = typename Vertex_base::template rebind<Pnt>;
    using V_other = typename Rebind_v::other;
    using Rebind_h = typename Halfedge_base::template rebind<Xcv>;
    using H_other = typename Rebind_h::other;

    typedef Envelope_pm_dcel<T, Dcel_data, V_other, H_other, Face_base> other;
  };

  /*! Constructor */
  Envelope_pm_dcel() {}
};

} // namespace Envelope_3

} // namespace CGAL

#endif
