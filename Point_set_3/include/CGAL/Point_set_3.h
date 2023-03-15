// Copyright (c) 2016  GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Simon Giraudot

#ifndef CGAL_POINT_SET_3_H
#define CGAL_POINT_SET_3_H

#include <CGAL/license/Point_set_3.h>

#include <typeindex>

#include <CGAL/Surface_mesh/Properties.h>

#include <CGAL/Point_set_3/IO.h>

#include <CGAL/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/demangle.h>
#include <CGAL/assertions.h>

#include <algorithm>
#include <iterator>
#include <sstream>
#include <stack>
#include <utility>
#include <vector>

namespace CGAL {

template <typename Point,
          typename Vector = typename Kernel_traits<Point>::Kernel::Vector_3>
class Point_set_3;

/// \cond SKIP_IN_MANUAL
namespace internal {
  template<typename Point, typename Vector>
  class Point_set_3_index
  {
  public:
#ifdef CGAL_POINT_SET_3_USE_STD_SIZE_T_AS_SIZE_TYPE
    typedef std::size_t size_type;
#else
    typedef boost::uint32_t size_type;
#endif
    typedef CGAL::Point_set_3<Point, Vector> Point_set_3;

  private:
    friend class CGAL::Point_set_3<Point, Vector>;
    friend class Properties::Property_container<Point_set_3, Point_set_3_index>;
    template <class> friend class Properties::Property_array;
    template <class> friend struct Property_map;
    friend class std::vector<Point_set_3_index>;
    size_type value;

  public:
    Point_set_3_index(const Point_set_3_index& index) : value(static_cast<size_type>(index)) { }
    Point_set_3_index(const std::size_t& value) : value(static_cast<size_type>(value)) { }
    Point_set_3_index() : value(static_cast<size_type>(-1)) { }
    Point_set_3_index operator= (const Point_set_3_index& index) { value = index.value; return *this; }

    operator std::size_t() const { return static_cast<std::size_t>(value); }
    bool operator== (const Point_set_3_index& index) const { return value == index.value; }
    bool operator!= (const Point_set_3_index& index) const { return value != index.value; }
    bool operator<  (const Point_set_3_index& index) const { return value < index.value; }
    Point_set_3_index& operator++ () { ++value; return *this; }
    Point_set_3_index& operator-- () { --value; return *this; }
    Point_set_3_index operator++ (int) { Point_set_3_index tmp(*this); ++value; return tmp; }
    Point_set_3_index operator-- (int) { Point_set_3_index tmp(*this); --value; return tmp; }
};
/// \endcond
} // namespace internal

/*!

  \ingroup PkgPointSet3Ref

  \brief A collection of points with dynamically associated
  properties.

  An instance of this class stores a set of indices of type `Index`,
  each representing a point.  Properties can be associated to each
  point and can be retrieved using the index of the point.  There are
  two particular properties that are hard coded by this class: the
  coordinates of the points and the normal vectors.

  The coordinates of a point can be accessed using the index of the
  point and the member function `point()`. This property is always
  present. The normal vector of a point can be accessed using the
  index of the point and the `normal()` method. This property must
  be explicitly created.

  All properties can be accessed as a range using the methods
  `points()`, `normals()`, and `range()` for point coordinates,
  normal vectors, and other properties respectively.

  Removing a point with properties is achieved by moving its index
  at the end of the container and keeping track of the number of
  removed elements. A garbage collection method must be called to
  really remove it from memory.

  For convenience, all functions of the package \ref
  PkgPointSetProcessing3 automatically create the right named
  parameters if called with a `CGAL::Point_set_3` object as argument.


  \tparam Point Point type
  \tparam Vector Normal vector type

  \cgalModels `Range`
 */

template <typename Point,
          typename Vector>
class Point_set_3
{
public:

  /// \cond SKIP_IN_MANUAL
  typedef Point Point_type;
  typedef Vector Vector_type;
  typedef Point_set_3<Point, Vector> Point_set;

  using Index = internal::Point_set_3_index<Point, Vector>;

  typedef typename Properties::Property_container<Point_set, Index> Base;

  template <class Type>
  struct Property_map
    : public Properties::Property_map_base<Index, Type, Property_map<Type> >
  {
    typedef Properties::Property_map_base<Index, Type, Property_map<Type> > Base;
    Property_map() : Base() {}
    Property_map(const Base& pm): Base(pm) {}
  };
  typedef Property_map<Index> Index_map;

  template <typename Key, typename T>
  struct Get_property_map {
    typedef Property_map<T> type;
  };
  /// \endcond

#ifdef DOXYGEN_RUNNING
   /*!
  \brief This represents a point with associated properties.
  \cgalModels `::Index`
  \cgalModels `LessThanComparable`
  \cgalModels `Hashable`
  */
  class Index;
#endif

  typedef Point Point_3; ///< The point type
  typedef Vector Vector_3; ///< The vector type

#ifdef DOXYGEN_RUNNING
  typedef unspecified_type iterator; ///< Iterator type of the point set with value type `Index` \cgalModels RandomAccessIterator
  typedef unspecified_type const_iterator; ///< Constant iterator type of the point set with value type `Index` \cgalModels RandomAccessIterator
#else
  typedef typename Index_map::iterator iterator; ///< Iterator type of the point set
  typedef typename Index_map::const_iterator const_iterator; ///< Constant iterator type of the point set
#endif

  typedef Property_map<Point> Point_map; ///< Property map of points
  typedef Property_map<Vector> Vector_map; ///< Property map of vectors

  /// \cond SKIP_IN_MANUAL
  template <class Type>
  class Property_range
  {
  public:
    typedef CGAL::Property_map_to_unary_function<Property_map<Type> > Unary_function;
    typedef boost::transform_iterator<Unary_function,
                                      typename Point_set::iterator> iterator; // for NP helper
    typedef boost::transform_iterator<Unary_function,
                                      typename Point_set::const_iterator> const_iterator;

  private:
    const_iterator m_begin;
    const_iterator m_end;
    std::size_t m_size;

  public:
    Property_range (const Property_map<Type>& pmap,
                    typename Point_set::const_iterator begin,
                    typename Point_set::const_iterator end,
                    std::size_t size)
    {
      m_begin = boost::make_transform_iterator (begin, Unary_function(pmap));
      m_end = boost::make_transform_iterator (end, Unary_function(pmap));
      m_size = size;
    }

    const_iterator begin() const { return m_begin; }
    const_iterator end() const { return m_end; }
    std::size_t size() const { return m_size; }
    bool empty() const { return (m_size == 0); }
  };
  /// \endcond

  typedef Property_range<Point> Point_range; ///< Constant range of points
  typedef Property_range<Vector> Vector_range; ///< Constant range of vectors

protected:

  /// \cond SKIP_IN_MANUAL
  Base m_base;
  Index_map m_indices;
  Point_map m_points;
  Vector_map m_normals;
  std::size_t m_nb_removed;
  /// \endcond

public:

  /// \name Construction, Destruction, Assignment
  /// @{


  /*!
    \brief creates an empty point set with no additional property.

    \param with_normal_map `true` if the normal map should be
    added. If `false` (default value), the normal map can still be
    added later on (see `add_normal_map()`).
   */
  Point_set_3 (bool with_normal_map = false) : m_base()
  {
    clear();
    if (with_normal_map)
      add_normal_map();
  }

  /*!
    \brief Assignment operator, all properties with their content are copied.
   */
  Point_set_3& operator= (const Point_set_3& ps)
  {
    m_base = ps.m_base;
    m_indices = this->property_map<Index> ("index").first;
    m_points = this->property_map<Point> ("point").first;
    m_normals = this->property_map<Vector> ("normal").first;
    m_nb_removed = ps.m_nb_removed;
    return *this;
  }

  /// \cond SKIP_IN_MANUAL
  // copy constructor (same as assignment)
  Point_set_3 (const Point_set_3& ps)
  {
    m_base = ps.m_base;
    m_indices = this->property_map<Index> ("index").first;
    m_points = this->property_map<Point> ("point").first;
    m_normals = this->property_map<Vector> ("normal").first;
    m_nb_removed = ps.m_nb_removed;
  }
  /// \endcond

  /// @}

  /// \cond SKIP_IN_MANUAL
  const Base& base() const { return m_base; }
  /// \endcond


  /// \name Memory Management
  /// @{


  /*!

    \brief returns `true` if the number of elements not marked as
    removed is 0, `false` otherwise.

    \note This does not count the removed elements.

    \note The method `empty()` is also available (see `Range`) and
    does the same thing.
  */
  bool is_empty() const { return (m_base.size() == m_nb_removed); }
  /// \cond SKIP_IN_MANUAL
  bool empty() const { return is_empty(); }
  /// \endcond
  /*!
    \brief returns the number of elements (not counting elements marked as removed).

    \note See `number_of_removed_points()` for getting the number of elements marked as removed.

    \note The method `size()` is also available (see `Range`) and
    does the same thing.
  */
  std::size_t number_of_points () const { return m_base.size() - m_nb_removed; }
  /// \cond SKIP_IN_MANUAL
  std::size_t size () const { return number_of_points(); }
  /// \endcond

  /*!
    \brief merges `other` in the point set.

    Shifts the indices of points of `other` by `number_of_points() +
    other.number_of_points()`.

    Copies entries of all property maps which have the same name in
    the point set and `other`.  Property maps which are only in
    `other` are ignored.

    \note If `copy_properties()` with `other` as argument is called
    before calling this method, then all the content of `other` will
    be copied and no property will be lost in the process.

    \note Garbage is collected in both point sets when calling this method.
   */
  bool join (Point_set_3& other)
  {
    collect_garbage();
    other.collect_garbage();
    resize (number_of_points() + other.number_of_points());
    m_base.transfer (other.m_base);

    // Reset indices
    for (std::size_t i = 0; i < this->m_base.size(); ++ i)
      this->m_indices[i] = i;

    return true;
  }

  /*!
    \brief Clears the point set properties and content.

    After calling this method, the object is the same as a newly
    constructed object. The additional properties (such as normal
    vectors) are also removed and must thus be re-added if needed.
   */
  void clear()
  {
    m_base.clear();
    m_indices = this->add_property_map<Index>("index", typename Index::size_type(-1)).first;
    m_points = this->add_property_map<Point>("point", CGAL::ORIGIN).first;
    m_nb_removed = 0;
  }

  /*!
    \brief Clears all properties created.

    After calling this method, all properties are removed. The points
    are left unchanged.
   */
  void clear_properties()
  {
    Base other;
    other.template add<Index>("index", typename Index::size_type(-1));
    other.template add<Point>("point", CGAL::ORIGIN);
    other.resize(m_base.size());
    other.transfer(m_base);
    m_base.swap(other);
    m_indices = this->property_map<Index>("index").first;
    m_points = this->property_map<Point>("point").first;
  }

  /*!
    \brief increases the capacity of internal containers to be able to
    efficiently accommodate at least `s` elements

    \param s Expected final number of elements.

    \note This method does not change the content of the point set and
    is only used for optimization.
   */
  void reserve (std::size_t s) { m_base.reserve (s); }

  /*!
    \brief changes size of the point set.

    \param s Target size of the point set

    \note If the given size is larger than the current size, the
    capacity of the internal container is extended. If there are
    element marked as removed, they may be overwritten. If the given
    size is smaller than the current size, garbage is collected and
    the container is resized.
   */
  void resize (std::size_t s)
  {
    if (s < number_of_points() + number_of_removed_points())
      {
        collect_garbage();
        if (s < number_of_points())
          m_base.resize (s);
        else
          {
            std::size_t prev_s = number_of_points();
            m_base.resize (s);
            for (std::size_t i = prev_s; i < s; ++ i)
              m_indices[i] = i;
          }
      }
    else
      {
        std::size_t prev_s = number_of_points() + number_of_removed_points();
        m_base.resize (s);
        for (std::size_t i = prev_s; i < s; ++ i)
          m_indices[i] = i;
      }
  }

  /// @}

  /// \name Adding Points and Normals
  /// @{

  /*!
    \brief inserts a new element with default property values.

    \return The iterator on the newly added element.

    \note If a reallocation happens, all iterators, pointers and
    references related to the container are invalidated.  Otherwise,
    only the end iterator is invalidated, and all iterators, pointers
    and references to elements are guaranteed to keep referring to the
    same elements they were referring to before the call.
   */
  iterator insert ()
  {
    if (m_nb_removed == 0)
      {
        m_base.push_back();
        m_indices[size()-1] = size()-1;
        return m_indices.end() - 1;
      }
    else
      {
        -- m_nb_removed;
        iterator out = m_indices.end() - m_nb_removed - 1;
        Index idx = *out;
        m_base.reset(*out);
        *out = idx;
        return out;
      }
  }

  /*!
    \brief inserts new point with default property values.

    \param p Point to insert

    \note Properties of the added point are initialized to their
    default value.

    \note If a reallocation happens, all iterators, pointers and
    references related to the container are invalidated.  Otherwise,
    only the end iterator is invalidated, and all iterators, pointers
    and references to elements are guaranteed to keep referring to the
    same elements they were referring to before the call.

    \return The iterator on the newly added element.
   */
  iterator insert (const Point& p)
  {
    iterator out = insert();
    m_points[size()-1] = p;
    return out;
  }

  /*!
    \brief Convenience method to add a point with a normal vector.

    \param p Point to insert
    \param n Associated normal vector

    \note Properties of the added point other than its normal vector
    are initialized to their default value.

    \note If not already added, a normal property is automatically
    added to the point set when using this method. The default value
    for normal vectors is `CGAL::NULL_VECTOR`.

    \note If a reallocation happens, all iterators, pointers and
    references related to the container are invalidated.  Otherwise,
    only the end iterator is invalidated, and all iterators, pointers
    and references to elements are guaranteed to keep referring to the
    same elements they were referring to before the call.

    \return The iterator on the newly added element.
   */
  iterator insert (const Point& p, const Vector& n)
  {
    iterator out = insert (p);
    normal_map()[size()-1] = n;
    return out;
  }

  /*!
    \brief Convenience method to copy a point with all its properties
    from another point set.

    In the case where two point sets have the same properties, this
    method allows the user to easily copy one point (along with the
    values of all its properties) from one point set to another.

    \param other Point set to which the point to copy belongs
    \param idx Index of the point to copy in `other`

    \warning This point set and `other` must have the exact same
    properties, with the exact same names and types in the exact same
    order.

    \note If a reallocation happens, all iterators, pointers and
    references related to the container are invalidated.  Otherwise,
    only the end iterator is invalidated, and all iterators, pointers
    and references to elements are guaranteed to keep referring to the
    same elements they were referring to before the call.

    \return The iterator on the newly added element.
   */
  iterator insert (const Point_set_3& other, const Index& idx)
  {
    iterator out = insert();
    Index new_idx = *out;
    m_base.transfer(other.base(), idx, new_idx);
    *out = new_idx; // Do not copy index from other point set
    return out;
  }

  /// @}

  /// \name Accessors and Iterators
  /// @{

  /*!
    \brief returns the begin iterator.
  */
  iterator begin() { return m_indices.begin(); }
  /*!
    \brief returns the past-the-end iterator.
    \note The returned value is the same as `garbage_begin()`.
  */
  iterator end() { return m_indices.end() - m_nb_removed; }
  /*!
    \brief returns the begin constant iterator.
  */
  const_iterator begin() const { return m_indices.begin(); }
  /*!
    \brief returns the past-the-end constant iterator.
    \note The returned value is the same as `garbage_begin()`.
  */
  const_iterator end() const { return m_indices.end() - m_nb_removed; }
  /*!
    \brief returns a reference to the point corresponding to `index`.
  */
  Point& point (const Index& index) { return m_points[index]; }
  /*!
    \brief returns a constant reference to the point corresponding to `index`.
  */
  const Point& point (const Index& index) const { return m_points[index]; }
  /*!
    \brief returns a reference to the normal corresponding to `index`.

    \note If not already added, a normal property is automatically
    added to the point set (see `add_normal_map()`).
  */
  Vector& normal (const Index& index) { return normal_map()[index]; }
  /*!
    \brief returns a constant reference to the normal corresponding to `index`.

    \note If not already added, a normal property is automatically
    added to the point set (see `add_normal_map()`).
  */
  const Vector& normal (const Index& index) const { return normal_map()[index]; }

  /// @}

  /// \name Removal Methods
  /// @{

  /*!
    \brief marks all elements between `first` and `last` as removed.

    \note The elements are just marked as removed and are not erased
    from the memory. `collect_garbage()` should be called if the
    memory needs to be disallocated. Elements can be recovered with
    `cancel_removals()`.

    \note All iterators, pointers and references related to the
    container are invalidated.
  */
  void remove (iterator first, iterator last)
  {
    if (std::distance (last, end()) < 0)
      last = end();

    if (last == end())
      m_nb_removed += static_cast<std::size_t>(std::distance (first, end()));
    if (std::distance (first, end()) > 0)
      {
        iterator source = first;
        iterator dest = end() - 1;
        m_nb_removed += static_cast<std::size_t>(std::distance (first, last));
        while (source != last // All elements have been moved
               && dest != last - 1) // All elements are at the end of the container
          {
            std::cerr << "Swapping " << *source << " and " << *dest << std::endl;
            std::swap (*(source ++), *(dest --));
          }
      }
  }

  /// \cond SKIP_IN_MANUAL
  void remove_from (iterator first)
  {
    remove (first, end());
  }
  /// \endcond

  /*!
    \brief marks element specified by iterator as removed.

    \note The element is just marked as removed and is not erased from
    the memory. `collect_garbage()` should be called if the memory
    needs to be freed. The element can be recovered with
    `cancel_removals()`.

    \note The `end()` iterator is invalidated, `it` dereferences to a
    different element.
  */
  void remove (iterator it)
  {
    std::iter_swap (it, (end() - 1));
    ++ m_nb_removed;
  }

  /*!
    \brief marks element specified by `Index` as removed.

    \note The element is just marked as removed and is not erased from
    the memory. `collect_garbage()` should be called if the memory
    needs to be freed. The element can be recovered with
    `cancel_removals()`.

    \note The `end()` iterator is invalidated, `it` dereferences to a
    different element.
  */
  void remove (const Index& index)
  {
    iterator it = m_indices.begin() + index;
    while (*it != index)
      it = m_indices.begin() + *it;
    remove (it);
  }


  /// @}

  /// \name Garbage Management
  /// @{

  /*!
    \brief returns `true` if the element is marked as removed, `false`
    otherwise.

    \note When iterating between `begin()` and `end()`, no element
    marked as removed can be found.
  */
  bool is_removed (const_iterator it) const
  {
    return (std::distance (it, garbage_begin()) <= 0);
  }

  /*!
    \brief returns `true` if the element is marked as removed, `false`
    otherwise.

    \note When iterating between `begin()` and `end()`, no element
    marked as removed can be found.
  */
  bool is_removed (const Index& index) const
  {
    const_iterator it = m_indices.begin() + index;
    while (*it != index)
      it = m_indices.begin() + *it;
    return is_removed (it);
  }

  /*!
    \brief returns the constant iterator to the first element marked as removed
    (equal to `garbage_end()` if no elements are marked as removed.
  */
  const_iterator garbage_begin () const { return m_indices.end() - m_nb_removed; }
  /*!
    \brief returns the past-the-end constant iterator of the elements marked as removed.
  */
  const_iterator garbage_end () const { return m_indices.end(); }
  /*!
    \brief Number of removed points.
    \sa `garbage_size()`
  */
  std::size_t number_of_removed_points () const { return m_nb_removed; }
  /*!
    \brief Number of removed points.
    \sa `number_of_removed_points()`
  */
  std::size_t garbage_size () const { return number_of_removed_points(); }

  /*!  \brief returns `true` if there are elements marked as removed,
    `false` otherwise.
  */
  bool has_garbage () const { return (m_nb_removed != 0); }

  /*!
    \brief erases from memory the elements marked as removed.
  */
  void collect_garbage ()
  {
    // Indices indicate where to get the properties
    std::vector<std::size_t> indices (m_base.size());
    for (std::size_t i = 0; i < m_base.size(); ++ i)
      indices[m_indices[i]] = i;

    // Indices now indicate where to put the properties
    for (std::size_t i = 0; i < m_base.size(); ++ i)
      m_indices[i] = indices[i];

    // Sorting based on the indices reorders the point set correctly
    quick_sort_on_indices ((std::ptrdiff_t)0, (std::ptrdiff_t)(m_base.size() - 1));

    m_base.resize (size ());
    m_base.shrink_to_fit ();
    m_nb_removed = 0;
  }

  /*!
    \brief restores all removed points.

    After removing one or several points, calling this method restores
    the point set to its initial state: points that were removed (and
    their associated properties) are restored.

    \note This method is only guaranteed to work if no point was
    inserted after the removal: otherwise, some points might not be
    restored.

    \note If `collect_garbage()` was called after removal, the points
    are irremediably lost and nothing will be restored.
  */
  void cancel_removals()
  {
    m_nb_removed = 0;
    for (std::size_t i = 0; i < this->m_base.size(); ++ i)
      this->m_indices[i] = i;
  }

  /// @}


  /*! \name Property Handling

    A property `Property_map<Type>` allows to associate properties of
    type `Type` to a point. Properties can be added, looked up with a
    string and removed at runtime.
  */

  /// @{

#ifdef DOXYGEN_RUNNING
  /// Model of `LvaluePropertyMap` with `Index` as a key type and `Type`
  /// as value type.
  template <class Type>
  using Property_map = unspecified_type;
#endif


  /*!
    \brief tests whether property `name` of type `T` already exists.

    \tparam T type of the property.

    \param name Name of the property.
  */
  template <typename T>
  bool has_property_map (const std::string& name) const
  {
    std::pair<Property_map<T>, bool>
      pm = m_base.template get<T> (name);
    return pm.second;
  }

  /*!
    \brief adds a new property `name` of type `T` with given default value.

    \tparam T type of the property.

    \param name Name of the property.

    \param t Value taken by the property on already created elements.

    \return Returns a pair containing the property map and a Boolean
    that is `true` if the property was added and `false` if it already
    exists (and was therefore not added but only returned).
  */
  template <class T>
  std::pair<Property_map<T>, bool>
  add_property_map (const std::string& name, const T t=T())
  {
    Property_map<T> pm;
    bool added = false;
    std::tie (pm, added) = m_base.template add<T> (name, t);
    return std::make_pair (pm, added);
  }

  /*!
    \brief returns the property `name` of type `T`.

    \tparam T type of the property.

    \param name Name of the property.

    \return Returns a pair containing: the specified property map and a
    Boolean set to `true` or an empty property map and a Boolean set
    to `false` (if the property was not found).
  */
  template <class T>
  std::pair<Property_map<T>,bool>
  property_map (const std::string& name) const
  {
    Property_map<T> pm;
    bool okay = false;
    std::tie (pm, okay) = m_base.template get<T>(name);
    return std::make_pair (pm, okay);
  }

  /*!
    \brief removes the specified property.

    \tparam T type of the property.

    \param prop The property.

    \return Returns `true` if the property was removed and `false` if
    the property was not found.
  */
  template <class T>
  bool remove_property_map (Property_map<T>& prop)
  {
    return m_base.template remove<T> (prop);
  }

  /*!
    \brief Convenience method that tests whether the point set has normals.

    This method tests whether a property of type `Vector` and named
    `normal` exists.
  */
  bool has_normal_map() const
  {
    std::pair<Vector_map, bool> pm = this->property_map<Vector> ("normal");
    return pm.second;
  }
  /*!
    \brief Convenience method that adds a normal property.

    This method adds a property of type `Vector` and named
    `normal`.

    \return Returns a pair containing the normal map and a Boolean
    that is `true` if the property was added and `false` if it already
    exists (and was therefore not added but only returned).
  */
  std::pair<Vector_map, bool> add_normal_map (const Vector& default_value = CGAL::NULL_VECTOR)
  {
    bool out = false;
    std::tie (m_normals, out) = this->add_property_map<Vector> ("normal", default_value);
    return std::make_pair (m_normals, out);
  }
  /*!
    \brief returns the property map of the normal property.

    \note If the normal property has not been added yet to the point set
    before calling this method, the property map is automatically added
    with `add_normal_map()`.
  */
  Vector_map normal_map ()
  {
    if (!m_normals)
      add_normal_map();
    return m_normals;
  }
  /*!
    \brief returns the property map of the normal property (constant version).

    \note The normal property must have been added to the point set
    before calling this method (see `add_normal_map()`).
  */
  const Vector_map normal_map () const
  {
    return m_normals;
  }
  /*!
    \brief Convenience method that removes the normal property.

    \return Returns `true` if the property was removed and `false` if
    the property was not found.
  */
  bool remove_normal_map()
  {
    return m_base.template remove<Vector> (m_normals);
  }
  /*!
    \brief returns the property map of the point property.
  */
  Point_map point_map()
  {
    return m_points;
  }

  /*!
    \brief returns the property map of the point property (constant version).
  */
  const Point_map point_map() const
  {
    return m_points;
  }

  /*!
    \brief Copies the properties from another point set.

    All properties from `other` that do not already exist in this
    point set are added and filled to their default values. Properties
    that exist in both point sets are left unchanged.
  */
  void copy_properties (const Point_set_3& other)
  {
    m_base.copy_properties (other.base());

    m_normals = this->property_map<Vector> ("normal").first; // In case normal was added
  }


  /*!
    \brief returns a vector with all strings that describe properties.
  */
  std::vector<std::string> properties() const
  {
    std::vector<std::string> out = m_base.properties();
    out.erase (out.begin()); // remove "index"
    out.erase (out.begin()); // remove "point"
    return out;
  }

  /*!
    \brief returns a vector of pairs that describe properties and associated types.
  */
  std::vector<std::pair<std::string, std::type_index> > properties_and_types() const
  {
    std::vector<std::string> prop = m_base.properties();
    prop.erase (prop.begin()); // remove "index"
    prop.erase (prop.begin()); // remove "point"

    std::vector<std::pair<std::string, std::type_index> > out; out.reserve (prop.size());
    for (std::size_t i = 0; i < prop.size(); ++ i)
      out.push_back (std::make_pair (prop[i], std::type_index(m_base.get_type(prop[i]))));
    return out;
  }


  /*!
    \brief returns a sequence of \ref bgl_namedparameters "Named Parameters" to be used in Point Set Processing algorithms.

    The following named parameters are used:
    - `point_map`: contains the point property map (see `point_map()`)
    - `normal_map`: contains the normal map (see `normal_map()`)
    - `geom_traits`: contains the kernel `typename Kernel_traits<Point>`::`Kernel`

    \warning this method does not check if the normal map was
    instantiated or not. The normal map named parameter should not be
    used if this property was not instantiated first.
  */
#ifdef DOXYGEN_RUNNING
  unspecified_type
#else
  Named_function_parameters
  <typename Kernel_traits<Point>::Kernel,
   internal_np::geom_traits_t,
   Named_function_parameters
   <typename CGAL::Point_set_3<Point, Vector>::template Property_map<Vector>,
    internal_np::normal_t,
    Named_function_parameters
    <typename CGAL::Point_set_3<Point, Vector>::template Property_map<Point>,
     internal_np::point_t> > >
#endif
  inline parameters() const
  {
    return CGAL::parameters::point_map (m_points).
      normal_map (normal_map()).
      geom_traits (typename Kernel_traits<Point>::Kernel());
  }


  /// \cond SKIP_IN_MANUAL
  std::string info() const
  {
    std::ostringstream oss;
    oss << "CGAL::Point_set_3<" << CGAL::demangle(typeid(Point).name())
        << "> with " << size() << " point(s) ("
        << number_of_removed_points() << " removed point(s) waiting to be deleted)" << std::endl;
    std::vector<std::string> prop = m_base.properties();
    for (std::size_t i = 0; i < prop.size(); ++ i)
      oss << " * \"" << prop[i] << "\" property of type "
          << CGAL::demangle(m_base.get_type(prop[i]).name()) << std::endl;

    return oss.str();
  }
  /// \endcond
  /// @}

  /// \name Ranges
  /// @{

#ifdef DOXYGEN_RUNNING
  /// Model of `ConstRange` that handles constant ranges for property
  /// maps with value type `Type`.
  template <class Type>
  using Property_range = unspecified_type;
#endif

  /*!
    \brief returns a property as a range.
  */
  template <class T>
  Property_range<T> range (const Property_map<T>& pmap) const
  {
    return Property_range<T> (pmap, begin(), end(), number_of_points());
  }

  /*!
    \brief returns a constant range of points.
  */
  Point_range points () const
  {
    return this->range<Point> (m_points);
  }

  /*!
    \brief returns a constant range of normals.
  */
  Vector_range normals () const
  {
    return this->range<Vector> (normal_map());
  }

  /// @}

  /*!
    \name Push Property Maps and Inserters (Advanced)

    \cgalAdvancedBegin
    The following method are specifically designed to make
    `CGAL::Point_set_3` usable with \cgal input/output functions.
    \cgalAdvancedEnd
  */

  /// @{


#ifdef DOXYGEN_RUNNING
  /// \cgalAdvancedType
  /// \cgalAdvancedBegin
  /// Model of `OutputIterator` used to insert elements by defining
  /// the value of the property `Property`.
  /// \cgalAdvancedEnd
  template <class Property>
  using Property_back_inserter = unspecified_type;

  /// \cgalAdvancedType
  /// \cgalAdvancedBegin
  /// Model of `WritablePropertyMap` based on `Property` and that
  /// is allowed to push new items to the point set if needed.
  /// \cgalAdvancedEnd
  template <class Property>
  using Push_property_map = unspecified_type;
#endif

  /// \cond SKIP_IN_MANUAL
  template <typename Property>
  class Property_back_inserter {

  public:
    typedef std::output_iterator_tag iterator_category;
    typedef typename Property::value_type value_type;
    typedef std::ptrdiff_t           difference_type;
    typedef void                     pointer;
    typedef void                     reference;

  private:

    Point_set* ps;
    Property* prop;
    Index ind;

  public:

    Property_back_inserter(Point_set* ps, Property* prop, Index ind=Index())
      : ps(ps), prop (prop), ind(ind) {}
    Property_back_inserter& operator++() { return *this; }
    Property_back_inserter& operator++(int) { return *this; }
    Property_back_inserter& operator*() { return *this; }
    Property_back_inserter& operator= (const value_type& p)
    {
      if(ps->size() <= ind)
        ps->insert();
      put(*prop, ind, p);
      ++ ind;
      return *this;
    }

  };

  template <typename Property>
  class Push_property_map
  {

  public:
    typedef Index key_type;
    typedef typename Property::value_type value_type;
    typedef value_type& reference;
    typedef boost::read_write_property_map_tag category;

    Point_set* ps;
    Property* prop;
    mutable Index ind;

    Push_property_map(Point_set* ps = nullptr,
                      Property* prop = nullptr,
                      Index ind=Index())
      : ps(ps), prop(prop), ind(ind) {}

    friend void put(const Push_property_map& pm, Index& i, const value_type& t)
    {
      if(pm.ps->size() <= (pm.ind))
        pm.ps->insert();
      put(*(pm.prop), pm.ind, t);
      i = pm.ind;
      ++pm.ind;
    }

    friend reference get (const Push_property_map& pm, const Index& i)
    {
      return ((*(pm.prop))[i]);
    }

  };
  /// \endcond

  /// \cgalAdvancedType
  /// \cgalAdvancedBegin
  /// Back inserter on indices
  /// \cgalAdvancedEnd
  typedef Property_back_inserter<Index_map> Index_back_inserter;
  /// \cgalAdvancedType
  /// \cgalAdvancedBegin
  /// Back inserter on points
  /// \cgalAdvancedEnd
  typedef Property_back_inserter<Point_map> Point_back_inserter;
  /// \cgalAdvancedType
  /// \cgalAdvancedBegin
  /// Property map for pushing new points
  /// \cgalAdvancedEnd
  typedef Push_property_map<Point_map> Point_push_map;
  /// \cgalAdvancedType
  /// \cgalAdvancedBegin
  /// Property map for pushing new vectors
  /// \cgalAdvancedEnd
  typedef Push_property_map<Vector_map> Vector_push_map;

  /*!
    \cgalAdvancedFunction
    \cgalAdvancedBegin
    \brief returns the push property map of the given property.

    \tparam T type of the property.

    \param prop The property map.

    \return Returns a pair containing: the specified property map and a
    Boolean set to `true` or an empty property map and a Boolean set
    to `false` (if the property was not found).
    \cgalAdvancedEnd
  */
  template <class T>
  Push_property_map<Property_map<T> >
  push_property_map (Property_map<T>& prop)
  {
    return Push_property_map<Property_map<T> > (this, &prop, size());
  }
  /*!
    \cgalAdvancedFunction
    \cgalAdvancedBegin
    \brief returns the push property map of the point property.
    \cgalAdvancedEnd
  */
  Point_push_map point_push_map ()
  {
    return Point_push_map (this, &m_points, size());
  }
  /*!
    \cgalAdvancedFunction
    \cgalAdvancedBegin
    \brief returns the push property map of the normal property.

    \note The normal property must have been added to the point set
    before calling this method (see `add_normal_map()`).
    \cgalAdvancedEnd
  */
  Vector_push_map normal_push_map ()
  {
    return Vector_push_map (this, &m_normals, size());
  }
  /*!
    \cgalAdvancedFunction
    \cgalAdvancedBegin
    \brief returns the back inserter on the index property.
    \cgalAdvancedEnd
  */
  Index_back_inserter index_back_inserter ()
  {
    return Index_back_inserter (this, &m_indices, size());
  }
  /*!
    \cgalAdvancedFunction
    \cgalAdvancedBegin
    \brief returns the back inserter on the point property.
    \cgalAdvancedEnd
  */
  Point_back_inserter point_back_inserter ()
  {
    return Point_back_inserter (this, &m_points, size());
  }

  /// @}


private:
  /// \cond SKIP_IN_MANUAL
  void quick_sort_on_indices (std::ptrdiff_t begin, std::ptrdiff_t end)
  {
    std::stack<std::pair<std::ptrdiff_t, std::ptrdiff_t> >
      todo;
    todo.push (std::make_pair (begin, end));

    while (!(todo.empty()))
      {
        std::pair<std::ptrdiff_t, std::ptrdiff_t>
          current = todo.top();
        todo.pop();

        if (current.first < current.second)
          {
            std::ptrdiff_t p = current.first + (rand() % (current.second - current.first));
            p = quick_sort_partition (current.first, current.second, p);
            todo.push (std::make_pair (current.first, p-1));
            todo.push (std::make_pair (p+1, current.second));
          }
      }
  }

  std::ptrdiff_t quick_sort_partition (std::ptrdiff_t begin, std::ptrdiff_t end, std::ptrdiff_t p)
  {
    m_base.swap (p, end);
    std::ptrdiff_t j = begin;
    for (std::ptrdiff_t i = begin; i < end; ++ i)
      if (m_indices[i] <= m_indices[end])
        {
          m_base.swap (i, j);
          j ++;
        }
    m_base.swap (end, j);
    return j;
  }
  /// \endcond


}; // end of class Point_set_3

/*!

  \brief Append `other` at the end of `ps`.

  \relates Point_set_3

   Shifts the indices of points of `other` by `ps.number_of_points() +
   other.number_of_points()`.

   Copies entries of all property maps which have the same name in `ps` and `other`.
   Property maps which are only in `other` are ignored.

   \note Garbage is collected in both point sets when calling this method.

*/
template <typename Point, typename Vector>
Point_set_3<Point, Vector>& operator+=(Point_set_3<Point, Vector>& ps,
                                       Point_set_3<Point, Vector>& other)
{
  ps.join(other);
  return ps;
}

/// \cond SKIP_IN_MANUAL
// specialization for default named parameters
template <typename Point, typename Vector, typename NamedParameters, typename DPM, typename DVM>
struct Point_set_processing_3_np_helper<Point_set_3<Point, Vector>, NamedParameters, DPM, DVM>
{
  typedef typename std::iterator_traits<typename Point_set_3<Point, Vector>::iterator>::value_type Value_type;

  typedef typename Kernel_traits<Point>::Kernel Default_geom_traits;
  typedef typename Point_set_3<Point, Vector>::template Property_map<Vector> DefaultNMap;
  typedef typename Point_set_3<Point, Vector>::template Property_map<Point> DefaultPMap;
  typedef const typename Point_set_3<Point, Vector>::template Property_map<Point> DefaultConstPMap;

  typedef typename internal_np::Lookup_named_param_def<internal_np::point_t,
    NamedParameters,DefaultPMap> ::type  Point_map; // public
  typedef typename internal_np::Lookup_named_param_def<internal_np::point_t,
    NamedParameters,DefaultConstPMap> ::type  Const_point_map; // public

  typedef typename internal_np::Lookup_named_param_def <
      internal_np::geom_traits_t,
      NamedParameters,
      Default_geom_traits
    > ::type  Geom_traits; // public

  typedef typename Geom_traits::FT FT; // public

  typedef typename internal_np::Lookup_named_param_def<
    internal_np::normal_t,
    NamedParameters,
    DefaultNMap
    > ::type  Normal_map; // public

  static Point_map get_point_map(Point_set_3<Point, Vector>& ps, const NamedParameters& np)
  {
    return parameters::choose_parameter(parameters::get_parameter(np, internal_np::point_map), ps.point_map());
  }

  static Const_point_map get_const_point_map(const Point_set_3<Point, Vector>& ps, const NamedParameters& np)
  {
    return parameters::choose_parameter(parameters::get_parameter(np, internal_np::point_map), ps.point_map());
  }

  static const Normal_map get_normal_map(const Point_set_3<Point, Vector>& ps, const NamedParameters& np)
  {
    return parameters::choose_parameter(parameters::get_parameter(np, internal_np::normal_map), ps.normal_map());
  }

  static Normal_map get_normal_map(Point_set_3<Point, Vector>& ps, const NamedParameters& np)
  {
    return parameters::choose_parameter(parameters::get_parameter(np, internal_np::normal_map), ps.normal_map());
  }

  static Geom_traits get_geom_traits(const Point_set_3<Point, Vector>&, const NamedParameters& np)
  {
    return parameters::choose_parameter<Geom_traits>(parameters::get_parameter(np, internal_np::geom_traits));
  }

  static bool has_normal_map(const Point_set_3<Point, Vector>& ps, const NamedParameters&)
  {
    if (ps.has_normal_map())
      return true;
    using CGAL::parameters::is_default_parameter;
    return !(is_default_parameter<NamedParameters, internal_np::normal_t>::value);
  }

  static constexpr bool has_normal_map(Point_set_3<Point, Vector>&, const NamedParameters&)
  {
    return true; // either available in named parameters, and always available in Point_set_3 otherwise
  }
};
/// \endcond

} // namespace CGAL

namespace boost {
template <typename Point, typename Vector>
std::size_t hash_value(const typename CGAL::internal::Point_set_3_index<Point, Vector>& i)
{
  return i;
}
} // namespace boost

namespace std {
template <typename Point, typename Vector>
struct hash< CGAL::internal::Point_set_3_index<Point, Vector> > {
  std::size_t operator()(const CGAL::internal::Point_set_3_index<Point, Vector>& i) const
  {
    std::size_t ret = i;
    return ret;
  }
};
} // namespace std

#endif // CGAL_POINT_SET_3_H
