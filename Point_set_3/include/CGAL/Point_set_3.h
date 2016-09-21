// Copyright (c) 2016  GeometryFactory Sarl (France).
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
//
// Author(s)     : Simon Giraudot

#ifndef CGAL_POINT_SET_3_H
#define CGAL_POINT_SET_3_H

#include <stack>

#include <CGAL/Surface_mesh/Properties.h>

namespace CGAL {



/*!

  \ingroup PkgPointSet3

  \brief A collection of points with dynamically associated
  properties.

  This class provides the user with a flexible way to store and access
  a point set:

  - it can embed an arbitrary number of additional attributes such as
    normal vectors, colors, indices, etc.;

  - all functions of the package \ref PkgPointSetProcessing are
    provided with an overload that take a `Point_set_3` object as an
    argument.

  \tparam Point Point type.
  \tparam Vector Normal vector type.
 */

template <typename Point,
          typename Vector = typename Kernel_traits<Point>::Kernel::Vector_3>
class Point_set_3
{
public:

  /// \cond SKIP_IN_MANUAL
  typedef Point_set_3<Point, Vector> Point_set;
  /// \endcond
  
  typedef typename std::size_t Index; ///< Items are indices

  /// \cond SKIP_IN_MANUAL
  typedef typename Properties::Property_container<Index> Base;
  /// \endcond
  
  /*!
    \brief Property map used to associate attributes to the items of the point set.

    This class is a model of `LValuePropertyMap`.

    \tparam Type The type of the property.
  */
  template <class Type>
  struct Property_map
#ifndef DOXYGEN_RUNNING
    : public Properties::Property_map<Index, Type>
#endif
  {
  };

  /// \cond SKIP_IN_MANUAL
  typedef Property_map<std::size_t> Index_map;
  /// \endcond
  typedef Property_map<Point> Point_map; ///< Property map of points
  typedef Property_map<Vector> Vector_map; ///< Property map of vectors

  typedef typename Index_map::iterator iterator; ///< Iterator type of the point set
  typedef typename Index_map::const_iterator const_iterator; ///< Constant iterator type of the point set

  /*!
    \brief Class used to insert elements by defining the value of one of its properties.

    This class is a model of `OutputIterator`.

    \tparam Property The object `Property_map<Type>` that will be filled by the output iteration.
  */
  template <typename Property>
  class Property_back_inserter {
    /// \cond SKIP_IN_MANUAL
  public:
    typedef std::output_iterator_tag iterator_category;
    typedef typename Property::value_type value_type;
    typedef std::ptrdiff_t           difference_type;
    typedef void                     pointer;
    typedef void                     reference;

  private:

    Point_set* ps;
    Property* prop;
    std::size_t ind;
  
  public:
    
    Property_back_inserter(Point_set* ps, Property* prop, std::size_t ind=0)
      : ps(ps), prop (prop), ind(ind) {}
    Property_back_inserter& operator++() { return *this; }
    Property_back_inserter& operator++(int) { return *this; }
    Property_back_inserter& operator*() { return *this; }
    Property_back_inserter& operator= (const value_type& p)
    {
      if(ps->size() <= (typename Point_set::Index(ind)))
        ps->add_item();
      put(*prop, Point_set::Index(ind), p);
      ++ ind;
      return *this;
    }
    /// \endcond      
  };

  /*!
    \brief Property map that pushes a new item to the point set if needed.

    This class is a model of `WritablePropertyMap`.

    \tparam Property The object `Property_map<Type>` where to push new values.
  */
  template <typename Property>
  class Push_property_map
  {
    /// \cond SKIP_IN_MANUAL
  public:
    typedef std::size_t key_type;
    typedef typename Property::value_type value_type;
    typedef value_type& reference;
    typedef boost::lvalue_property_map_tag category;
    
    Point_set* ps;
    Property* prop;
    mutable std::size_t ind;

    Push_property_map(Point_set* ps = NULL,
              Property* prop = NULL,
              std::size_t ind=0)
      : ps(ps), prop(prop), ind(ind) {}

    friend void put(const Push_property_map& pm, std::size_t& i, const typename Property::value_type& t)
    {
      if(pm.ps->size() <= (pm.ind))
        pm.ps->add_item();
      put(*(pm.prop), Point_set::Index(pm.ind), t);
      i = pm.ind;
      ++pm.ind;
    }

    friend const reference get (const Push_property_map& pm, const std::size_t& i)
    {
      return ((*(pm.prop))[i]);
    }
    /// \endcond
  };

  typedef Property_back_inserter<Index_map> Index_back_inserter; ///< Back inserter on indices
  typedef Property_back_inserter<Point_map> Point_back_inserter; ///< Back inserter on points
  typedef Push_property_map<Point_map> Point_push_map; ///< Property map for pushing new points
  typedef Push_property_map<Vector_map> Vector_push_map; ///< Property map for pushing new vectors


protected:

  /// \cond SKIP_IN_MANUAL
  Base m_base;
  Index_map m_indices;
  Point_map m_points;
  Vector_map m_normals;
  std::size_t m_nb_removed;

  // Assignment operator not implemented and declared private to make
  // sure nobody uses the default one without knowing it
  Point_set_3& operator= (const Point_set_3&)
  {
    return *this;
  }

  // Copy constructor not implemented and declared private to make
  // sure nobody uses the default one without knowing it
  Point_set_3 (const Point_set_3& p)
  {
  }
  /// \endcond
  
public:

  /// \name Constructor
  /// @{


  /*!
    \brief Create an empty point set with no additional property.
   */
  Point_set_3 () : m_base()
  {
    clear();
  }

  /// @}

  /// \name Accessors and modifiers
  /// @{
  
  /*!
    \brief Clears the point set attributes and content.

    After calling this function, the object is the same as a newly
    constructed object. The additional attributes (such as normal vectors)
    are removed too and must be added again if they are still needed.
   */
  void clear()
  {
    m_base.clear();
    Properties::Property_map<Index,Index> pm = m_base.template add<std::size_t> ("index", (std::size_t)(-1)).first;
    m_indices = reinterpret_cast<Index_map&>(pm);
    Properties::Property_map<Index,Point> pm2 = m_base.template add<Point> ("point", Point (0., 0., 0.)).first;
    m_points = reinterpret_cast<Point_map&>(pm2);
    m_nb_removed = 0;
  }

  /*!
    \brief Memory management: reserve size to make the following insertions quicker.

    \param s Expected number of element to be inserted next.

    \note This method does not change the content of the point set and
    is only used for optimization.
   */
  void reserve (std::size_t s) { m_base.reserve (s); }
  /*!
    \brief Change size of the point set.

    \param s Target size of the point set.

    \note If `s` is smaller than the current size, the last elements
    of the point set are removed. If it is bigger, new elements are
    added with default values.
   */
  void resize (std::size_t s)
  {
    m_base.resize (s);
    if (s > size ())
      m_nb_removed += s;
    else if (m_base.size() - s > m_nb_removed)
      m_nb_removed = 0;
    else
      m_nb_removed -= s;
  }

  /// \cond SKIP_IN_MANUAL
  iterator add_item ()
  {
    m_base.push_back();
    m_indices[size()-1] = size()-1;
    return m_indices.end() - 1;
  }
  /// \endcond

  /*!
    \brief Convenience function to add a point.

    \param p Point to insert

    \note Properties of the added point are initialized to their
    default value.
   */
  iterator insert (const Point& p)
  {
    iterator out = add_item();
    m_points[size()-1] = p;
    return out;
  }

  /*!
    \brief Convenience function to add a point with a normal vector.

    \param p Point to insert
    \param n Associated normal vector

    \note Properties of the added point other than its normal vector
    are initialized to their default value.

    \note A normal property must have been added to the point set
    before using this method.
   */
  iterator insert (const Point& p, const Vector& n)
  {
    iterator out = insert (p);
    assert (has_normals());
    m_normals[size()-1] = n;
    return out;
  }

  /*!
    \brief Return the begin iterator.
  */
  iterator begin() { return m_indices.begin(); }
  /*!
    \brief Return the past-the-end iterator.
    \note The returned value is the same as `removed_begin()`.
  */
  iterator end() { return m_indices.end() - m_nb_removed; }
  /*!
    \brief Return the begin constant iterator.
  */
  const_iterator begin() const { return m_indices.begin(); }
  /*!
    \brief Return the past-the-end constant iterator.
    \note The returned value is the same as `removed_begin()`.
  */
  const_iterator end() const { return m_indices.end() - m_nb_removed; }
  /*!
    \brief Returns `true` if the number of non-removed element is 0.

    \note This does not count the removed elements.
  */
  bool empty() const { return (m_base.size() == m_nb_removed); }
  /*!
    \brief Returns the number of elements (not counting removed one).

    \note See `removed_size()` for getting the number of removed elements.
  */
  std::size_t size () const { return m_base.size() - m_nb_removed; }

  /*!
    \brief Get a reference to the wanted indexed point.

    \param index Index of the wanted point.
  */
  Point& point (Index index) { return m_points[m_indices[index]]; }
  /*!
    \brief Get a constant reference to the wanted indexed point.

    \param index Index of the wanted point.
  */
  const Point& point (Index index) const { return m_points[m_indices[index]]; }
  /*!
    \brief Get a reference to the wanted indexed point.

    \param it Iterator of the wanted item.
  */
  Point& point (iterator it) { return m_points[*it]; }
  /*!
    \brief Get a constant reference to the wanted indexed point.

    \param it Iterator of the wanted item.
  */
  const Point& point (const_iterator it) const { return m_points[*it]; }
  /*!
    \brief Get a reference to the wanted indexed point (convenience method).

    \param index Index of the wanted point.
  */
  Point& operator[] (Index index) { return m_points[m_indices[index]]; }
  /*!
    \brief Get a const reference the wanted indexed point (convenience method).

    \param index Index of the wanted point.
  */
  const Point& operator[] (Index index) const { return m_points[m_indices[index]]; }
  
  /*!
    \brief Get a reference to the wanted indexed normal.

    \param index Index of the wanted normal.

    \note The normal property must have been added to the point set
    before calling this method (see `add_normal_property()`).
  */
  Vector& normal (Index index) { return m_normals[m_indices[index]]; }
  /*!
    \brief Get a constant reference to the wanted indexed normal.

    \param index Index of the wanted normal.

    \note The normal property must have been added to the point set
    before calling this method (see `add_normal_property()`).
  */
  const Vector& normal (Index index) const { return m_normals[m_indices[index]]; }
  /*!
    \brief Get a reference to the wanted indexed normal.

    \param it Iterator of the wanted item.

    \note The normal property must have been added to the point set
    before calling this method (see `add_normal_property()`).
  */
  Vector& normal (iterator it) { return m_normals[*it]; }
  /*!
    \brief Get a constant reference to the wanted indexed normal.

    \param it Iterator of the wanted item.

    \note The normal property must have been added to the point set
    before calling this method (see `add_normal_property()`).
  */
  const Vector& normal (const_iterator it) const { return m_normals[*it]; }

  /// @}

  /// \name Garbage management
  /// @{

  /*!
    \brief Mark all elements between `it` and the last element as removed.

    \note The elements are just marked as removed and are not erased
    from the memory. `collect_garbage()` should be called if the
    memory needs to be freed.
  */
  void remove_from (iterator it)
  {
    m_nb_removed = static_cast<std::size_t>(std::distance (it, removed_end()));
  }
  
  /*!
    \brief Mark element as removed.

    \note The element is just marked as removed and is not erased from
    the memory. `collect_garbage()` should be called if the memory
    needs to be freed.
  */
  void remove (iterator it)
  {
    std::swap (*it, *(removed_begin() - 1));
    ++ m_nb_removed;
  }

  /*!
    \brief Iterator to the first removed element (equal to `removed_end()` if
    no elements are marked as removed.
  */
  iterator removed_begin () { return m_indices.end() - m_nb_removed; }
  /*!
    \brief Past-the-end iterator of the removed elements.
  */
  iterator removed_end () { return m_indices.end(); }
  /*!
    \brief Constant iterator to the first removed element (equal to `removed_end()` if
    no elements are marked as removed.
  */
  const_iterator removed_begin () const { return m_indices.end() - m_nb_removed; }
  /*!
    \brief Past-the-end constant iterator of the removed elements.
  */
  const_iterator removed_end () const { return m_indices.end(); }
  /*!
    \brief Number of removed elements.
  */
  std::size_t removed_size () const { return m_nb_removed; }

  /*!
    \brief Returns `true` if there are still removed elements in memory.
  */
  bool has_garbage () const { return (m_nb_removed != 0); }  

  /*!
    \brief Erase from memory the elements marked as removed.
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

    // for (std::size_t i = 0; i < 10; ++ i)
    //   std::cerr << m_indices[i] << " ";
    // std::cerr << std::endl;

    // Sorting based on the indices reorders the point set correctly
    quick_sort_on_indices ((std::ptrdiff_t)0, (std::ptrdiff_t)(m_base.size() - 1));

    // for (std::size_t i = 0; i < 10; ++ i)
    //   std::cerr << m_indices[i] << " ";
    // std::cerr << std::endl;

    m_base.resize (size ());
    m_base.shrink_to_fit ();
    m_nb_removed = 0;
  }

  /// @}


  /// \name Properties
  /// @{

  /*!
    \brief Test if property `name` of type `T` already exists.

    \tparam T type of the property.

    \param name Name of the property.
  */
  template <typename T>
  bool has_property (const std::string& name) const
  {
    std::pair<typename Properties::template Property_map<Index, T>, bool>
      pm = m_base.template get<T> (name);
    return pm.second;
  }
  
  /*!
    \brief Adds a new property `name` of type `T` with given default value.

    \tparam T type of the property.

    \param name Name of the property.

    \param t Value taken by the property on already created elements.

    \return Returns a pair containing the property map and a Boolean
    that is `true` if the property was added and `false` if it already
    exists (and was therefore not added but only returned).
  */
  template <class T>
  std::pair<Property_map<T>, bool>
  add_property (const std::string& name, const T t=T())
  {
    Properties::Property_map<Index,T> pm;
    bool added = false;
    boost::tie (pm, added) = m_base.template add<T> (name, t);
    return std::make_pair (reinterpret_cast<Property_map<T>&>(pm), added);
  }
  
  /*!
    \brief Gets the property `name` of type `T`.

    \tparam T type of the property.

    \param name Name of the property.

    \return Returns a pair containing: the wanted property map and a
    Boolean set to `true` or an empty property map and a Boolean set
    to `false` (if the property was not found).
  */
  template <class T> 
  std::pair<Property_map<T>,bool>
  property (const std::string& name) const
  {
    Properties::Property_map<Index,T> pm;
    bool okay = false;
    boost::tie (pm, okay) = m_base.template get<T>(name);
    return std::make_pair (reinterpret_cast<Property_map<T>&>(pm), okay);
  }

  /*!
    \brief Removes the wanted property.

    \tparam T type of the property.

    \param prop The property.

    \return Returns `true` if the property was removed and `false` if
    the property was not found.
  */
  template <class T> 
  bool remove_property (Property_map<T>& prop)
  {
    return m_base.remove (reinterpret_cast<Properties::Property_map<Index,T>&>(prop));
  }

  /*!
    \brief Convenience method that tests if the point set has normals.

    This method tests if a property of type `CGAL::Vector_3<Gt>` and
    named `normal` exists.
  */
  bool has_normals() const
  {
    std::pair<Vector_map, bool> pm = this->property<Vector> ("normal");
    return pm.second;
  }
  /*!
    \brief Convenience method that adds a normal property.

    This method adds a property of type `CGAL::Vector_3<Gt>` and named
    `normal`.

    \return `true` if the property was added, `false` if it already
    existed.
  */
  bool add_normal_property(const Vector& default_value = Vector(0., 0., 0.))
  {
    bool out = false;
    boost::tie (m_normals, out) = this->add_property<Vector> ("normal", default_value);
    return out;
  }
  /*!
    \brief Convenience method that removes the normal property.

    \return Returns `true` if the property was removed and `false` if
    the property was not found.
  */
  bool remove_normal_property()
  {
    return m_base.remove (m_normals);
  }

  /// @}


  /// \name Property Maps and Inserters
  /// @{

  /*!
    \brief Gets the push property map of the given property.

    \tparam T type of the property.

    \param prop The property map.

    \return Returns a pair containing: the wanted property map and a
    Boolean set to `true` or an empty property map and a Boolean set
    to `false` (if the property was not found).
  */
  template <class T>
  Push_property_map<Property_map<T> >
  push_property_map (Property_map<T>& prop)
  {
    return Push_property_map<Property_map<T> > (this, &prop, size());
  }

  /*!
    \brief Get the property map of the point attribute.
  */
  Point_map point_map()
  {
    return m_points;
  }

  /*!
    \brief Get the property map of the point attribute (constant version).
  */
  const Point_map point_map() const
  {
    return m_points;
  }

  /*!
    \brief Get the push property map of the point attribute.
  */
  Point_push_map point_push_map ()
  {
    return Point_push_map (this, &m_points, size());
  }

  /*!
    \brief Get the property map of the normal attribute.

    \note The normal property must have been added to the point set
    before calling this method (see `add_normal_property()`).
  */
  Vector_map normal_map ()
  {
    return m_normals;
  }
  /*!
    \brief Get the property map of the normal attribute (constant version).

    \note The normal property must have been added to the point set
    before calling this method (see `add_normal_property()`).
  */
  const Vector_map normal_map () const
  {
    return m_normals;
  }

  /*!
    \brief Get the push property map of the normal attribute.

    \note The normal property must have been added to the point set
    before calling this method (see `add_normal_property()`).
  */
  Vector_push_map normal_push_map ()
  {
    return Vector_push_map (this, &m_normals, size());
  }
  /*!
    \brief Get the back inserter on the index attribute.
  */
  Index_back_inserter index_back_inserter ()
  {
    return Index_back_inserter (this, &m_indices, size());
  }
  /*!
    \brief Get the back inserter on the point attribute.
  */
  Point_back_inserter point_back_inserter ()
  {
    return Point_back_inserter (this, &m_points, size());
  }


  /// @}
    

  /// \name Miscellaneous
  /// @{
  /*!
    \brief List properties with their types in a `std::string` object.
  */
  std::string properties() const
  {
    std::ostringstream oss;
    oss << "CGAL::Point_set_3<" << boost::core::demangle(typeid(Point).name())
        << "> with " << size() << " point(s) ("
        << removed_size() << " removed point(s) waiting to be deleted)" << std::endl;
    std::vector<std::string> prop = m_base.properties();
    for (std::size_t i = 0; i < prop.size(); ++ i)
      oss << " * \"" << prop[i] << "\" property of type "
          << boost::core::demangle(m_base.get_type(prop[i]).name()) << std::endl;

    return oss.str();
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

} // namespace CGAL


#endif // CGAL_POINT_SET_3_H
