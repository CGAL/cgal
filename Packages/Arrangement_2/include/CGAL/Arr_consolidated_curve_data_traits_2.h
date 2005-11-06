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
#ifndef CGAL_ARR_CONSOLIDATED_CURVE_DATA_TRAITS_2_H
#define CGAL_ARR_CONSOLIDATED_CURVE_DATA_TRAITS_2_H

/*! \file
 * Definition of the Arr_consolidated_curve_data_traits_2<Traits,Data> class.
 */

#include<CGAL/Arr_curve_data_traits_2.h>

CGAL_BEGIN_NAMESPACE

/*! \class
 * Representation of a set of data objects (without duplicates), using a
 * simple list.
 */
template <class Data_>
class _Unique_list
{
public:
  typedef Data_                Data;
  typedef _Unique_list<Data>   Self;

  typedef typename std::list<Data>::const_iterator  const_iterator;

private:
  
  std::list<Data>     m_list;

public:

  /*! Default constructor. */
  _Unique_list () :
    m_list()
  {}

  /*! Construct a singleton list. */
  _Unique_list (const Data& data) :
    m_list ()
  {
    m_list.push_back (data);
  }

  /*! Go over the data objects in list. */
  const_iterator begin () const
  {
    return (m_list.begin());
  }

  const_iterator end () const
  {
    return (m_list.end());
  }

  /*! Get the list size. */
  std::size_t size () const
  {
    return (m_list.size());
  }

  /*! Get the first (or last) data object. */
  const Data& front () const
  {
    return (m_list.front());
  }

  const Data& back () const
  {
    return (m_list.back());
  }

  /*! Equality operator. */
  bool operator== (const Self& other) const
  {
    if (size() != other.size())
      return (false);

    const_iterator    iter;

    for (iter = begin(); iter != end(); ++iter)
    {
      if (other.find (*iter) == other.end())
        return (false);
    }

    for (iter = other.begin(); iter != other.end(); ++iter)
    {
      if (find (*iter) == end())
        return (false);
    }

    return (true);
  }

  /*!
   * Find the given data object is contained in the list.
   * \param data The data object.
   * \return An iterator for the data object, or end() if it is not found.
   */
  const_iterator find (const Data& data) const
  {
    const_iterator   iter = m_list.begin();

    while (iter != m_list.end())
    {
      if (*iter == data)
        break;
      ++iter;
    }
    return (iter);
  }
  
  /*!
   * Insert an object into the list. 
   * \param data The data object.
   * \return (true) if the data object has been successfully inserted;
   *         (false) otherwise (if it already exists).
   */
  bool insert (const Data& data)
  {
    if (find (data) != m_list.end())
      return (false);

    m_list.push_back (data);
    return (true);
  }
  
  /*!
   * Erase an object from the list. 
   * \param data The data object.
   * \return (true) if the data object has been successfully erased;
   *         (false) otherwise (if it is not in the list).
   */
  bool erase (const Data& data)
  {
    typename std::list<Data>::iterator  iter = m_list.begin();

    while (iter != m_list.end())
    {
      if (*iter == data)
      {
        // Erase the current data object.
        m_list.erase (iter);
        return (true);
      }
      ++iter;
    }

    // The data object is not found in the list:
    return (false);
  }

  /*! Clear the list. */
  void clear ()
  {
    m_list.clear();
    return;
  }
};

/*!
 * \struct A functor for consolidating two unique lists.
 */
template <class Data>
struct _Consolidate_unique_lists
{
  _Unique_list<Data> operator() (const _Unique_list<Data>& list1,
                                 const _Unique_list<Data>& list2) const
  {
    _Unique_list<Data>  result = list1;

    typename _Unique_list<Data>::const_iterator  iter;

    for (iter = list2.begin(); iter != list2.end(); ++iter)
      result.insert (*iter);

    return (result);
  }
};

/*! \class
 * A generic traits class for maintaining an arrangement of curves that have
 * an extra data field. This traits class is templated with a Data class an
 * an ordinary traits class which is also used as a based traits class to
 * inherit from. It extracts the original Curve_2 and X_monotone_curve_2 types
 * from the ordinary traits class, and redefines them to have Data as an extra
 * field in the Curve_2 type, and a container of Data objects for the extended
 * X_monotone_curve_2 type.
 * The Data field is updated when the curves are converted from Curve_2 to
 * X_monotone_curve_2, and when the X_monotone_curve_2 curves are split.
 * When two x-monotone curves overlap, their data containers are consolidated
 * and attached to the resulting subcurve.
 * All other functors are inherited from the base ordinary traits class.
 */
template <class Traits_, class Data_>
class Arr_consolidated_curve_data_traits_2 :
  public Arr_curve_data_traits_2<Traits_,
                                 _Unique_list<Data_>, 
                                 _Consolidate_unique_lists<Data_>,
                                 Data_>
{
public:

  typedef Traits_                                   Base_traits;
  typedef Data_                                     Data;
  typedef _Unique_list<Data_>                       Data_container;
  typedef typename Data_container::const_iterator   Data_iterator;
  typedef typename Data_container::const_iterator   Data_const_iterator;
  typedef typename Base_traits::Curve_2             Base_curve_2;
  typedef typename Base_traits::X_monotone_curve_2  Base_x_monotone_curve_2;
  typedef typename Base_traits::Point_2             Point_2;

  typedef typename Base_traits::Has_left_category   Has_left_category;
  typedef typename Base_traits::Has_merge_category  Base_has_merge_category;
  typedef Tag_true                                  Has_merge_category;

};

CGAL_END_NAMESPACE

#endif

