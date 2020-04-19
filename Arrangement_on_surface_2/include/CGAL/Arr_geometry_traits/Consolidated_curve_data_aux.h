// Copyright (c) 2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Ron Wein          <wein@post.tau.ac.il>
//                 Efi Fogel         <efif@post.tau.ac.il>
#ifndef CGAL_CONSOLIDATED_CURVE_DATA_AUX_H
#define CGAL_CONSOLIDATED_CURVE_DATA_AUX_H

#include <CGAL/license/Arrangement_on_surface_2.h>


/*! \file
 * Definition of auxiliary classes for the usage of the
 * Arr_consolidated_curve_data_traits_2<Traits,Data> class.
 */

namespace CGAL {

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

/*! \struct
 * A functor for consolidating two unique lists.
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

} //namespace CGAL

#endif
