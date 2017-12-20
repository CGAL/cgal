// Copyright (c) 2005,2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
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
//
// Author(s)     : Ron Wein          <wein@post.tau.ac.il>

#ifndef CGAL_ARRANGEMENT_2_ITERATORS_H
#define CGAL_ARRANGEMENT_2_ITERATORS_H

#include <CGAL/license/Arrangement_on_surface_2.h>


#include <functional>

/*! \file
 * Definitions of auxiliary iterator adaptors.
 */

namespace CGAL {

/*!
 * \class
 * An iterator adaptor for dereferencing the value-type of the iterator class
 * (given as Iterator_), which is supposed to be a pointer, and handle it as
 * the value-type given by Value_.
 */
template <class Iterator_, class Value_, class Diff_, class Category_>
class I_Dereference_iterator
{
public:

  // Type definitions:
  typedef Iterator_                               Iterator;
  typedef I_Dereference_iterator<Iterator_,
				 Value_,
				 Diff_,
				 Category_>       Self;

  typedef Category_                               iterator_category;
  typedef Value_                                  value_type;
  typedef value_type&                             reference;
  typedef value_type*                             pointer;
  typedef Diff_                                   difference_type;

protected:

  Iterator        iter;           // The internal iterator.

public:

  /// \name Construction
  //@{
  I_Dereference_iterator ()
  {}

  I_Dereference_iterator (Iterator it) :
    iter(it)
  {}
  //@}

  /// \name Basic operations.
  //@{
  bool operator== (const Self& it) const
  {
    return (iter == it.iter);
  }
    
  bool operator!= (const Self& it) const
  {
    return (!(iter == it.iter));
  }
    
  Iterator current_iterator () const
  {
    return (iter);
  }
    
  pointer ptr () const
  {
    return (static_cast<value_type *> (*iter));
  }

  reference operator* () const
  {
    return (*(ptr()));
  }

  pointer operator-> () const
  {
    return (ptr());
  }
  //@}
    
  /// \name Incremernt operations (forward category).
  //@{
  Self& operator++()
  {
    ++iter;
    return (*this);
  }

  Self operator++ (int )
  {
    Self tmp = *this;
    ++iter;
    return (tmp);
  }
  //@}

  /// \name Decremernt operations (bidirectional category).
  //@{
  Self& operator-- ()
  {
    --iter;
    return (*this);
  }

  Self operator-- (int )
  {
    Self tmp = *this;
    --iter;
    return (tmp);
  }
  //@}
};

/*!
 * \class
 * An iterator adaptor for dereferencing the value-type of the const iterator
 * class (given as CIterator_), which is supposed to be a pointer, and handle
 * it as the value-type given by Value_.
 */
template <class CIterator_, class MIterator_, 
	  class Value_, class Diff_, class Category_>
class I_Dereference_const_iterator
{
public:

  // Type definitions:
  typedef CIterator_                              Const_iterator;
  typedef MIterator_                              Mutable_iterator;
  typedef I_Dereference_const_iterator<CIterator_,
				       MIterator_,
				       Value_,
				       Diff_,
				       Category_> Self;

  typedef Category_                               iterator_category;
  typedef Value_                                  value_type;
  typedef const value_type&                       reference;
  typedef const value_type*                       pointer;
  typedef Diff_                                   difference_type;

protected:

  Const_iterator        iter;           // The internal iterator.

public:

  /// \name Construction
  //@{
  I_Dereference_const_iterator ()
  {}

  I_Dereference_const_iterator (Const_iterator it) :
    iter(it)
  {}

  I_Dereference_const_iterator (Mutable_iterator it) :
    iter (Const_iterator (&(*it)))
  {}

  //@}

  /// \name Basic operations.
  //@{
  bool operator== (const Self& it) const
  {
    return (iter == it.iter);
  }
    
  bool operator!= (const Self& it) const
  {
    return (!(iter == it.iter));
  }
    
  Const_iterator current_iterator () const
  {
    return (iter);
  }
    
  pointer ptr () const
  {
    return (static_cast<const value_type *> (*iter));
  }

  reference operator* () const
  {
    return (*(ptr()));
  }

  pointer operator-> () const
  {
    return (ptr());
  }
  //@}
    
  /// \name Incremernt operations (forward category).
  //@{
  Self& operator++()
  {
    ++iter;
    return (*this);
  }

  Self operator++ (int )
  {
    Self tmp = *this;
    ++iter;
    return (tmp);
  }
  //@}

  /// \name Decremernt operations (bidirectional category).
  //@{
  Self& operator-- ()
  {
    --iter;
    return (*this);
  }

  Self operator-- (int )
  {
    Self tmp = *this;
    --iter;
    return (tmp);
  }
  //@}
};

/*!
 * \class
 * An iterator adaptor for the filtering a DCEL iterator (given as Iterator_)
 * using a given filter functor (Filter_).
 */
template <class Iterator_, class Filter_,
          class Value_, class Diff_, class Category_>
class I_Filtered_iterator
{
public:

  typedef Iterator_                       Iterator;
  typedef Filter_                         Filter; 
  typedef I_Filtered_iterator<Iterator_,
                              Filter_,
                              Value_,
                              Diff_,
                              Category_>  Self;

  typedef Category_                       iterator_category;
  typedef Value_                          value_type;
  typedef value_type&                     reference;
  typedef value_type*                     pointer;
  typedef Diff_                           difference_type;

protected:

  Iterator        nt;       // The internal iterator (this member should not
                            // be renamed in order to comply with the
                            // HalfedgeDS circulators that refer to it).
  Iterator        iend;     // A past-the-end iterator.        
  Filter          filt;     // The filter functor.

public:

  /*! Constructors. */
  I_Filtered_iterator()
  {}

  I_Filtered_iterator (Iterator it) :
    nt (it),
    iend (nt)
  {}

  I_Filtered_iterator (Iterator it, Iterator end) :
    nt (it),
    iend (end)
  {
    while (nt != iend && ! filt (*nt))
      ++nt;
  }

  I_Filtered_iterator (Iterator it, Iterator end, Filter f) :
    nt (it),
    iend (end),
    filt (f)
  {
    while (nt != iend && ! filt (*nt))
      ++nt;
  }

  /*! Access operations. */
  Iterator current_iterator() const
  {
    return (nt);
  }

  Iterator past_the_end () const
  {
    return (iend);
  }

  Filter filter () const
  {
    return (filt);
  }

  pointer ptr() const
  {
    return static_cast<pointer>(&(*nt));
  }

  /*! Equality operators. */
  bool operator== (const Self& it) const
  {
    return (nt == it.nt);
  }
  
  bool operator!= (const Self& it) const 
  {
    return !(*this == it);
  }
  
  bool operator< (const Self& it) const
  {
    return &(**this) < (&*it);
  }

  /*! Dereferencing operators. */
  reference operator*() const
  {
    return (*(ptr()));
  }
  
  pointer operator->() const
  {
    return ptr();
  }

  /*! Increment operators. */
  Self& operator++ ()
  {
    do
    {
      ++nt;
    } while (!(nt == iend) && ! filt (*nt));

    return (*this);
  }
   
  Self operator++ (int)
  {
    Self tmp = *this;
    ++(*this);
    return tmp;
  }

  /*! Decrement operators. */
  Self& operator-- ()
  {
    do
    {
      --nt;
    } while (!(nt == iend) && ! filt (*nt));

    return (*this);
  }

  Self operator-- (int)
  {
    Self tmp = *this;
    --(*this);
    return tmp;
  }
};

/*!
 * \class
 * An iterator adaptor for the filtering a DCEL const iterator (given as
 * CIterator_) using a given filter functor (Filter_).
 */
template <class CIterator_, class Filter_, class MIterator_,
          class Value_, class Diff_, class Category_>
class I_Filtered_const_iterator
{
public:

  typedef CIterator_                             Iterator;
  typedef Filter_                                Filter;
  typedef I_Filtered_const_iterator<CIterator_,
                                    Filter_,
                                    MIterator_,
                                    Value_,
                                    Diff_,
                                    Category_>   Self;

  typedef Category_                              iterator_category;
  typedef Value_                                 value_type;
  typedef const value_type&                      reference;
  typedef const value_type*                      pointer;
  typedef Diff_                                  difference_type;

  typedef I_Filtered_iterator<MIterator_, Filter_,
                              Value_, Diff_,
                              Category_>         mutable_iterator;

protected:

  Iterator       nt;       // The internal iterator (this member should not
                           // be renamed in order to comply with the
                           // HalfedgeDS circulators that refer to it).
  Iterator       iend;     // A past-the-end iterator.        
  Filter         filt;     // The filter functor.

public:

  /*! Constructors. */
  I_Filtered_const_iterator()
  {}

  I_Filtered_const_iterator (Iterator it) :
    nt (it),
    iend (it)
  {}

  I_Filtered_const_iterator (Iterator it, Iterator end) :
    nt (it),
    iend (end)
  {
    while (nt != iend && ! filt (*nt))
      ++nt;
  }

  I_Filtered_const_iterator (Iterator it, Iterator end, Filter f) :
    nt (it),
    iend (end),
    filt (f)
  {
    while (nt != iend && ! filt (*nt))
      ++nt;
  }

  I_Filtered_const_iterator (mutable_iterator it) :
    nt (it.current_iterator()),
    iend (it.past_the_end()),
    filt (it.filter())
  {
    //    while (nt != iend && ! filt (*nt))
    //      ++nt;
  }

  /*! Access operations. */
  Iterator current_iterator() const
  {
    return (nt);
  }

  Iterator past_the_end () const
  {
    return (iend);
  }

  Filter filter () const
  {
    return (filt);
  }

  pointer ptr() const
  {
    return static_cast<pointer>(&(*nt));
  }

  /*! Equality operators. */
  bool operator== (const Self& it) const
  {
    return (nt == it.nt);
  }
  
  bool operator!= (const Self& it) const 
  {
    return !(*this == it);
  }
  
  bool operator< (const Self& it) const
  {
    return &(**this) < (&*it);
  }

  /*! Dereferencing operators. */
  reference operator*() const
  {
    return (*(ptr()));
  }
  
  pointer operator->() const
  {
    return ptr();
  }

  /*! Increment operators. */
  Self& operator++ ()
  {
    do
    {
      ++nt;
    } while (!(nt == iend) && ! filt (*nt));

    return (*this);
  }
   
  Self operator++ (int)
  {
    Self tmp = *this;
    ++(*this);
    return tmp;
  }

  /*! Decrement operators. */
  Self& operator-- ()
  {
    do
    {
      --nt;
    } while (!(nt == iend) && ! filt (*nt));

    return (*this);
  }

  Self operator-- (int)
  {
    Self tmp = *this;
    --(*this);
    return tmp;
  }
};

} //namespace CGAL

namespace std {

#if defined(BOOST_MSVC)
#  pragma warning(push)
#  pragma warning(disable:4099) // For VC10 it is class hash 
#endif

#ifndef CGAL_CFG_NO_STD_HASH
  
template <class CIterator_, class Filter_, class MIterator_,
          class Value_, class Diff_, class Category_>
struct hash<CGAL::I_Filtered_const_iterator<CIterator_,
                                            Filter_,
                                            MIterator_,
                                            Value_,
                                            Diff_,
                                            Category_> > {

  typedef CGAL::I_Filtered_const_iterator<CIterator_,
                                          Filter_,
                                          MIterator_,
                                          Value_,
                                          Diff_,
                                          Category_>  I;

    std::size_t operator()(const I& i) const
    {
      return reinterpret_cast<std::size_t>(&*i) / sizeof(Value_);
    }
  };

  template <class Iterator_, class Filter_,
          class Value_, class Diff_, class Category_>
  struct hash<CGAL::I_Filtered_iterator<Iterator_,
                                        Filter_,
                                        Value_,
                                        Diff_,
                                        Category_> > {
  typedef CGAL::I_Filtered_iterator<Iterator_,
                                    Filter_,
                                    Value_,
                                    Diff_,
                                    Category_>  I;

    std::size_t operator()(const I& i) const
    {
      return reinterpret_cast<std::size_t>(&*i) / sizeof(typename I::value_type);
    }
  };

#endif // CGAL_CFG_NO_STD_HASH

#if defined(BOOST_MSVC)
#  pragma warning(pop)
#endif

} // namespace std

namespace  boost {
  template <typename T> struct hash;

template <class CIterator_, class Filter_, class MIterator_,
          class Value_, class Diff_, class Category_>
struct hash<CGAL::I_Filtered_const_iterator<CIterator_,
                                            Filter_,
                                            MIterator_,
                                            Value_,
                                            Diff_,
                                            Category_> > {

  typedef CGAL::I_Filtered_const_iterator<CIterator_,
                                          Filter_,
                                          MIterator_,
                                          Value_,
                                          Diff_,
                                          Category_>  I;

    std::size_t operator()(const I& i) const
    {
      return reinterpret_cast<std::size_t>(&*i) / sizeof(Value_);
    }
  };

  template <class Iterator_, class Filter_,
          class Value_, class Diff_, class Category_>
  struct hash<CGAL::I_Filtered_iterator<Iterator_,
                                        Filter_,
                                        Value_,
                                        Diff_,
                                        Category_> > {
  typedef CGAL::I_Filtered_iterator<Iterator_,
                                    Filter_,
                                    Value_,
                                    Diff_,
                                    Category_>  I;

    std::size_t operator()(const I& i) const
    {
      return reinterpret_cast<std::size_t>(&*i) / sizeof(typename I::value_type);
    }
  };

} // namespace boost
#endif
