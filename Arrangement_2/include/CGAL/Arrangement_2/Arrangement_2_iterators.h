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

#ifndef CGAL_ARRANGEMENT_2_ITERATORS_H
#define CGAL_ARRANGEMENT_2_ITERATORS_H

/*! \file
 * Definitions of auxiliary iterator adaptors.
 */

CGAL_BEGIN_NAMESPACE

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
    return (reinterpret_cast<value_type *> (*iter));
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
    return (reinterpret_cast<const value_type *> (*iter));
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

CGAL_END_NAMESPACE

#endif
