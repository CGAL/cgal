// Copyright (c) 2003,2004,2007-2010  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:  $
// $Id:  $
//
// Author(s)     : Clement Jamin

#ifdef LINKED_WITH_TBB

#ifndef CGAL_CONCURRENT_COMPACT_CONTAINER_H
#define CGAL_CONCURRENT_COMPACT_CONTAINER_H

//#define USE_BOOST_ITERATOR

// CJTODO: can we remove some of these includes?
#include <algorithm>
#include <vector>

// Boost
//#include <boost/iterator/filter_iterator.hpp>

// CGAL
#include <CGAL/Default.h>
#ifdef USE_BOOST_ITERATOR
  #include <boost/iterator/filter_iterator.hpp>
#else
  #include <CGAL/iterator.h>
#endif

// TBB
#include <tbb/concurrent_vector.h>
#include <tbb/enumerable_thread_specific.h>

namespace CGAL {

// The traits class describes the way to access the pointer.
// It can be specialized.
// We use the same pointer than Compact_container, as the need is the same
template < class T >
struct Concurrent_compact_container_traits {
  static void *   pointer(const T &t) { return t.for_compact_container(); }
  static void * & pointer(T &t)       { return t.for_compact_container(); }
};

// A basic "do nothing" CCC_strategy_base
// One can inheritate from it for partial specialisation
template <typename Element>
class CCC_strategy_base {
public:
  // Do nothing
  static unsigned int get_erase_counter(const Element &) { return 0; }
  static void set_erase_counter(Element &, unsigned int) {}
  static void increment_erase_counter(Element &) {}
};


// A CCC_strategy managing an internal counter
template <typename Element>
class CCC_strategy_with_counter
{
public:    
  static unsigned int get_erase_counter(const Element &e) 
  {
    return e.get_erase_counter(); 
  }

  static void set_erase_counter(Element &e, unsigned int c) 
  {
    e.set_erase_counter(c);
  }

  static void increment_erase_counter(Element &e) 
  {
    e.increment_erase_counter();
  }
};

// Class Concurrent_compact_container
//
// Strategy_ is a functor which provides several functions
// See documentation
//
template < class T, class Strategy_ = Default >
class Concurrent_compact_container
  : private tbb::concurrent_vector<T>
{
  typedef Strategy_                                                 Strat;
  typedef typename Default::Get<Strat, CCC_strategy_base<T> >::type Strategy;
  typedef Concurrent_compact_container <T, Strat>                   Self;
  typedef Concurrent_compact_container_traits <T>                   Traits;
  typedef tbb::concurrent_vector<T>                                 Base;
  typedef typename Base::iterator                                   Base_iterator;
  typedef typename Base::const_iterator                             Base_const_iterator;
  typedef typename Base::reverse_iterator                           Base_reverse_iterator;
  typedef typename Base::const_reverse_iterator                     Base_const_reverse_iterator;
  
#ifdef USE_BOOST_ITERATOR
  struct Is_element_valid
  {
    bool operator()(const T &e) const
    { 
      return get_type(e) == USED;
    }
  };
#else
  struct Is_element_invalid
  {
    template <typename Iterator>
    bool operator()(const Iterator &it) const
    { 
      return get_type(*it) != USED;
    }
  };
#endif

public:
  typedef T                                                 value_type;
  typedef typename Base::allocator_type                     allocator_type;
  typedef typename Base::reference                          reference;
  typedef typename Base::const_reference                    const_reference;
  typedef typename Base::pointer                            pointer;
  typedef typename Base::const_pointer                      const_pointer;
  typedef typename Base::size_type                          size_type;
  typedef typename Base::difference_type                    difference_type;
#ifdef USE_BOOST_ITERATOR
  typedef boost::filter_iterator<
    Is_element_valid, Base_iterator>                        iterator;
  typedef boost::filter_iterator<
    Is_element_valid, Base_const_iterator>                  const_iterator;
  typedef boost::filter_iterator<
    Is_element_valid, Base_reverse_iterator>                reverse_iterator;
  typedef boost::filter_iterator<
    Is_element_valid, Base_const_reverse_iterator>          const_reverse_iterator;
#else
  typedef Filter_iterator<
    Base_iterator, Is_element_invalid>                      iterator;
  typedef Filter_iterator<
    Base_const_iterator, Is_element_invalid>                const_iterator;
  typedef Filter_iterator<
    Base_reverse_iterator, Is_element_invalid>              reverse_iterator;
  typedef Filter_iterator<
    Base_const_reverse_iterator, Is_element_invalid>        const_reverse_iterator;
#endif
  // CJTODO TEMP
  /*typedef typename Base::iterator               iterator;
  typedef typename Base::const_iterator         const_iterator;
  typedef typename Base::reverse_iterator       reverse_iterator;
  typedef typename Base::const_reverse_iterator const_reverse_iterator;*/

  Concurrent_compact_container()
  {
    init (); 
  }

  template < class InputIterator >
  Concurrent_compact_container(InputIterator first, InputIterator last)
  {
    init();
    std::copy(first, last, CGAL::inserter(*this));
  }

  // The copy constructor and assignment operator preserve the iterator order
  Concurrent_compact_container(const Concurrent_compact_container &c)
  : Base(c)
  {
    init();
    std::copy(c.begin(), c.end(), CGAL::inserter(*this));
  }

  Concurrent_compact_container & operator=(
    const Concurrent_compact_container &c)
  {
    if (&c != this) {
      Self tmp(c);
      swap(tmp);
    }
    return *this;
  }

  ~Concurrent_compact_container()
  {
    clear();
  }

  void init()
  {
  }

  void swap(Self &c)
  {
    std::swap(m_free_lists, c.m_free_lists);
    Base::swap(c);
  }

  //============= Iterators =============

#ifdef USE_BOOST_ITERATOR
  // CJTODO: remplacer Is_element_invalid() par un membre de la classe ?
  iterator begin()
  {
    return iterator(Base::begin(), Base::end());
  }
  iterator end()
  {
    return iterator(Base::end(), Base::end());
  }
  const_iterator begin() const
  {
    return const_iterator(Base::begin(), Base::end()); 
  }
  const_iterator end() const
  {
    return const_iterator(Base::end(), Base::end()); 
  }
  reverse_iterator rbegin()
  {
    return reverse_iterator(Base::rbegin(), Base::rend());
  }
  reverse_iterator rend()
  {
    return reverse_iterator(Base::rend(), Base::rend());
  }
  const_reverse_iterator rbegin() const 
  { 
    return const_reverse_iterator(Base::rbegin(), Base::rend()); 
  }
  const_reverse_iterator rend() const 
  {
    return const_reverse_iterator(Base::rend(), Base::rend()); 
  }
#else
  // CJTODO: remplacer Is_element_invalid() par un membre de la classe ?
  iterator begin()
  {
    return iterator(Base::begin(), Is_element_invalid(), Base::end());
  }
  iterator end()
  {
    return iterator(Base::end(), Is_element_invalid(), Base::end());
  }
  const_iterator begin() const
  {
    return const_iterator(Base::begin(), Is_element_invalid(), Base::end()); 
  }
  const_iterator end() const
  {
    return const_iterator(Base::end(), Is_element_invalid(), Base::end()); 
  }
  reverse_iterator rbegin()
  {
    return reverse_iterator(Base::rbegin(), Is_element_invalid(), Base::rend());
  }
  reverse_iterator rend()
  {
    return reverse_iterator(Base::rend(), Is_element_invalid(), Base::rend());
  }
  const_reverse_iterator rbegin() const 
  { 
    return const_reverse_iterator(Base::rbegin(), Is_element_invalid(), Base::rend()); 
  }
  const_reverse_iterator rend() const 
  {
    return const_reverse_iterator(Base::rend(), Is_element_invalid(), Base::rend()); 
  }
#endif

  //============== Emplace, insert, erase... ==============

  
  // Special insert methods that construct the objects in place
  // (just forward the arguments to the constructor, to optimize a copy).
#ifndef CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES
  template < typename... Args >
  iterator
  emplace(const Args&... args)
  {
    // Free list empty?
    if (m_free_lists.local().empty())
      return insert(value_type(args...));

    iterator it_free_element = m_free_lists.local().back();
    m_free_lists.local().pop_back();
    // Call free_element constructor
    new (&*it_free_element) value_type(args...);
    CGAL_assertion(get_type(*it_free_element) == USED);
    return it_free_element;
  }
#else
  // inserts a default constructed item.
  iterator emplace()
  {
    // Free list empty?
    if (m_free_lists.local().empty())
      return insert(value_type());

    iterator it_free_element = m_free_lists.local().back();
    m_free_lists.local().pop_back();
    // Call free_element constructor
    new (&*it_free_element) value_type();
    CGAL_assertion(get_type(*it_free_element) == USED);
    return it_free_element;
  }
  
  template < typename T1 >
  iterator
  emplace(const T1 &t1)
  {
    // Free list empty?
    if (m_free_lists.local().empty())
      return insert(value_type(t1));

    iterator it_free_element = m_free_lists.local().back();
    m_free_lists.local().pop_back();
    // Call free_element constructor
    new (&*it_free_element) value_type(t1);
    CGAL_assertion(get_type(*it_free_element) == USED);
    return it_free_element;
  }

  template < typename T1, typename T2 >
  iterator
  emplace(const T1 &t1, const T2 &t2)
  {
    // Free list empty?
    if (m_free_lists.local().empty())
      return insert(value_type(t1, t2));

    iterator it_free_element = m_free_lists.local().back();
    m_free_lists.local().pop_back();
    // Call free_element constructor
    new (&*it_free_element) value_type(t1, t2);
    CGAL_assertion(get_type(*it_free_element) == USED);
    return it_free_element;
  }

  template < typename T1, typename T2, typename T3 >
  iterator
  emplace(const T1 &t1, const T2 &t2, const T3 &t3)
  {
    // Free list empty?
    if (m_free_lists.local().empty())
      return insert(value_type(t1, t2, t3));

    iterator it_free_element = m_free_lists.local().back();
    m_free_lists.local().pop_back();
    // Call free_element constructor
    new (&*it_free_element) value_type(t1, t2, t3);
    CGAL_assertion(get_type(*it_free_element) == USED);
    return it_free_element;
  }

  template < typename T1, typename T2, typename T3, typename T4 >
  iterator
  emplace(const T1 &t1, const T2 &t2, const T3 &t3, const T4 &t4)
  {
    // Free list empty?
    if (m_free_lists.local().empty())
      return insert(value_type(t1, t2, t3, t4));

    iterator it_free_element = m_free_lists.local().back();
    m_free_lists.local().pop_back();
    // Call free_element constructor
    new (&*it_free_element) value_type(t1, t2, t3, t4);
    CGAL_assertion(get_type(*it_free_element) == USED);
    return it_free_element;
  }

  template < typename T1, typename T2, typename T3, typename T4, typename T5 >
  iterator
  emplace(const T1 &t1, const T2 &t2, const T3 &t3, const T4 &t4,
	  const T5 &t5)
  {
    // Free list empty?
    if (m_free_lists.local().empty())
      return insert(value_type(t1, t2, t3, t4, t5));

    iterator it_free_element = m_free_lists.local().back();
    m_free_lists.local().pop_back();
    // Call free_element constructor
    new (&*it_free_element) value_type(t1, t2, t3, t4, t5);
    CGAL_assertion(get_type(*it_free_element) == USED);
    return it_free_element;
  }

  template < typename T1, typename T2, typename T3, typename T4,
             typename T5, typename T6 >
  iterator
  emplace(const T1 &t1, const T2 &t2, const T3 &t3, const T4 &t4,
          const T5 &t5, const T6 &t6)
  {
    // Free list empty?
    if (m_free_lists.local().empty())
      return insert(value_type(t1, t2, t3, t4, t5, t6));

    iterator it_free_element = m_free_lists.local().back();
    m_free_lists.local().pop_back();
    // Call free_element constructor
    new (&*it_free_element) value_type(t1, t2, t3, t4, t5, t6);
    CGAL_assertion(get_type(*it_free_element) == USED);
    return it_free_element;
  }

  template < typename T1, typename T2, typename T3, typename T4,
             typename T5, typename T6, typename T7 >
  iterator
  emplace(const T1 &t1, const T2 &t2, const T3 &t3, const T4 &t4,
          const T5 &t5, const T6 &t6, const T7 &t7)
  {
    // Free list empty?
    if (m_free_lists.local().empty())
      return insert(value_type(t1, t2, t3, t4, t5, t6, t7));

    iterator it_free_element = m_free_lists.local().back();
    m_free_lists.local().pop_back();
    // Call free_element constructor
    new (&*it_free_element) value_type(t1, t2, t3, t4, t5, t6, t7);
    CGAL_assertion(get_type(*it_free_element) == USED);
    return it_free_element;
  }

  template < typename T1, typename T2, typename T3, typename T4,
             typename T5, typename T6, typename T7, typename T8 >
  iterator
  emplace(const T1 &t1, const T2 &t2, const T3 &t3, const T4 &t4,
          const T5 &t5, const T6 &t6, const T7 &t7, const T8 &t8)
  {
    // Free list empty?
    if (m_free_lists.local().empty())
      return insert(value_type(t1, t2, t3, t4, t5, t6, t7, t8));

    iterator it_free_element = m_free_lists.local().back();
    m_free_lists.local().pop_back();
    // Call free_element constructor
    new (&*it_free_element) value_type(t1, t2, t3, t4, t5, t6, t7, t8);
    CGAL_assertion(get_type(*it_free_element) == USED);
    return it_free_element;
  }
#endif // CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES


  iterator insert(const T &t)
  {
#ifdef USE_BOOST_ITERATOR
    return Base::push_back(t);
#else
    return iterator(Base::push_back(t), Is_element_invalid(), Base::end());
#endif
  }

  template < class InputIterator >
  void insert(InputIterator first, InputIterator last)
  {
    for (; first != last; ++first)
      insert(*first);
  }
  
  void erase(iterator x)
  {
    CGAL_precondition(get_type(*x) == USED);
    Strategy::increment_erase_counter(*x);
    get_allocator().destroy(&*x);
/*#ifndef CGAL_NO_ASSERTIONS
    std::memset(&*x, 0, sizeof(T));
#endif*/
    put_on_free_list(x);
  }

  void erase(iterator first, iterator last) {
    while (first != last)
      erase(first++);
  }
  
  // Warning: not concurrent with respect to insertion/removal
  size_type size() const
  {
    size_t size = Base::size();
    for( TLS_Free_lists::iterator it_list = m_free_lists.begin() ; 
         it_list != m_free_lists.end() ; 
         ++it_list )
    {
      size -= it_list->size();
    }
    return size;
  }
  
  // Returns whether the iterator "cit" is in the range [begin(), end()].
  // This function is mostly useful for purposes of efficient debugging at
  // higher levels.
  bool owns(const_iterator cit) const
  {
    if (cit == end())
      return true;

    if (get_type(*cit) == FREE)
      return false;

    const_iterator it = begin();
    const_iterator it_end = end();
    for( ; it != it_end ; ++it)
    {
      if (it == cit)
        return true;
    }

    return false;
  }

  bool owns_dereferencable(const_iterator cit) const
  {
    return cit != end() && owns(cit);
  }

  // Warning: not concurrent with respect to insertion/removal
  void clear()
  {
    Base::clear();
    for( TLS_Free_lists::iterator it_list = m_free_lists.begin() ; 
         it_list != m_free_lists.end() ; 
         ++it_list )
    {
      it_list->clear();
    }
  }

protected:
  
  void put_on_free_list(const iterator &it)
  {
    set_type(*it, FREE);
    m_free_lists.local().push_back(it);
  }

  // [Inspired from Compact_container
  // Definition of the bit squatting :
  // =================================
  // ptr is composed of a pointer part and the last 1 bit.
  // Here is the meaning of each of the 8 cases.
  //
  //                  value of the last 1 bit as "Type"
  // pointer part     0            1
  //         NULL     user elt     free_list end
  //      != NULL     user elt     free elt
  //
  // meaning of ptr : user stuff   ptr to the iterator of the next free element

  enum Type { USED = 0, FREE = 1 };

  // The bit squatting is implemented by casting pointers to (char *), then
  // subtracting to NULL, doing bit manipulations on the resulting integer,
  // and converting back.

  static char * clean_pointer(char * p)
  {
    return ((p - (char *) NULL) & ~ (std::ptrdiff_t) FREE) + (char *) NULL;
  }
    
  /*
  // Returns the pointee, cleaned up from the squatted bits.
  static iterator *clean_pointee(const iterator *p_it)
  {
    return (iterator *) clean_pointer((char *) Traits::pointer(**p_it));
  }
  */

  // Get the type of the pointee.
  static Type get_type(const_reference e)
  {
    char * p = (char *) Traits::pointer(e);
    return (Type) (p - clean_pointer(p));
  }

  // Sets the type of the pointee.
  static void set_type(reference elt, Type t)
  {
    CGAL_precondition(0 <= t && t < 2);
    Traits::pointer(elt) = (void *) ((char *) Traits::pointer(elt) + (int) t);
  }

  // ============= Protected member variables =============
  typedef std::vector<iterator>                       Free_list;
  typedef tbb::enumerable_thread_specific<Free_list>  TLS_Free_lists;

  TLS_Free_lists        m_free_lists;
}; // class Concurrent_compact_container

//namespace internal {
//
//  template < class Container, bool Const >
//  class CCC_iterator
//   : public Filter_iterator<Container::iterator, CCC_Validity_tester>
//  {
//    typedef Filter_iterator<Container::iterator, CCC_Validity_tester> Base;
//    typedef CCC_iterator<Container, Const>                            Self;
//  public:
//
//    CCC_iterator() : Base() {}
//    CCC_iterator(const Base &b) : Base(b) {}
//
//    Self & operator++() { Base::operator++(); return *this; }
//    Self & operator--() { Base::operator--(); return *this; }
//    Self operator++(int) { Self tmp(*this); ++(*this); return tmp; }
//    Self operator--(int) { Self tmp(*this); --(*this); return tmp; }
//  };
//} // namespace internal

} //namespace CGAL

#endif // CGAL_CONCURRENT_COMPACT_CONTAINER_H
#endif // TBB_LINKED
