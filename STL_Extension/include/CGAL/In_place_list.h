// Copyright (c) 2003
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>
//                 Lutz Kettner <kettner@mpi-sb.mpg.de>
//                 Sylvain Pion

#ifndef CGAL_IN_PLACE_LIST_H
#define CGAL_IN_PLACE_LIST_H 1

#include <CGAL/disable_warnings.h>

#include <CGAL/basic.h>
#include <cstddef>
#include <iterator>
#include <functional>
#include <algorithm>
#include <CGAL/memory.h>
#include <boost/functional/hash.hpp>

namespace CGAL {

// Forward declarations
namespace internal {
  template <class T, class Alloc> class In_place_list_iterator;
  template <class T, class Alloc> class In_place_list_const_iterator;
}

template <class T, bool managed, class Alloc = CGAL_ALLOCATOR(T)>
class In_place_list;

template < class T >
class In_place_sl_list_base {
public:
  T* next_link;        // forward pointer
};

template < class T >
class In_place_list_base {
public:
  In_place_list_base()
    : next_link(nullptr), prev_link(nullptr)
  {}

  T* next_link;        // forward pointer
  T* prev_link;        // backwards pointer
  //friend  class internal::In_place_list_iterator<T, Alloc>;
  //friend  class internal::In_place_list_const_iterator<T, Alloc>;
  //friend  class In_place_list<T,false, Alloc>;
  //friend  class In_place_list<T,true, Alloc>;
};


namespace internal {
  template <class T, class Alloc>
  class In_place_list_iterator {
  protected:
    T* node;
  public:
    friend  class In_place_list<T,false, Alloc>;
    friend  class In_place_list<T,true, Alloc>;

    typedef In_place_list_iterator<T, Alloc>  Self;
    typedef In_place_list_base<T>      Base;

    typedef T               value_type;
    typedef T*              pointer;
    typedef T&              reference;
    typedef std::size_t     size_type;
    typedef std::ptrdiff_t  difference_type;
    typedef std::bidirectional_iterator_tag   iterator_category;

    In_place_list_iterator() : node(0) {}
    In_place_list_iterator(T* x) : node(x) {}

    bool  operator==( const Self& x) const { return node == x.node; }
    bool  operator!=( const Self& x) const { return node != x.node; }
    bool  operator< ( const Self& x) const { return node< x.node;   }
    bool  operator<=( const Self& x) const { return node<= x.node;  }
    bool  operator> ( const Self& x) const { return node> x.node;   }
    bool  operator>=( const Self& x) const { return node>= x.node;  }
    T&    operator*()  const { return *node; }
    T*    operator->() const { return  node; }
    Self& operator++() {
      node = ((Base*)(node))->next_link;
      return *this;
    }
    Self  operator++(int) {
      Self tmp = *this;
      ++*this;
      return tmp;
    }
    Self& operator--() {
      node = ((Base*)(node))->prev_link;
      return *this;
    }
    Self  operator--(int) {
      Self tmp = *this;
      --*this;
      return tmp;
    }
  };
}

namespace internal {
  template <class T, class Alloc>
  class In_place_list_const_iterator {
  protected:
    const T* node;  // It's not Ptr. Otherwise traversal won't work.
  public:
    friend  class In_place_list<T,false, Alloc>;
    friend  class In_place_list<T,true, Alloc>;

    typedef In_place_list_const_iterator<T, Alloc> Self;
    typedef In_place_list_iterator<T, Alloc>       Iterator;
    typedef In_place_list_base<T>           Base;

    typedef T               value_type;
    typedef const T*        pointer;
    typedef const T&        reference;
    typedef std::size_t     size_type;
    typedef std::ptrdiff_t  difference_type;
    typedef std::bidirectional_iterator_tag   iterator_category;

    In_place_list_const_iterator() : node(0) {}
    In_place_list_const_iterator(Iterator i) : node(i.operator->()) {}
    In_place_list_const_iterator(const T* x) : node(x) {}

    bool     operator==( const Self& x) const { return node == x.node; }
    bool     operator!=( const Self& x) const { return node != x.node; }
    bool     operator< ( const Self& x) const { return node< x.node;   }
    bool     operator<=( const Self& x) const { return node<= x.node;  }
    bool     operator> ( const Self& x) const { return node> x.node;   }
    bool     operator>=( const Self& x) const { return node>= x.node;  }
    const T& operator*()  const { return *node; }
    const T* operator->() const { return  node; }
    Self& operator++() {
      node = ((const Base*)(node))->next_link;
      return *this;
    }
    Self  operator++(int) {
      Self tmp = *this;
      ++*this;
      return tmp;
    }
    Self& operator--() {
      node = ((const Base*)(node))->prev_link;
      return *this;
    }
    Self  operator--(int) {
      Self tmp = *this;
      --*this;
      return tmp;
    }
    In_place_list_iterator<T,Alloc>
    remove_const() const
    {
      return In_place_list_iterator<T,Alloc>(const_cast<T*>(node));
    }
  };



template <class T, class Alloc>
  std::size_t hash_value(const In_place_list_iterator<T,Alloc>&  i)
  {
    T* ptr = &*i;
    return reinterpret_cast<std::size_t>(ptr)/ sizeof(T);
  }


template <class T, class Alloc>
  std::size_t hash_value(const In_place_list_const_iterator<T,Alloc>&  i)
  {
    const T* ptr = &*i;
    return reinterpret_cast<std::size_t>(ptr)/ sizeof(T);
   }

}


template <class T, bool managed, class Alloc>
class In_place_list {

  // Bidirectional List Managing Objects in Place
  // --------------------------------------------
  //
  // DEFINITION An object of the class In_place_list<T,bool> is a
  // sequence that supports bidirectional iterators and allows constant time
  // insert and erase operations anywhere within the sequence. The
  // functionality is similar to the `list<T>' in the STL.
  //
  // The In_place_list<T,bool> manages element items in place. Two
  // pointers `T*' are expected in the class. For example the base class
  // `In_place_list_base<T>' can be used.
  //
  // The In_place_list<T,bool> does not copy element items during
  // insertion (unless otherwise stated for a function). On removal or
  // destruction of the list the element items are not deleted by default.
  // The second template parameter `bool' has to be set to `false' in this
  // case. If the In_place_list<T,bool> should take the responsibility
  // for the stored objects the `bool' parameter could be set to `true', in
  // which case the list will delete removed items and will delete all
  // remaining items on destruction. In any case, the `destroy()' member
  // function deletes all elements.
  //
  // On purpose, these two possible versions of In_place_list<T,bool>
  // are not assignment compatible to avoid confusions between the different
  // storage responsibilities.
  //
  // PARAMETERS
  //
  // The full classname is `In_place_list<T,bool managed = false, Alloc
  // = CGAL_ALLOCATOR(T)>'.
  //
  // TYPES

public:
  typedef Alloc           Allocator;
  typedef Alloc           allocator_type; // STL compliant

  // Note: the standard requires the following types to be equivalent
  // to T, T*, const T*, T&, const T&, size_t, and ptrdiff_t, respectively.
  // So we don't pass these types to the iterators explicitly.

  typedef typename std::allocator_traits<Allocator>::value_type            value_type;
  typedef typename std::allocator_traits<Allocator>::pointer               pointer;
  typedef typename std::allocator_traits<Allocator>::const_pointer         const_pointer;
  typedef typename std::allocator_traits<Allocator>::size_type             size_type;
  typedef typename std::allocator_traits<Allocator>::difference_type       difference_type;

  typedef value_type&       reference;
  typedef const value_type& const_reference;

  typedef internal::In_place_list_iterator<T, Alloc> iterator;
  typedef internal::In_place_list_const_iterator<T, Alloc> const_iterator;

  typedef std::reverse_iterator<iterator>         reverse_iterator;
  typedef std::reverse_iterator<const_iterator>   const_reverse_iterator;

  typedef In_place_list<T,managed,Alloc>          Self;

protected:
  Allocator allocator;

  pointer      node;
  size_type    length;

  // These are the only places where the allocator gets called.
  pointer get_node() {
    pointer p = allocator.allocate(1);
#ifdef CGAL_USE_ALLOCATOR_CONSTRUCT_DESTROY
    allocator.construct(p, value_type());
#else
    new (p) value_type;
#endif
    return p;
  }
  pointer get_node( const T& t) {
    pointer p = allocator.allocate(1);
#ifdef CGAL_USE_ALLOCATOR_CONSTRUCT_DESTROY
    std::allocator_traits<Allocator>::construct(allocator, p, t);
#else
    new (p) value_type(t);
#endif
    return p;
  }
  void put_node( pointer p) {
#ifdef CGAL_USE_ALLOCATOR_CONSTRUCT_DESTROY
    std::allocator_traits<Allocator>::destroy(allocator, p);
#else // not CGAL_USE_ALLOCATOR_CONSTRUCT_DESTROY
   p->~value_type();
#endif
    allocator.deallocate( p, 1);
  }

public:
  // CREATION
  //
  // New creation variable is: `l'

  explicit In_place_list() : length(0) {
    // introduces an empty list.
    node = get_node();
    (*node).next_link = node;
    (*node).prev_link = node;
  }
  void swap(Self& x) {
    std::swap(node, x.node);
    std::swap(length, x.length);
  }

  // ACCESS MEMBER FUNCTIONS

  allocator_type get_allocator() const { return allocator; }

  iterator       begin() { return (*node).next_link; }
  const_iterator begin() const { return (*node).next_link; }
  iterator       end() { return node; }
  const_iterator end() const { return node; }

  reverse_iterator       rbegin() { return reverse_iterator(end()); }
  const_reverse_iterator rbegin() const {
    return const_reverse_iterator(end());
  }
  reverse_iterator       rend() { return reverse_iterator(begin()); }
  const_reverse_iterator rend() const {
    return const_reverse_iterator(begin());
  }

  bool            empty() const    { return length == 0; }
  size_type       size() const     { return length; }
  size_type       max_size() const { return size_type(-1); }

  reference       front()          { return *begin(); }
  const_reference front() const    { return *begin(); }
  reference       back()           { return *(--end()); }
  const_reference back() const     { return *(--end()); }

  // INSERTION

  iterator insert(iterator position, T& x) {
    // inserts `t' in front of iterator `pos'. The return value points
    // to the inserted item.
    x.next_link = position.node;
    x.prev_link = (*position.node).prev_link;
    (*((*position.node).prev_link)).next_link = &x;
    (*position.node).prev_link = &x;
    ++length;
    return &x;
  }
  iterator insert(T* pos, T& x) {
    return insert( iterator(pos), x);
  }
  void push_front(T& x) { insert(begin(), x); }
  // inserts an item in front of list `l'.

  void push_back(T& x)  { insert(end(), x); }
  // inserts an item at the back of list `l'.

  void insert(iterator position, size_type n);
  // inserts n copies of `T()' in front of iterator `pos'.

  void insert(iterator position, size_type n, const T& x);
  // inserts n copies of `t' in front of iterator `pos'.

  void insert( T* pos, size_type n) { insert( iterator(pos), n); }
  void insert( T* pos, size_type n, const T& x) {
    insert( iterator(pos), n, x);
  }

  template <class InputIterator>
  void insert(iterator pos, InputIterator first, InputIterator last) {
    // inserts the range [`first, last') in front of iterator `pos'.
    while (first != last)
      insert(pos, *get_node(*first++));
  }

  template <class InputIterator>
  void insert(T* pos, InputIterator first, InputIterator last) {
    // inserts the range [`first, last') in front of iterator `pos'.
    while (first != last)
      insert(pos, *get_node(*first++));
  }

  void insert(T* pos, const T* first, const T* last) {
    insert( iterator(pos), const_iterator(first),
            const_iterator(last));
  }


  // REMOVAL

  void erase(iterator i) {
    // removes the item from list `l', where `pos' refers to.
    CGAL_assertion( length > 0);
    (*((*i.node).prev_link)).next_link = (*i.node).next_link;
    (*((*i.node).next_link)).prev_link = (*i.node).prev_link;
    if (managed)
      put_node(i.node);
    --length;
  }
  void erase(T* pos)  { erase( iterator( pos)); }

  void pop_front() { erase(begin()); }
  // removes the first item from list `l'.

  void pop_back() {
    // removes the last item from list `l'.
    iterator tmp = end();
    erase(--tmp);
  }

  void erase(iterator first, iterator last);
  // removes the items in the range [`first, last') from list `l'.

  void erase(T* first, T* last) {
    erase( iterator(first), iterator(last));
  }

  void clear() { erase( begin(), end()); }

  // CREATION (Continued)

  explicit In_place_list(size_type n, const T& value = T()) : length(0) {
    // introduces a list with n items, all initialized with copies of
    // value.
    node = get_node();
    (*node).next_link = node;
    (*node).prev_link = node;
    insert(begin(), n, value);
  }

  template <class InputIterator>
  In_place_list( InputIterator first, InputIterator last) : length(0) {
    // a list with copies from the range [`first,last').
    node = get_node();
    (*node).next_link = node;
    (*node).prev_link = node;
    insert( begin(), first, last);
  }

  In_place_list(const T* first, const T* last) : length(0) {
    // a list with copies from the range [`first,last').
    node = get_node();
    (*node).next_link = node;
    (*node).prev_link = node;
    insert(begin(), first, last);
  }
  In_place_list(const Self& x) : length(0) {
    // copy constructor. Each item in `l1' is copied.
    node = get_node();
    (*node).next_link = node;
    (*node).prev_link = node;
    insert(begin(), x.begin(), x.end());
  }
  ~In_place_list() {
    erase(begin(), end());
    put_node(node);
  }

  Self& operator=(const Self& x);

  void destroy();

  template <class InputIterator>
  void assign( InputIterator first, InputIterator last) {
    erase( begin(), end());
    insert( begin(), first, last);
  }

  void assign( size_type n, const T& t) {
    erase( begin(), end());
    insert( begin(), n, t);
  }

  void resize( size_type sz, const T& c = T()) {
    if ( sz > size())
      insert( end(), sz - size(), c);
    else if ( sz < size()) {
      iterator i = begin();
      while ( sz-- > 0)
        ++i;
      erase( i, end());
    }  // else do nothing
  }

  // COMPARISON OPERATIONS

  bool operator==( const Self& y) const {
    return size() == y.size() && std::equal(begin(), end(), y.begin());
  }

  bool operator!=( const Self& y) const {
    return size() != y.size() || ! std::equal(begin(),end(),y.begin());
  }

  bool operator<(const Self& y) const {
    return std::lexicographical_compare( begin(),end(),
                                         y.begin(),y.end());
  }
  bool operator> ( const Self& i) const { return i < *this; }
  bool operator<=( const Self& i) const { return !(i < *this); }
  bool operator>=( const Self& i) const { return !(*this < i); }

  // SPECIAL LIST OPERATIONS

protected:
  void transfer(iterator position, iterator first, iterator last) {
    // move the range [`first, last') before the position.
    (*((*last.node).prev_link)).next_link = position.node;
    (*((*first.node).prev_link)).next_link = last.node;
    (*((*position.node).prev_link)).next_link = first.node;
    T* tmp = (*position.node).prev_link;
    (*position.node).prev_link = (*last.node).prev_link;
    (*last.node).prev_link = (*first.node).prev_link;
    (*first.node).prev_link = tmp;
  }

public:
  void splice(iterator position, Self& x) {
    // inserts the list x before position `pos' and x becomes empty.
    // It takes constant time. Precondition: `&l != &x'.
    if (!x.empty()) {
      transfer(position, x.begin(), x.end());
      length += x.length;
      x.length = 0;
    }
  }
  void splice(T* position, Self& x) {
    splice( iterator(position), x);
  }
  void splice( iterator position, Self& x, iterator i) {
    // inserts an element pointed to by i from list x before position
    // `pos' and removes the element from x. It takes constant time. i
    // is a valid dereferenceable iterator of x. The result is
    // unchanged if `pos == i' or `pos == ++i'.
    iterator j = i;
    if (position == i || position == ++j) return;
    transfer(position, i, j);
    ++length;
    --x.length;
  }
  void splice(T* position, Self& x, T* i) {
    splice( iterator(position), x, iterator(i));
  }
  void splice(iterator pos, Self& x, iterator first, iterator last) {
    // inserts elements in the range [`first, last') before position
    // `pos' and removes the elements from x. It takes constant time
    // if `&x == $l'; otherwise, it takes linear time. [`first,
    // last') is a valid range in x. Precondition: `pos' is not in the
    // range [`first, last').
    if (first != last) {
      if (&x != this) {
        difference_type n = std::distance(first, last);
        x.length -= n;
        length += n;
      }
      transfer(pos, first, last);
    }
  }
  void splice(T* p, Self& x, T* first, T* last) {
    splice( iterator(p), x, iterator(first), iterator(last));
  }

  void remove(const T& value);
  // erases all elements e in the list l for which `e == value'.
  // It is stable. Precondition: a suitable `operator==' for the
  // type T.

  void reverse();
  // reverses the order of the elements in `l' in linear time.

  void unique();
  // erases all but the first element from every consecutive group
  // of equal elements in the list `l'. Precondition: a suitable
  // `operator==' for the type T.

  void merge(Self& x);
  // merges the list x into the list `l' and x becomes empty. It is
  // stable. Precondition: Both lists are increasingly sorted. A
  // suitable `operator<' for the type T.

  template < class StrictWeakOrdering >
  void merge(Self& x, StrictWeakOrdering ord)
  // merges the list x into the list `l' and x becomes empty.
  // It is stable.
  // Precondition: Both lists are increasingly sorted wrt. ord.
  {
    iterator first1 = begin();
    iterator last1 = end();
    iterator first2 = x.begin();
    iterator last2 = x.end();
    while (first1 != last1 && first2 != last2)
      if (ord(*first2, *first1)) {
        iterator next = first2;
        transfer(first1, first2, ++next);
        first2 = next;
      } else
        ++first1;
    if (first2 != last2)
      transfer(last1, first2, last2);
    length += x.length;
    x.length= 0;
  }

  void sort();
  // sorts the list `l' according to the `operator<' in time O(n
  // log n) where `n = size()'. It is stable. Precondition: a
  // suitable `operator<' for the type T.

  template < class StrictWeakOrdering >
  void sort(StrictWeakOrdering ord)
  // sorts the list `l' according to ord in time O(n log n)
  // where `n = size()'. It is stable.
  {
    if (size() < 2) return;
    In_place_list<T,managed,Alloc> carry;
    In_place_list<T,managed,Alloc> counter[64];
    int fill = 0;
    while (!empty()) {
      carry.splice(carry.begin(), *this, begin());
      int i = 0;
      while(i < fill && !counter[i].empty()) {
        counter[i].merge(carry, ord);
        carry.swap(counter[i++]);
      }
      carry.swap(counter[i]);
      if (i == fill)
        ++fill;
    }
    for (int i = 1; i < fill; ++i)
      counter[i].merge(counter[i-1], ord);
    swap(counter[fill-1]);
  }

};

template <class T, bool managed, class Alloc>
void In_place_list<T,managed,Alloc>::
insert(internal::In_place_list_iterator<T, Alloc> position, size_type n) {
  while (n--)
    insert(position, *get_node());
}

template <class T, bool managed, class Alloc>
void In_place_list<T,managed,Alloc>::
insert(internal::In_place_list_iterator<T, Alloc> position, size_type n, const T& x) {
  while (n--)
    insert(position, *get_node(x));
}

template <class T, bool managed, class Alloc>
void In_place_list<T,managed,Alloc>::
erase(internal::In_place_list_iterator<T, Alloc> first,
      internal::In_place_list_iterator<T, Alloc> last)
{
  while (first != last)
    erase(first++);
}

template <class T, bool managed, class Alloc>
In_place_list<T,managed,Alloc>&
In_place_list<T,managed,Alloc>::
operator=(const In_place_list<T,managed,Alloc>& x) {
  if (this != &x) {
    iterator first1 = begin();
    iterator last1  = end();
    const_iterator first2 = x.begin();
    const_iterator last2  = x.end();
    while (first1 != last1 && first2 != last2) {
      // Save the pointer values before assignment.
      // Assignment avoids unneccassary delete's and new's.
      T* tmp1 = (*first1).next_link;
      T* tmp2 = (*first1).prev_link;
      *first1 = *first2++;
      (*first1).next_link = tmp1;
      (*first1).prev_link = tmp2;
      ++first1;
    }
    if (first2 == last2)
      erase(first1, last1);
    else
      insert(last1, first2, last2);
  }
  return *this;
}

template <class T, bool managed, class Alloc>
void In_place_list<T,managed,Alloc>::
destroy() {
  iterator first = begin();
  iterator last  = end();
  while( first != last) {
    iterator i = first++;
    put_node(i.node);
  }
  length = 0;
  (*node).next_link = node;
  (*node).prev_link = node;
}

template <class T, bool managed, class Alloc>
void In_place_list<T,managed,Alloc>::remove(const T& value) {
  iterator first = begin();
  iterator last = end();
  while (first != last) {
    iterator next = first;
    ++next;
    if (*first == value)
      erase(first);
    first = next;
  }
}

template <class T, bool managed, class Alloc>
void In_place_list<T,managed,Alloc>::reverse() {
  if (size() < 2) return;
  for (iterator first = ++begin(); first != end();) {
    iterator old = first++;
    transfer(begin(), old, first);
  }
}

template <class T, bool managed, class Alloc>
void In_place_list<T,managed,Alloc>::unique() {
  iterator first = begin();
  iterator last = end();
  if (first == last) return;
  iterator next = first;
  while (++next != last) {
    if (*first == *next)
      erase(next);
    else
      first = next;
    next = first;
  }
}

template <class T, bool managed, class Alloc>
void In_place_list<T,managed,Alloc>::merge(In_place_list<T,managed,Alloc>& x) {
  iterator first1 = begin();
  iterator last1 = end();
  iterator first2 = x.begin();
  iterator last2 = x.end();
  while (first1 != last1 && first2 != last2)
    if (*first2 < *first1) {
      iterator next = first2;
      transfer(first1, first2, ++next);
      first2 = next;
    } else
      ++first1;
  if (first2 != last2)
    transfer(last1, first2, last2);
  length += x.length;
  x.length= 0;
}

template <class T, bool managed, class Alloc>
void In_place_list<T,managed,Alloc>::sort() {
  if (size() < 2) return;
  In_place_list<T,managed,Alloc> carry;
  In_place_list<T,managed,Alloc> counter[64];
  int fill = 0;
  while (!empty()) {
    carry.splice(carry.begin(), *this, begin());
    int i = 0;
    while(i < fill && !counter[i].empty()) {
      counter[i].merge(carry);
      carry.swap(counter[i++]);
    }
    carry.swap(counter[i]);
    if (i == fill)
      ++fill;
  }
  for (int i = 1; i < fill; ++i)
    counter[i].merge(counter[i-1]);
  swap(counter[fill-1]);
}

} //namespace CGAL

namespace std {

#if defined(BOOST_MSVC)
#  pragma warning(push)
#  pragma warning(disable:4099) // For VC10 it is class hash
#endif

#ifndef CGAL_CFG_NO_STD_HASH

  template < class T, class Alloc >
  struct hash<CGAL::internal::In_place_list_iterator<T, Alloc> >
    : public CGAL::cpp98::unary_function<CGAL::internal::In_place_list_iterator<T, Alloc>, std::size_t>  {

    std::size_t operator()(const CGAL::internal::In_place_list_iterator<T, Alloc>& i) const
    {
      const T* ptr = &*i;
      return reinterpret_cast<std::size_t>(ptr)/ sizeof(T);
    }
  };

  template < class T, class Alloc >
  struct hash<CGAL::internal::In_place_list_const_iterator<T, Alloc> >
    : public CGAL::cpp98::unary_function<CGAL::internal::In_place_list_const_iterator<T, Alloc>, std::size_t> {

    std::size_t operator()(const CGAL::internal::In_place_list_const_iterator<T, Alloc>& i) const
    {
      const T* ptr = &*i;
      return reinterpret_cast<std::size_t>(ptr)/ sizeof(T);
    }
  };
#endif // CGAL_CFG_NO_STD_HASH

#if defined(BOOST_MSVC)
#  pragma warning(pop)
#endif

} // namespace std

#include <CGAL/enable_warnings.h>

#endif // CGAL_IN_PLACE_LIST_H
