// Copyright (c) 2022 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
#ifndef CGAL_COMPACT_CONTAINER_WITH_INDEX_H
#define CGAL_COMPACT_CONTAINER_WITH_INDEX_H

#include <CGAL/Compact_container.h>
#include <CGAL/tags.h>
#include <limits>
#include <algorithm>
#include <stack>

// An STL like container, similar to Compact_container, but uses indices
// instead of handles.
// - free list can be stored in a std::deque or in a data member stored in T
// - Boolean used to mark used/free elements can be stored in a std::bitset
//   or in the most significant bit of the data member.

// TODO low priority :
// - rebind<> the allocator
// - Exception safety guarantees
// - Thread safety guarantees
// - std requirements on iterators says all defined operations are constant
//   time amortized (it's not true here, maybe it could be with some work...)
// - all this is expected especially when there are not so many free objects
//   compared to the allocated elements.
// - Currently, end() can be invalidated on erase() if the upper bound has changed.

namespace CGAL {

// default policy: each time the array of element is increased, the size
// is multiply by 2 (like for std::vector).
template<unsigned int k>
struct Multiply_by_two_policy_for_cc_with_size
{
  static const unsigned int first_block_size=k;

  template<typename Compact_container>
  static void increase_size(Compact_container& cc)
  { cc.next_delta=cc.capacity_; }
};

// constant size policy: the size of the array is always extended by the same
// number of elements.
template<unsigned int k>
struct Constant_size_policy_for_cc_with_size
{
  static const unsigned int first_block_size=k;

  template<typename Compact_container>
  static void increase_size(Compact_container& /*cc*/)
  {}
};

// The traits class describes the way to access the size_type.
// It can be specialized.
template < class T, class size_type >
struct Compact_container_with_index_traits {
  static size_type size_t(const T &t)
  { return t.for_compact_container(); }
  static void set_size_t(T &t, size_type v)
  { t.for_compact_container(v); }
};

namespace internal {
  template < class DSC, bool Const>
  class CC_iterator_with_index;

  template < class T, class ST >
  class MyIndex;
}

// Free list management: three versions:
// * free elements stored in a std::deque and used Booleans in a std::vector
// * free elements stored in a data member of T and used Booleans in a std::vector
// * free elements stored in a data member of T, and used Booleans in its most significant bit
// Note that version deque and data member does not exist (no interest?)
// (1) By default, use a deque and a vector.
template<typename CC_with_index, class Use_deque, class Use_vector>
class Free_list_management
{};

// (1) Case with a deque for the freelist and a vector for Booleans.
template<typename CC_with_index>
class Free_list_management<CC_with_index, CGAL::Tag_true, CGAL::Tag_true>
{
  using Self=Free_list_management<CC_with_index, CGAL::Tag_true, CGAL::Tag_true>;
  using T=typename CC_with_index::value_type;
  using size_type=typename CC_with_index::size_type;

public:
  static const size_type null_descriptor=(std::numeric_limits<size_type>::max)();

  Free_list_management(CC_with_index* cc_with_index):
    m_cc_with_index(cc_with_index)
  {}

  void init()
  {
    m_free_list=std::stack<size_type>();
    assert(m_free_list.empty());
    if(m_cc_with_index->capacity()>0)
    { m_used.assign(m_cc_with_index->capacity(), false); }
    else { m_used.clear(); }
    m_cc_with_index->upper_bound_=0;
  }

  void increase_to(size_type old_size)
  {
    CGAL_USE(old_size);
    CGAL_assertion(m_cc_with_index->capacity()>old_size);
    m_used.resize(m_cc_with_index->capacity(), false);
    // m_cc_with_index->upper_bound_ does not change.
  }

  void swap(Self& other)
  {
    // We should not swap m_cc_with_index, but only the content of the free list
    m_free_list.swap(other.m_free_list);
    m_used.swap(other.m_used);
  }

  bool is_empty() const
  { return m_free_list.empty() &&
        m_cc_with_index->upper_bound()>=m_cc_with_index->capacity(); }

  bool is_used(size_type i) const
  {
    CGAL_assertion(i<m_cc_with_index->capacity() && i!=null_descriptor);
    return m_used[i];
  }

  // Push the ith element on the free list: it becomes free
  void push(size_type i)
  {
    CGAL_assertion(i<m_cc_with_index->capacity() && i!=null_descriptor);
    CGAL_assertion(is_used(i));
    m_used[i]=false;
    if(i+1==m_cc_with_index->upper_bound())
    { --m_cc_with_index->upper_bound_; }
    else
    { m_free_list.push(i); }
  }

  size_type top() const
  {
    CGAL_assertion(!is_empty());
    if(m_cc_with_index->upper_bound()<m_cc_with_index->capacity())
    { return m_cc_with_index->upper_bound(); }
    return m_free_list.top();
  }

  // Pop one element from the free list (the top): it becomes used
  size_type pop()
  {
    CGAL_assertion(!is_empty());
    CGAL_assertion(!is_used(top()));
    size_type res=m_cc_with_index->upper_bound();
    if(m_cc_with_index->upper_bound()<m_cc_with_index->capacity())
    {
      res=m_cc_with_index->upper_bound();
      ++m_cc_with_index->upper_bound_;
    }
    else
    {
      res=m_free_list.top();
      m_free_list.pop();
    }
    CGAL_assertion(res<m_cc_with_index->upper_bound());
    CGAL_assertion(res!=null_descriptor);
    CGAL_assertion(!is_used(res));
    m_used[res]=true;
    return res;
  }

  void copy_special_data(const T& /*src*/, T& /*dest*/)
  {}

protected:
  CC_with_index* const  m_cc_with_index;
  std::stack<size_type> m_free_list;
  std::vector<bool>     m_used;
};

// (2) Case with the "in place" free list, and a vector for Booleans.
template<typename CC_with_index>
class Free_list_management<CC_with_index, CGAL::Tag_false, CGAL::Tag_true>
{
  using Self=Free_list_management<CC_with_index, CGAL::Tag_false, CGAL::Tag_true>;
  using T=typename CC_with_index::value_type;
  using size_type=typename CC_with_index::size_type;
  using Traits=Compact_container_with_index_traits <T, size_type>;

public:
  static const size_type null_descriptor=(std::numeric_limits<size_type>::max)();

  Free_list_management(CC_with_index* cc_with_index):
    m_cc_with_index(cc_with_index)
  {}

  void init()
  {
    if(m_cc_with_index->capacity()>0)
    {
      m_used.assign(m_cc_with_index->capacity(), false);
      m_free_list=0;
      for(size_type i=0;
          i<static_cast<size_type>(m_cc_with_index->capacity()-1); ++i)
      { Traits::set_size_t((*m_cc_with_index)[i], i+1); }
      // Next of the last element is null_descriptor.
      Traits::set_size_t((*m_cc_with_index)[m_cc_with_index->capacity()-1],
          null_descriptor);
    }
    else
    {
      m_free_list=null_descriptor;
      m_used.clear();
    }
    m_cc_with_index->upper_bound_=0;
  }

  void increase_to(size_type old_size)
  {
    CGAL_assertion(m_cc_with_index->capacity()>old_size);
    m_used.resize(m_cc_with_index->capacity(), false);
    for(size_type i=old_size;
        i<static_cast<size_type>(m_cc_with_index->capacity()-1); ++i)
    { Traits::set_size_t((*m_cc_with_index)[i], i+1); }
    // Next of the last element is null_descriptor.
    Traits::set_size_t((*m_cc_with_index)[m_cc_with_index->capacity()-1],
        m_free_list);
    m_free_list=old_size;
    // m_cc_with_index->upper_bound_ does not change.
  }

  void swap(Self& other)
  {
    // We should not swap m_cc_with_index, but only the content of the free list
    std::swap(m_free_list, other.m_free_list);
    m_used.swap(other.m_used);
  }

  bool is_empty() const
  { return m_free_list==null_descriptor; }

  bool is_used(size_type i) const
  {
    CGAL_assertion(i<m_cc_with_index->capacity() && i!=null_descriptor);
    return m_used[i];
  }

  // Push the ith element on the free list: it becomes free
  void push(size_type i)
  {
    CGAL_assertion(i<m_cc_with_index->capacity() && i!=null_descriptor);
    CGAL_assertion(is_used(i));
    m_used[i]=false;
    if(i+1==m_cc_with_index->upper_bound())
    { --m_cc_with_index->upper_bound_; }
    Traits::set_size_t((*m_cc_with_index)[i], m_free_list);
    m_free_list=i;
  }

  size_type top() const
  {
    CGAL_assertion(!is_empty());
    return m_free_list;
  }

  // Pop one element from the free list (the top): it becomes used
  size_type pop()
  {
    CGAL_assertion(!is_empty());
    CGAL_assertion(!is_used(top()));
    size_type res=m_free_list;
    if(res>=m_cc_with_index->upper_bound())
    { m_cc_with_index->upper_bound_=res+1; }
    CGAL_assertion(res<m_cc_with_index->upper_bound());
    CGAL_assertion(res!=null_descriptor);
    CGAL_assertion(!is_used(res));
    m_free_list=Traits::size_t((*m_cc_with_index)[res]);
    m_used[res]=true;
    return res;
  }

  void copy_special_data(const T& src, T& dest)
  { Traits::set_size_t(dest, Traits::size_t(src)); }

protected:
  CC_with_index* const m_cc_with_index;
  size_type            m_free_list; // First free element, null_descriptor if no free
  std::vector<bool>    m_used;
};

// (3) Case with the "in place" free list, and "in place" Booleans.
template<typename CC_with_index>
class Free_list_management<CC_with_index, CGAL::Tag_false, CGAL::Tag_false>
{
  using Self=Free_list_management<CC_with_index, CGAL::Tag_false, CGAL::Tag_false>;
  using T=typename CC_with_index::value_type;
  using size_type=typename CC_with_index::size_type;
  using Traits=Compact_container_with_index_traits <T, size_type>;

public:
  static const size_type null_descriptor=(std::numeric_limits<size_type>::max)()/2;

  Free_list_management(CC_with_index* cc_with_index):
    m_cc_with_index(cc_with_index)
  {}

  void init()
  {
    if(m_cc_with_index->capacity()>0)
    {
      m_free_list=0;
      for(size_type i=0;
          i<static_cast<size_type>(m_cc_with_index->capacity()-1); ++i)
      { static_set_val((*m_cc_with_index)[i], i+1, FREE); }
      // Next of the last element is null_descriptor.
      static_set_val((*m_cc_with_index)[m_cc_with_index->capacity()-1],
          null_descriptor, FREE);
    }
    else
    {
      m_free_list=null_descriptor;
    }
    m_cc_with_index->upper_bound_=0;
  }

  void increase_to(size_type old_size)
  {
    CGAL_assertion(m_cc_with_index->capacity()>old_size);
    for(size_type i=old_size;
        i<static_cast<size_type>(m_cc_with_index->capacity()-1); ++i)
    { static_set_val((*m_cc_with_index)[i], i+1, FREE); }
    static_set_val((*m_cc_with_index)[m_cc_with_index->capacity()-1],
        m_free_list, FREE);
    m_free_list=old_size;
    // m_cc_with_index->upper_bound_ does not change.
  }

  void swap(Self& other)
  {
    // We should not swap m_cc_with_index, but only the content of the free list
    std::swap(m_free_list, other.m_free_list);
  }

  bool is_empty() const
  { return m_free_list==null_descriptor; }

  bool is_used(size_type i) const
  {
    CGAL_assertion(i<m_cc_with_index->capacity() && i!=null_descriptor);
    return static_type((*m_cc_with_index)[i])==USED;
  }

  // Push the ith element on the free list: it becomes free
  void push(size_type i)
  {
    CGAL_assertion(i<m_cc_with_index->capacity() && i!=null_descriptor);
    CGAL_assertion(is_used(i));
    static_set_val((*m_cc_with_index)[i], m_free_list, FREE);
    m_free_list=i;
  }

  size_type top() const
  {
    CGAL_assertion(!is_empty());
    return m_free_list;
  }

  // Pop one element from the free list (the top): it becomes used
  size_type pop()
  {
    CGAL_assertion(!is_empty());
    CGAL_assertion(!is_used(top()));
    size_type res=m_free_list;
    if(res>=m_cc_with_index->upper_bound())
    { m_cc_with_index->upper_bound_=res+1; }
    CGAL_assertion(res<m_cc_with_index->upper_bound());
    CGAL_assertion(res!=null_descriptor);
    CGAL_assertion(!is_used(res));
    static_set_type((*m_cc_with_index)[res], USED);
    m_free_list=static_get_val((*m_cc_with_index)[res]);
    return res;
  }

  void copy_special_data(const T& src, T& dest)
  { Traits::set_size_t(dest, Traits::size_t(src)); }

protected:
  // Definition of the bit squatting :
  // =================================
  // e is composed of a size_t and the big 1 bit.
  // value of the last bit as "Type" : 0  == reserved element; 1==free element.
  // When an element is free, the other bits represent the index of the
  // next free element.

  enum Type { USED = 0, FREE = 1 };

  static const int nbbits_size_type_m1 = sizeof(size_type)*8 - 1;
  static const size_type mask_type = ((size_type)-1)-(((size_type)-1)/2);

  // Get the type of the pointee.
  static Type static_type(const T& e)
  // TODO check if this is ok for little and big endian
  { return (Type) ((Traits::size_t(e) & mask_type)>>(nbbits_size_type_m1)); }

  // get the value of the element (removing the used bit)
  static size_type static_get_val(const T& e)
  { return (Traits::size_t(e) & ~mask_type); }

  // set the value of the element and its type
  static void static_set_type(T& e, Type t)
  { Traits::set_size_t(e, static_get_val(e) |
                       ( ((size_type)t) <<(nbbits_size_type_m1))); }

  // set the value of the element and its type
  static void static_set_val(T& e, size_type v, Type t)
  { Traits::set_size_t(e, v | ( ((size_type)t) <<(nbbits_size_type_m1))); }

protected:
  CC_with_index* m_cc_with_index;
  size_type      m_free_list;
};

// Index class
// templated by a compact container with index, to differentiate the different class
// (for example to allow two methods with same name and two different Index as parameter)
template<class CCWI, class Index_type, Index_type null_descriptor>
class Index_for_cc_with_index
{
public:
  using Self=Index_for_cc_with_index<CCWI, Index_type, null_descriptor>;
  using size_type=Index_type;

  /// Constructor.
  Index_for_cc_with_index(size_type idx=null_descriptor)
    : m_idx(idx)
  {}

  /// Get the underlying index
  operator size_type() const
  { return m_idx; }

  // Constructor allowing to transform an index from one container to another
  template<typename Index2>
  Index_for_cc_with_index(const Index2& idx): m_idx(static_cast<size_t>(idx))
  {}

  /// Increment the internal index. This operations does not
  /// guarantee that the index is valid or undeleted after the
  /// increment.
  Self& operator++() { ++m_idx; return *this; }
  /// Decrement the internal index. This operations does not
  /// guarantee that the index is valid or undeleted after the
  /// decrement.
  Self& operator--() { --m_idx; return *this; }

  /// Increment the internal index. This operations does not
  /// guarantee that the index is valid or undeleted after the
  /// increment.
  Self operator++(int) { Self tmp(*this); ++m_idx; return tmp; }
  /// Decrement the internal index. This operations does not
  /// guarantee that the index is valid or undeleted after the
  /// decrement.
  Self operator--(int) { Self tmp(*this); --m_idx; return tmp; }

  size_type for_compact_container() const
  { return m_idx; }
  void for_compact_container(size_type v)
  { m_idx=v; }

private:
  size_type m_idx;
};

namespace internal
{
struct Index_hash_function {
  typedef std::size_t result_type;
  template <class H>
  std::size_t operator() (const H& h) const {
    return h;
  }
};
}

template < class T, class Allocator_, class Increment_policy,
           class IndexType = std::size_t >
class Compact_container_with_index
{
  typedef Allocator_                                Al;
  typedef Increment_policy                          Incr_policy;
  typedef typename Default::Get< Al, CGAL_ALLOCATOR(T) >::type Allocator;
  typedef Compact_container_with_index<T, Al, Increment_policy, IndexType> Self;
public:
  typedef T                                         value_type;
  typedef IndexType                                 size_type;
  typedef Allocator                                 allocator_type;
  typedef value_type&                               reference;
  typedef const value_type&                         const_reference;

  typedef typename std::allocator_traits<Allocator>::pointer         pointer;
  typedef typename std::allocator_traits<Allocator>::const_pointer   const_pointer;
  typedef typename std::allocator_traits<Allocator>::difference_type difference_type;

  typedef internal::CC_iterator_with_index<Self, false> iterator;
  typedef internal::CC_iterator_with_index<Self, true>  const_iterator;
  typedef std::reverse_iterator<iterator>           reverse_iterator;
  typedef std::reverse_iterator<const_iterator>     const_reverse_iterator;

  using TFree_list_management=Free_list_management
                             <Self, CGAL::Tag_false, CGAL::Tag_false>;
  static const size_type null_descriptor=TFree_list_management::null_descriptor;
  using Index=Index_for_cc_with_index<Self, IndexType, null_descriptor>;

  friend class internal::CC_iterator_with_index<Self, false>;
  friend class internal::CC_iterator_with_index<Self, true>;
  friend TFree_list_management;

  template<unsigned int first_block_size_, unsigned int block_size_increment>
  friend struct Addition_size_policy;
  template<unsigned int k> friend struct Constant_size_policy_for_cc_with_size;
  template<unsigned int k>
  friend struct Multiply_by_two_policy_for_cc_with_size;

  explicit Compact_container_with_index(const Allocator &a = Allocator())
    : alloc(a),
      free_list(this)
  { init(); }

  template < class InputIterator >
  Compact_container_with_index(InputIterator first, InputIterator last,
                               const Allocator & a = Allocator())
  : alloc(a),
    free_list(this)
  {
    init();
    std::copy(first, last, CGAL::inserter(*this));
  }

  // The copy constructor and assignment operator preserve the iterator order
  Compact_container_with_index(const Compact_container_with_index &c)
  : alloc(c.get_allocator()),
    free_list(this)
  {
    init();
    next_delta=c.next_delta;
    std::copy(c.begin(), c.end(), CGAL::inserter(*this));
  }

  Compact_container_with_index(Compact_container_with_index&& c) noexcept
    : alloc(c.get_allocator()),
      free_list(*this)
  {  c.swap(*this); }

  Compact_container_with_index &
  operator=(const Compact_container_with_index &c)
  {
    if (&c != this) {
      Self tmp(c);
      swap(tmp);
    }
    return *this;
  }

  Compact_container_with_index & operator=(Compact_container_with_index&& c) noexcept
  {
    Self tmp(std::move(c));
    tmp.swap(*this);
    return *this;
  }

  ~Compact_container_with_index()
  { clear(); }

  void swap(Self &c)
  {
    std::swap(alloc, c.alloc);
    std::swap(capacity_, c.capacity_);
    std::swap(upper_bound_, c.upper_bound_);
    std::swap(size_, c.size_);
    std::swap(next_delta, c.next_delta);
    std::swap(all_items, c.all_items);
    free_list.swap(c.free_list);
  }

  bool is_used(size_type i) const
  { return free_list.is_used(i); }

  const T& operator[] (size_type i) const
  {
    CGAL_assertion(all_items!=nullptr && i<capacity_);
    return all_items[i];
  }

  T& operator[] (size_type i)
  {
    CGAL_assertion(all_items!=nullptr && i<capacity_);
    return all_items[i];
  }

  iterator begin() { return empty()?end():iterator(this, 0, 0); }
  iterator end()   { return iterator(this, upper_bound()); }

  const_iterator begin() const { return empty()?end():const_iterator(this, 0, 0); }
  const_iterator end()   const { return const_iterator(this, upper_bound()); }

  reverse_iterator rbegin() { return reverse_iterator(end()); }
  reverse_iterator rend()   { return reverse_iterator(begin()); }

  const_reverse_iterator
  rbegin() const { return const_reverse_iterator(end()); }
  const_reverse_iterator
  rend()   const { return const_reverse_iterator(begin()); }

  // Compute the index of a given pointer to an element of the compact container.
  size_type compute_index(const_pointer value) const
  {
    if (value >=all_items && value < (all_items+capacity_))
    {
      return (value-all_items);
    }
    return 0;
  }

  iterator index_to(size_type value) {
    return iterator(this, value);
  }
  const_iterator index_to(size_type value) const {
    return const_iterator(this, value);
  }

  // Boost.Intrusive interface
  iterator iterator_to(reference value) {
    return iterator(this, compute_index(&value));
  }
  const_iterator iterator_to(const_reference value) const {
    return const_iterator(this, compute_index(&value));
  }

  // Special insert methods that construct the objects in place
  // (just forward the arguments to the constructor, to optimize a copy).
  template < typename... Args >
  Index emplace(const Args&... args)
  {
    if (free_list.is_empty())
    { increase_size(); }

    Index ret=free_list.pop();
    T& e=operator[](ret);
    //std::allocator_traits<allocator_type>::construct(alloc, &e, args...);
    new (&e) value_type(args...);
    ++size_;
    return ret;
  }

  Index insert(const T &t)
  {
    if (free_list.is_empty())
    { increase_size(); }

    Index ret=free_list.pop();
    T& e=operator[](ret);
    //std::allocator_traits<allocator_type>::construct(alloc, &e, t);
    new (&e) value_type(t);
    ++size_;
    return ret;
  }

  template < class InputIterator >
  void insert(InputIterator first, InputIterator last)
  {
    for (; first != last; ++first)
    { insert(*first); }
  }

  template < class InputIterator >
  void assign(InputIterator first, InputIterator last)
  {
    clear();
    insert(first, last);
  }

  void erase(Index i)
  {
    CGAL_precondition(is_used(i));
    T& e=operator[](i);
    std::allocator_traits<allocator_type>::destroy(alloc, &e);
    //e.~T();
    free_list.push(i);
    --size_;
  }

  void erase(iterator first, iterator last) {
    while (first != last)
    { erase(first++); }
  }

  void clear();

  // Merge the content of d into *this. d gets cleared.
  // TODO? void merge(Self &d);

  size_type size() const
  {
    CGAL_expensive_assertion(size_ ==
                             (size_type) std::distance(begin(), end()));
    return size_;
  }

  size_type max_size() const
  { return std::allocator_traits<allocator_type>::max_size(alloc); }

  size_type capacity() const
  { return capacity_; }

  size_type upper_bound() const
  { return upper_bound_; }

  // void resize(size_type sz, T c = T()); // TODO  makes sense ???

  bool empty() const
  { return size_==0; }

  allocator_type get_allocator() const
  { return alloc; }

  size_type index(const_iterator cit) const
  { return static_cast<size_type>(cit); }

  size_type index(Index idx) const
  { return static_cast<size_type>(idx); }

  // Returns whether the iterator "cit" is in the range [begin(), end()].
  // This function is mostly useful for purposes of efficient debugging at
  // higher levels.
  bool owns(const_iterator cit) const
  {
     if (cit==end())
     { return true; }

    const_pointer c=&*cit;
    if (c>=all_items && c<(all_items+capacity_))
    { return is_used(cit); }
    return false;
  }

  bool owns(size_type i) const
  { return i<capacity() && is_used(i); }

  bool owns_dereferenceable(const_iterator cit) const
  { return cit!=end() && owns(cit); }

  CGAL_DEPRECATED bool owns_dereferencable(const_iterator cit) const
  { return owns_dereferenceable(cit); }

  /** Reserve method to ensure that the capacity of the Compact_container be
   * greater or equal than a given value n.
   */
  void reserve(size_type n)
  {
    if(capacity_>=n) return;
    increase_size();
  }

private:

  void increase_size();

  void init()
  {
    next_delta  =Incr_policy::first_block_size;
    capacity_   =0;
    upper_bound_=0;
    size_       =0;
    all_items   =nullptr;
    free_list.init();
  }

  allocator_type        alloc;        // allocator
  size_type             capacity_;    // capacity of the block
  size_type             upper_bound_; // <=capacity_ all indices >=upper_bound are unused
  size_type             size_;        // Number of elements allocated (used)
  size_type             next_delta;   // Size of the delta for the block size
  pointer               all_items;    // pointer to the block
  TFree_list_management free_list;    // free list management
};

template < class T, class Allocator, class Increment_policy, class IndexType >
void Compact_container_with_index<T, Allocator, Increment_policy, IndexType>::
clear()
{
  for (size_type i=0; i<capacity_; ++i)
  { if ( is_used(i) )
    { std::allocator_traits<allocator_type>::destroy(alloc, &operator[](i)); }
  }
  std::allocator_traits<allocator_type>::deallocate(alloc, all_items, capacity_);
  all_items=nullptr;

  init();
}

template < class T, class Allocator, class Increment_policy, class IndexType >
void Compact_container_with_index<T, Allocator, Increment_policy, IndexType>::
increase_size()
{
  size_type oldcapacity=capacity_;
  capacity_+=next_delta;

  pointer all_items2=
      std::allocator_traits<allocator_type>::allocate(alloc, capacity_);
  for (size_type index=0; index<oldcapacity; ++index)
  {
    if(is_used(index))
    {
      std::allocator_traits<allocator_type>::construct
          (alloc, &(all_items2[index]), std::move(all_items[index]));
      //new (&all_items2[index]) value_type(all_items[index]);
      std::allocator_traits<allocator_type>::destroy(alloc, &(all_items[index]));
    }
    else
    { free_list.copy_special_data(all_items[index], all_items2[index]); }
  }
  std::swap(all_items, all_items2);
  std::allocator_traits<allocator_type>::deallocate(alloc,
                                                    all_items2, oldcapacity);
  free_list.increase_to(oldcapacity);

  // Increase the block_size for the next time.
  Increment_policy::increase_size(*this);
}

namespace internal {

  // **********************************************************************
  // Iterator if we use index: in this case, an iterator use one more
  // data member of type size_type: memory footprint is more important than
  // iterator without index. However such iterator is not supposed to be used
  // in function parameters, to store handles through elements...
  // We must use indices for that.
  template < class DSC, bool Const >
  class CC_iterator_with_index
  {
    typedef typename DSC::iterator                    iterator;
    typedef CC_iterator_with_index<DSC, Const>        Self;

    friend class CC_iterator_with_index<DSC, true>;
    friend class CC_iterator_with_index<DSC, false>;

  public:
    typedef typename DSC::value_type                  value_type;
    typedef typename DSC::size_type                   size_type;
    typedef typename DSC::difference_type             difference_type;
    typedef std::conditional_t< Const, const value_type*,
                                       value_type*>   pointer;
    typedef std::conditional_t< Const, const value_type&,
                                       value_type&>  reference;
    typedef std::bidirectional_iterator_tag           iterator_category;

    typedef std::conditional_t< Const, const DSC*, DSC*> cc_pointer;

    CC_iterator_with_index(): m_ptr_to_cc(nullptr),
      m_index(0)
    {}

    // Either a harmless copy-ctor,
    // or a conversion from iterator to const_iterator.
    CC_iterator_with_index(const iterator &it): m_ptr_to_cc(it.m_ptr_to_cc),
      m_index(it.m_index)
    {
      CGAL_assertion(m_index<=m_ptr_to_cc->upper_bound());
    }

    // Same for assignment operator
    CC_iterator_with_index& operator=(const iterator &it)
    {
      m_ptr_to_cc=it.m_ptr_to_cc;
      m_index=it.m_index;
      return *this;
    }

    operator size_type() const
    { return m_index; }

    size_type get_current() const
    { return m_index; }

  protected:
    void set_current(size_type dh)
    { m_index=dh; }

  protected:

    // Only Compact_container should access these constructors.
    //template<class,class,class>
    friend class Compact_container_with_index<value_type, typename DSC::Al,
        typename DSC::Incr_policy,
        typename DSC::size_type>;
    cc_pointer m_ptr_to_cc;
    size_type m_index;

    // For begin()
    CC_iterator_with_index(cc_pointer ptr, int, int): m_ptr_to_cc(ptr),
      m_index(0)
    {
      if(!m_ptr_to_cc->is_used(m_index))
      { increment(); }
    }

    // Construction from raw pointer and for end().
    CC_iterator_with_index(cc_pointer ptr, size_type index): m_ptr_to_cc(ptr),
      m_index(index)
    {}

    // NB: in case empty container, begin==end==upper_bound.
    void increment()
    {
      // It's either pointing to end(), or valid.
      CGAL_assertion_msg(m_ptr_to_cc!=nullptr,
         "Incrementing a singular iterator or an empty container iterator ?");
      CGAL_assertion_msg(m_index<m_ptr_to_cc->upper_bound(),
         "Incrementing end() ?");

      // If it's not end(), then it's valid, we can do ++.
      do
      { ++m_index; }
      while(m_index<m_ptr_to_cc->upper_bound() &&
            (!m_ptr_to_cc->is_used(m_index)));
    }

    void decrement()
    {
      // It's either pointing to end(), or valid.
      CGAL_assertion_msg(m_ptr_to_cc!=nullptr,
         "Decrementing a singular iterator or an empty container iterator ?");
      CGAL_assertion_msg(m_index>0, "Decrementing begin() ?");

      // If it's not begin(), then it's valid, we can do --.
      do
      {
        CGAL_assertion(m_index>0);
        --m_index;
      }
      while(!m_ptr_to_cc->is_used(m_index));
    }

  public:

    Self & operator++()
    { increment(); return *this; }

    Self & operator--()
    { decrement(); return *this; }

    Self operator++(int) { Self tmp(*this); ++(*this); return tmp; }
    Self operator--(int) { Self tmp(*this); --(*this); return tmp; }

    reference operator*() const { return ((*m_ptr_to_cc)[m_index]); }

    pointer   operator->() const { return &((*m_ptr_to_cc)[m_index]); }

    bool is_end() const { return m_index>=m_ptr_to_cc->upper_bound(); }

    // Can itself be used for bit-squatting.
    size_type for_compact_container() const
    { return m_index; }
    void for_compact_container(size_type v)
    { m_index=v; }

    template<class ADSC,bool AC1,bool AC2>
    friend bool operator==(const CC_iterator_with_index<ADSC,AC1>&,
                           const CC_iterator_with_index<ADSC,AC2>&);

    template<class ADSC,bool AC1,bool AC2>
    friend bool operator!=(const CC_iterator_with_index<ADSC,AC1>&,
                           const CC_iterator_with_index<ADSC,AC2>&);
  };

  template < class DSC, bool Const1, bool Const2 >
  inline
  bool operator==(const CC_iterator_with_index<DSC, Const1> &rhs,
                  const CC_iterator_with_index<DSC, Const2> &lhs)
  {
    return rhs.m_ptr_to_cc==lhs.m_ptr_to_cc &&
      (rhs.m_index==lhs.m_index ||
       (rhs.is_end() && lhs.is_end()));
  }

  template < class DSC, bool Const1, bool Const2 >
  inline
  bool operator!=(const CC_iterator_with_index<DSC, Const1> &rhs,
                  const CC_iterator_with_index<DSC, Const2> &lhs)
  {
    return rhs.m_ptr_to_cc!=lhs.m_ptr_to_cc ||
                            (rhs.m_index!=lhs.m_index &&
                                          (!rhs.is_end() || !lhs.is_end()));
  }

} // namespace internal

} //namespace CGAL

namespace std
{

#ifndef CGAL_CFG_NO_STD_HASH

template <class CCWI, class Index_type, Index_type null_descriptor>
struct hash<CGAL::Index_for_cc_with_index<CCWI, Index_type, null_descriptor>>:
    public CGAL::cpp98::unary_function<CGAL::Index_for_cc_with_index<CCWI, Index_type, null_descriptor>,
    std::size_t>
{
  std::size_t operator()(const CGAL::Index_for_cc_with_index<CCWI, Index_type, null_descriptor>& idx) const
  { return CGAL::internal::Index_hash_function()(idx); }
};
#endif // CGAL_CFG_NO_STD_HASH

} // namespace std

#endif // CGAL_COMPACT_CONTAINER_WITH_INDEX_H
