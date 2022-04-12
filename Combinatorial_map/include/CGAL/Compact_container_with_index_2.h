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
// $URL$
// $Id$
//
// Author(s)     : Sylvain Pion

#ifndef CGAL_COMPACT_CONTAINER_WITH_INDEX_2_H
#define CGAL_COMPACT_CONTAINER_WITH_INDEX_2_H

#include <CGAL/Compact_container_with_index.h>

// An STL like container with the following properties :
// - to achieve compactness, it requires access to a pointer stored in T,
//   specified by a traits.  This pointer is supposed to be 4 bytes aligned
//   when the object is alive, otherwise, the container uses the 2 least
//   significant bits to store information in the pointer.
// - Ts are allocated in arrays of increasing size, which are linked together
//   by their first and last element.
// - the iterator looks at the famous 2 bits to know if it has to deal with
//   a free/used/boundary element.

// TODO :
// - Add .resize() (and proper copy of capacity_).
// - Add preconditions in input that real pointers need to have clean bits.
//   Also for the allocated memory alignment, and sizeof().
// - Do a benchmark before/after.
// - Check the end result with Valgrind.
// - The bit squatting mechanism will be reused for the conflict flag, maybe
//   it could be put out of the class.

// TODO low priority :
// - rebind<> the allocator
// - Exception safety guarantees
// - Thread safety guarantees
// - std requirements on iterators says all defined operations are constant
//   time amortized (it's not true here, maybe it could be with some work...)
// - all this is expected especially when there are not so many free objects
//   compared to the allocated elements.
// - Should block_size be selectable/hintable by .reserve() ?
// - would be nice to have a temporary_free_list (still active elements, but
//   which are going to be freed soon).  Probably it prevents compactness.
// - eventually something to copy this data structure, providing a way to
//   update the pointers (give access to a hash_map, at least a function that
//   converts an old pointer to the new one ?).  Actually it doesn't have to
//   be stuck to a particular DS, because for a list it's useful too...
// - Currently, end() can be invalidated on insert() if a new block is added.
//   It would be nice to fix this.  We could insert the new block at the
//   beginning instead ?  That would drop the property that iterator order
//   is preserved.  Maybe it's not a problem if end() is not preserved, after
//   all nothing is going to dereference it, it's just for comparing with
//   end() that it can be a problem.
//   Another way would be to have end() point to the end of an always
//   empty block (containing no usable element), and insert new blocks just
//   before this one.
//   Instead of having the blocks linked between them, the start/end pointers
//   could point back to the container, so that we can do more interesting
//   things (e.g. freeing empty blocks automatically) ?

namespace CGAL {

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
class Compact_container_with_index_2
{
  typedef Allocator_                                Al;
  typedef Increment_policy                          Incr_policy;
  typedef typename Default::Get< Al, CGAL_ALLOCATOR(T) >::type Allocator;
  typedef Compact_container_with_index_2 <T, Al, Increment_policy, IndexType> Self;
  typedef Compact_container_with_index_traits <T, IndexType>   Traits;
public:
  typedef T                                         value_type;
  typedef IndexType                                 size_type;
  typedef Allocator                                 allocator_type;
  typedef typename Allocator::reference             reference;
  typedef typename Allocator::const_reference       const_reference;
  typedef typename Allocator::pointer               pointer;
  typedef typename Allocator::const_pointer         const_pointer;
  typedef typename Allocator::difference_type       difference_type;
  typedef internal::CC_iterator_with_index<Self, false> iterator;
  typedef internal::CC_iterator_with_index<Self, true>  const_iterator;
  typedef std::reverse_iterator<iterator>           reverse_iterator;
  typedef std::reverse_iterator<const_iterator>     const_reverse_iterator;

  static const size_type bottom;

  class Index : public internal::MyIndex<T,size_type>
  {
  public:
    typedef typename Compact_container_with_index_2::size_type size_type;
    typedef internal::MyIndex<T,size_type> Base;

    Index(size_type idx=(std::numeric_limits<size_type>::max)()/2)
      : Base(idx)
    {}

    Index(const Index& idx): Base(idx)
    {}

    Index(const const_iterator& it) : Base(it)
    {}

    Index(const iterator& it) : Base(it)
    {}

    // Constructor allowing to transform an index from one container to another
    template<typename Index2>
    Index(const Index2& idx): Base(static_cast<size_t>(idx))
    {}
  };
  friend class internal::CC_iterator_with_index<Self, false>;
  friend class internal::CC_iterator_with_index<Self, true>;

  template<unsigned int first_block_size_, unsigned int block_size_increment>
  friend struct Addition_size_policy;
  template<unsigned int k> friend struct Constant_size_policy_for_cc_with_size;
  template<unsigned int k>
  friend struct Multiply_by_two_policy_for_cc_with_size;

  explicit Compact_container_with_index_2(const Allocator &a = Allocator())
  : alloc(a)
  {
    init();
  }

  template < class InputIterator >
  Compact_container_with_index_2(InputIterator first, InputIterator last,
                               const Allocator & a = Allocator())
  : alloc(a)
  {
    init();
    std::copy(first, last, CGAL::inserter(*this));
  }

  // The copy constructor and assignment operator preserve the iterator order
  Compact_container_with_index_2(const Compact_container_with_index_2 &c)
  : alloc(c.get_allocator())
  {
    init();
    block_size = c.block_size;
    std::copy(c.begin(), c.end(), CGAL::inserter(*this));
  }

  Compact_container_with_index_2(Compact_container_with_index_2&& c) noexcept
  : alloc(c.get_allocator())
  {
    c.swap(*this);
  }

  Compact_container_with_index_2 &
  operator=(const Compact_container_with_index_2 &c)
  {
    if (&c != this) {
      Self tmp(c);
      swap(tmp);
    }
    return *this;
  }

  Compact_container_with_index_2 & operator=(Compact_container_with_index_2&& c) noexcept
  {
    Self tmp(std::move(c));
    tmp.swap(*this);
    return *this;
  }

  ~Compact_container_with_index_2()
  {
    clear();
  }

  bool is_used(size_type i) const
  {
    return (type(i)==USED);
  }

  const T& operator[] (size_type i) const
  {
    CGAL_assertion(all_items!=NULL && i<capacity_);
    return all_items[i];
  }

  T& operator[] (size_type i)
  {
    CGAL_assertion(all_items!=NULL && i<capacity_);
    return all_items[i];
  }

  void swap(Self &c)
  {
    std::swap(alloc, c.alloc);
    std::swap(capacity_, c.capacity_);
    std::swap(size_, c.size_);
    std::swap(block_size, c.block_size);
    std::swap(free_list, c.free_list);
    std::swap(all_items, c.all_items);
  }

  iterator begin() { return iterator(this, 0, 0); }
  iterator end()   { return iterator(this, capacity_); }

  const_iterator begin() const { return const_iterator(this, 0, 0); }
  const_iterator end()   const { return const_iterator(this, capacity_); }

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
    if (free_list == bottom)
      allocate_new_block();

    Index ret = free_list;
    T& e = operator[](free_list);
    static_set_type(e, USED);
    free_list = static_get_val(e);
    //std::allocator_traits<allocator_type>::construct(alloc, &e, args...);
    new (&e) value_type(args...);
    ++size_;
    return ret;
  }

  Index insert(const T &t)
  {
    if (free_list == bottom)
      allocate_new_block();

    Index ret(free_list);
    T& e = operator[](free_list);
    static_set_type(e, USED);
    free_list = static_get_val(e);
    //std::allocator_traits<allocator_type>::construct(alloc, &e, t);
    new (&e) value_type(t);
    ++size_;
    return ret;
  }

  template < class InputIterator >
  void insert(InputIterator first, InputIterator last)
  {
    for (; first != last; ++first)
      insert(*first);
  }

  template < class InputIterator >
  void assign(InputIterator first, InputIterator last)
  {
    clear(); // erase(begin(), end()); // ?
    insert(first, last);
  }

  void erase(Index x)
  {
    CGAL_precondition(type(x) == USED);
    T& e = operator[](x);
    std::allocator_traits<allocator_type>::destroy(alloc, &e);
    //e.~T();
#ifndef CGAL_NO_ASSERTIONS
    std::memset(&e, 0, sizeof(T));
#endif
    put_on_free_list(x);
    --size_;
  }

  void erase(iterator first, iterator last) {
    while (first != last)
      erase(first++);
  }

  void clear();

  // Merge the content of d into *this.  d gets cleared.
  // The complexity is O(size(free list = capacity-size)).
  void merge(Self &d);

  size_type size() const
  {
    CGAL_expensive_assertion(size_ ==
                             (size_type) std::distance(begin(), end()));
    return size_;
  }

  size_type max_size() const
  {
    return std::allocator_traits<allocator_type>::max_size(alloc);
  }

  size_type capacity() const
  {
    return capacity_;
  }

  // void resize(size_type sz, T c = T()); // TODO  makes sense ???

  bool empty() const
  {
    return size_ == 0;
  }

  allocator_type get_allocator() const
  {
    return alloc;
  }

  size_type index(const_iterator cit) const
  { return static_cast<size_type>(cit); }

  size_type index(Index idx) const
  { return static_cast<size_type>(idx); }

  // Returns whether the iterator "cit" is in the range [begin(), end()].
  // Complexity : O(#blocks) = O(sqrt(capacity())).
  // This function is mostly useful for purposes of efficient debugging at
  // higher levels.
  bool owns(const_iterator cit) const
  {
     if (cit == end())
        return true;

    const_pointer c = &*cit;

    if (c >=all_items && c < (all_items+capacity_))
    {
        return type(cit) == USED;
    }
    return false;
  }

  bool owns_dereferencable(const_iterator cit) const
  {
    return cit != end() && owns(cit);
  }

  /** Reserve method to ensure that the capacity of the Compact_container be
   * greater or equal than a given value n.
   */
  void reserve(size_type n)
  {
    if ( capacity_>=n ) return;
// TODO
 /*   size_type lastblock = all_items.size();

    while ( capacity_<n )
    {
      pointer new_block = alloc.allocate(block_size);
      all_items.push_back(std::make_pair(new_block, block_size));
      capacity_ += block_size;
      // Increase the block_size for the next time.
      Increment_policy::increase_size(*this);
    }

    // Now we put all the new elements on freelist, starting from the last block
    // inserted and mark them free in reverse order, so that the insertion order
    // will correspond to the iterator order...
    // We don't touch the first and the last one.
    size_type curblock=all_items.size();
    size_type index = capacity_-1;
    do
    {
      --curblock; // We are sure we have at least create a new block
      for (size_type i = all_items[curblock].second-1; i >= 0; --i, --index)
        put_on_free_list(index);
    }
    while ( curblock>lastblock );*/
  }

private:

  void allocate_new_block();

  // Definition of the bit squatting :
  // =================================
  // e is composed of a size_t and the big 1 bit.
  // Here is the meaning of each of the 4 cases.
  //
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

  Type type(size_type e) const
  { return static_type(operator[](e)); }

  // Sets the pointer part and the type of the pointee.
  static void static_set_type(T& e, Type t)
  {
    // This out of range compare is always true and causes lots of
    // unnecessary warnings.
    // CGAL_precondition(0 <= t && t < 2);
    Traits::size_t(e) &= ( ~mask_type | ( ((size_type)t) <<(nbbits_size_type_m1) ) );
  }

  // get the value of the element (removing the two bits)
  static size_type static_get_val(const T& e)
  { return (Traits::size_t(e) & ~mask_type); }

  size_type get_val(size_type e) const
  { return static_get_val(operator[](e)); }

  // set the value of the element and its type
  static void static_set_val(T& e, size_type v, Type t)
  { Traits::size_t(e)=v | ( ((size_type)t) <<(nbbits_size_type_m1)); }

  void set_val(size_type e, size_type v, Type t)
  { static_set_val(operator[](e), v, t); }

  void put_on_free_list(size_type x)
  {
    set_val(x, free_list, FREE);
    free_list = x;
  }

  void init()
  {
    block_size = Incr_policy::first_block_size;
    capacity_  = 0;
    size_      = 0;
    free_list  = bottom;
    all_items  = NULL;
  }

  allocator_type   alloc;
  size_type        capacity_;
  size_type        size_;
  size_type        block_size;
  size_type        free_list;
  pointer          all_items;
};

template < class T, class Allocator, class Increment_policy, class IndexType >
const typename Compact_container_with_index_2<T, Allocator, Increment_policy, IndexType>::size_type
Compact_container_with_index_2<T, Allocator, Increment_policy, IndexType>::bottom =  (std::numeric_limits<typename Compact_container_with_index_2<T, Allocator, Increment_policy, IndexType>::size_type>::max)()/2;

/*template < class T, class Allocator, class Increment_policy, class IndexType >
void Compact_container_with_index<T, Allocator, Increment_policy, IndexType>::merge(Self &d)
{
  CGAL_precondition(&d != this);

  // Allocators must be "compatible" :
  CGAL_precondition(get_allocator() == d.get_allocator());

  // Concatenate the free_lists.
  if (free_list == bottom) {
    free_list = d.free_list;
  } else if (d.free_list != 0) {
    size_type e = free_list;
    while (get_val(e) != 0)
      e = get_val(e);
    set_val(e, d.free_list, FREE);
  }
  // Add the sizes.
  size_ += d.size_;
  // Add the capacities.
  capacity_ += d.capacity_;
  // It seems reasonnable to take the max of the block sizes.
  block_size = (std::max)(block_size, d.block_size);
  // Clear d.
  d.init();
}*/

template < class T, class Allocator, class Increment_policy, class IndexType >
void Compact_container_with_index_2<T, Allocator, Increment_policy, IndexType>::clear()
{
  for (size_type i=0; i<capacity_; ++i)
  {
    if ( is_used(i) ) alloc.destroy(&operator[](i));
  }

  std::allocator_traits<allocator_type>::deallocate(alloc, all_items, capacity_);
  all_items=nullptr;

  init();
}

template < class T, class Allocator, class Increment_policy, class IndexType >
void Compact_container_with_index_2<T, Allocator, Increment_policy, IndexType>::allocate_new_block()
{
  size_type oldcapacity=capacity_;
  capacity_ += block_size;

  pointer all_items2=
      std::allocator_traits<allocator_type>::allocate(alloc, capacity_);
  for (size_type index=0; index<oldcapacity; ++index)
  {
    if(is_used(index))
    {
      std::allocator_traits<allocator_type>::construct(alloc, &(all_items2[index]),
                                                       std::move(all_items[index]));
      alloc.destroy(&(all_items[index]));
    }
    else
    {
      static_set_val(all_items2[index], static_get_val(all_items[index]), FREE);
    }
  }
  std::swap(all_items, all_items2);
  std::allocator_traits<allocator_type>::deallocate(alloc, all_items2, oldcapacity);

  // We mark them free in reverse order, so that the insertion order
  // will correspond to the iterator order...
  for (size_type index = capacity_-1; index>oldcapacity; --index)
  { put_on_free_list(index); }
  put_on_free_list(oldcapacity);

  // Increase the block_size for the next time.
  Increment_policy::increase_size(*this);
}

} //namespace CGAL

#endif // CGAL_COMPACT_CONTAINER_WITH_INDEX_2_H
