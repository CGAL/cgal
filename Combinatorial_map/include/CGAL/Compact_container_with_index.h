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

#ifndef CGAL_COMPACT_CONTAINER_WITH_INDEX_H
#define CGAL_COMPACT_CONTAINER_WITH_INDEX_H

#include <CGAL/Compact_container.h>

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

template<unsigned int k>
struct Multiply_by_two_policy_for_cc_with_size
{
  static const unsigned int first_block_size = k;

  template<typename Compact_container>
  static void increase_size(Compact_container& cc)
  { cc.block_size=cc.capacity_; }
};

template<unsigned int k>
struct Constant_size_policy_for_cc_with_size
{
  static const unsigned int first_block_size = k;

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

template <class, class, class, class>
class Compact_container_with_index_2;

namespace internal {
  template < class DSC, bool Const>
  class CC_iterator_with_index;

  template < class T, class ST >
  class MyIndex;
}

template < class T, class Allocator_, class Increment_policy, class IndexType = std::size_t >
class Compact_container_with_index
{
  typedef Allocator_                                Al;
  typedef Increment_policy                          Incr_policy;
  typedef typename Default::Get< Al, CGAL_ALLOCATOR(T) >::type Allocator;
  typedef Compact_container_with_index <T, Al, Increment_policy, IndexType> Self;
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
    typedef typename Compact_container_with_index::size_type size_type;
    typedef internal::MyIndex<T,size_type> Base;

    explicit Index(size_type idx=(std::numeric_limits<size_type>::max)()/2)
      : Base(idx)
    {}

    Index(const const_iterator& it) : Base(it)
    {}

    Index(const iterator& it) : Base(it)
    {}
  };

  friend class internal::CC_iterator_with_index<Self, false>;
  friend class internal::CC_iterator_with_index<Self, true>;

  template<unsigned int first_block_size_, unsigned int block_size_increment>
    friend struct Addition_size_policy;
  template<unsigned int k> friend struct Constant_size_policy_for_cc_with_size;

  explicit Compact_container_with_index(const Allocator &a = Allocator())
  : alloc(a)
  {
    init();
  }

  template < class InputIterator >
  Compact_container_with_index(InputIterator first, InputIterator last,
                               const Allocator & a = Allocator())
  : alloc(a)
  {
    init();
    std::copy(first, last, CGAL::inserter(*this));
  }

  // The copy constructor and assignment operator preserve the iterator order
  Compact_container_with_index(const Compact_container_with_index &c)
  : alloc(c.get_allocator())
  {
    init();
    block_size = c.block_size;
    std::copy(c.begin(), c.end(), CGAL::inserter(*this));
  }

  Compact_container_with_index &
  operator=(const Compact_container_with_index &c)
  {
    if (&c != this) {
      Self tmp(c);
      swap(tmp);
    }
    return *this;
  }

  ~Compact_container_with_index()
  {
    clear();
  }

  bool is_used(size_type i) const
  {
    return (type(i)==USED);
  }

  const T& operator[] (size_type i) const
  {
    typename Self::size_type block_number, index_in_block;
    Increment_policy::template get_index_and_block<Self>(i,
                                                         index_in_block,
                                                         block_number);
    return all_items[block_number].first[index_in_block];
  }

  T& operator[] (size_type i)
  {
    typename Self::size_type block_number, index_in_block;
    Increment_policy::template get_index_and_block<Self>(i,
                                                         index_in_block,
                                                         block_number);
    return all_items[block_number].first[index_in_block];
  }

  void swap(Self &c)
  {
    std::swap(alloc, c.alloc);
    std::swap(capacity_, c.capacity_);
    std::swap(size_, c.size_);
    std::swap(block_size, c.block_size);
    std::swap(free_list, c.free_list);
    all_items.swap(c.all_items);
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
    size_type res = 0;
    for (typename All_items::const_iterator it = all_items.begin(),
         itend = all_items.end(); it != itend; ++it)
    {
      const_pointer p = it->first;
      size_type     s = it->second;
      if (value >= p && value < (p+s))
      {
        return res + (value-p);
      }
      res += s;
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
#ifndef CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES
  template < typename... Args >
  size_type emplace(const Args&... args)
  {
    if (free_list == bottom)
      allocate_new_block();

    size_type ret = free_list;
    T& e = operator[](free_list);
    static_set_type(e, USED);
    free_list = static_get_val(e);
    new (&e) value_type(args...);
    ++size_;
    return ret;
  }
#else
  // inserts a default constructed item.
  size_type emplace()
  {
    if (free_list == bottom)
      allocate_new_block();

    size_type ret = free_list;
    T& e = operator[](free_list);
    static_set_type(e, USED);
    free_list = static_get_val(e);
    new (&e) value_type();
    ++size_;

    return ret;
  }

  template < typename T1 >
  size_type
  emplace(const T1 &t1)
  {
    if (free_list == bottom)
      allocate_new_block();

    size_type ret = free_list;
    T& e = operator[](free_list);
    static_set_type(e, USED);
    free_list = static_get_val(e);
    new (&e) value_type(t1);
    ++size_;
    return ret;
  }

  template < typename T1, typename T2 >
  size_type
  emplace(const T1 &t1, const T2 &t2)
  {
    if (free_list == bottom)
      allocate_new_block();

    size_type ret = free_list;
    T& e = operator[](free_list);
    static_set_type(e, USED);
    free_list = static_get_val(e);
    new (&e) value_type(t1, t2);
    ++size_;
    return ret;
  }

  template < typename T1, typename T2, typename T3 >
  size_type
  emplace(const T1 &t1, const T2 &t2, const T3 &t3)
  {
    if (free_list == bottom)
      allocate_new_block();

    size_type ret = free_list;
    T& e = operator[](free_list);
    static_set_type(e, USED);
    free_list = static_get_val(e);
    new (&e) value_type(t1, t2, t3);
    ++size_;
    return ret;
  }

  template < typename T1, typename T2, typename T3, typename T4 >
  size_type
  emplace(const T1 &t1, const T2 &t2, const T3 &t3, const T4 &t4)
  {
    if (free_list == bottom)
      allocate_new_block();

    size_type ret = free_list;
    T& e = operator[](free_list);
    static_set_type(e, USED);
    free_list = static_get_val(e);
    new (&e) value_type(t1, t2, t3, t4);
    ++size_;
    return ret;
   }

  template < typename T1, typename T2, typename T3, typename T4, typename T5 >
  size_type
  emplace(const T1 &t1, const T2 &t2, const T3 &t3, const T4 &t4,
          const T5 &t5)
  {
    if (free_list == bottom)
      allocate_new_block();

    size_type ret = free_list;
    T& e = operator[](free_list);
    static_set_type(e, USED);
    free_list = static_get_val(e);
    new (&e) value_type(t1, t2, t3, t4, t5);
    ++size_;
    return ret;
  }

  template < typename T1, typename T2, typename T3, typename T4,
             typename T5, typename T6 >
  size_type
  emplace(const T1 &t1, const T2 &t2, const T3 &t3, const T4 &t4,
          const T5 &t5, const T6 &t6)
  {
    if (free_list == bottom)
      allocate_new_block();

    size_type ret = free_list;
    T& e = operator[](free_list);
    static_set_type(e, USED);
    free_list = static_get_val(e);
    new (&e) value_type(t1, t2, t3, t4, t5, t6);
    ++size_;
    return ret;
 }

  template < typename T1, typename T2, typename T3, typename T4,
             typename T5, typename T6, typename T7 >
  size_type
  emplace(const T1 &t1, const T2 &t2, const T3 &t3, const T4 &t4,
          const T5 &t5, const T6 &t6, const T7 &t7)
  {
    if (free_list == bottom)
      allocate_new_block();

    size_type ret = free_list;
    T& e = operator[](free_list);
    static_set_type(e, USED);
    free_list = static_get_val(e);
    new (&e) value_type(t1, t2, t3, t4, t5, t6, t7);
    ++size_;
    return ret;
  }

  template < typename T1, typename T2, typename T3, typename T4,
             typename T5, typename T6, typename T7, typename T8 >
  size_type
  emplace(const T1 &t1, const T2 &t2, const T3 &t3, const T4 &t4,
          const T5 &t5, const T6 &t6, const T7 &t7, const T8 &t8)
  {
    if (free_list == bottom)
      allocate_new_block();

    size_type ret = free_list;
    T& e = operator[](free_list);
    static_set_type(e, USED);
    free_list = static_get_val(e);
    new (&e) value_type(t1, t2, t3, t4, t5, t6, t7, t8);
    ++size_;
    return ret;
  }
#endif // CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES

  size_type insert(const T &t)
  {
    if (free_list == bottom)
      allocate_new_block();

    size_type ret = free_list;
    T& e = operator[](free_list);
    static_set_type(e, USED);
    free_list = static_get_val(e);
    alloc.construct(&e, t);
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

  void erase(size_type x)
  {
    CGAL_precondition(type(x) == USED);
    T& e = operator[](x);
    alloc.destroy(&e);
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
    return alloc.max_size();
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

  // Returns whether the iterator "cit" is in the range [begin(), end()].
  // Complexity : O(#blocks) = O(sqrt(capacity())).
  // This function is mostly useful for purposes of efficient debugging at
  // higher levels.
  bool owns(const_iterator cit) const
  {
    // We use the block structure to provide an efficient version :
    // we check if the address is in the range of each block,
    // and then test whether it is valid (not a free element).

    if (cit == end())
      return true;

    const_pointer c = &*cit;

    for (typename All_items::const_iterator it = all_items.begin(),
         itend = all_items.end(); it != itend; ++it) {
      const_pointer p = it->first;
      size_type s = it->second;

      // Are we in the address range of this block (excluding first and last
      // elements) ?
      if (c>=p && c<(p+s))
      {
        return type(cit) == USED;
      }
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

    size_type lastblock = all_items.size();

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
    while ( curblock>lastblock );
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
  // TODO check if this is ok for little and big endian
  static Type static_type(const T& e)
  {
    return (Type) ((Traits::size_t(e).get_idx() & mask_type)>>(nbbits_size_type_m1));
  }

  Type type(size_type e) const
  { return static_type(operator[](e)); }

  // get the value of the element (removing the used bit)
  static size_type static_get_val(const T& e)
  { return (Traits::size_t(e).get_idx() & ~mask_type); }

  size_type get_val(size_type e) const
  { return static_get_val(operator[](e)); }

  // set the value of the element and its type
  static void static_set_val(T& e, size_type v, Type t)
  { Traits::set_size_t(e, v | ( ((size_type)t) <<(nbbits_size_type_m1))); }

  // Sets the pointer part and the type of the pointee.
  static void static_set_type(T& e, Type t)
  { static_set_val(e, Traits::size_t(e)&~mask_type, t); }

  void set_val(size_type e, size_type v, Type t)
  { static_set_val(operator[](e), v, t); }

  void put_on_free_list(size_type x)
  {
    set_val(x, free_list, FREE);
    free_list = x;
  }

  // We store a vector of pointers to all allocated blocks and their sizes.
  // Knowing all pointers, we don't have to walk to the end of a block to reach
  // the pointer to the next block.
  // Knowing the sizes allows to deallocate() without having to compute the size
  // by walking through the block till its end.
  // This opens up the possibility for the compiler to optimize the clear()
  // function considerably when has_trivial_destructor<T>.
  typedef std::vector<std::pair<pointer, size_type> >  All_items;

  void init()
  {
    block_size = Incr_policy::first_block_size;
    capacity_  = 0;
    size_      = 0;
    free_list  = bottom;
    all_items  = All_items();
  }

  allocator_type   alloc;
  size_type        capacity_;
  size_type        size_;
  size_type        block_size;
  size_type        free_list;
  All_items        all_items;
};
 template < class T, class Allocator, class Increment_policy, class IndexType >
 const typename Compact_container_with_index<T, Allocator, Increment_policy, IndexType>::size_type Compact_container_with_index<T, Allocator, Increment_policy, IndexType>::bottom = (std::numeric_limits<typename Compact_container_with_index<T, Allocator, Increment_policy, IndexType>::size_type>::max)()/2;

/*template < class T, class Allocator, class Increment_policy >
void Compact_container_with_index<T, Allocator, Increment_policy>::merge(Self &d)
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
void Compact_container_with_index<T, Allocator, Increment_policy, IndexType>::clear()
{
  for (size_type i=0; i<capacity_; ++i)
  {
    if ( is_used(i) ) alloc.destroy(&operator[](i));
  }

  for (typename All_items::iterator it = all_items.begin(),
       itend = all_items.end(); it != itend; ++it)
  {
    alloc.deallocate(it->first, it->second);
  }

  init();
}

template < class T, class Allocator, class Increment_policy, class IndexType >
void Compact_container_with_index<T, Allocator, Increment_policy, IndexType>::allocate_new_block()
{
  pointer new_block = alloc.allocate(block_size);
  all_items.push_back(std::make_pair(new_block, block_size));
  // We mark them free in reverse order, so that the insertion order
  // will correspond to the iterator order...
  for (size_type index = capacity_+block_size-1; index>capacity_; --index)
    put_on_free_list(index);

  put_on_free_list(capacity_);
  capacity_ += block_size;
  // Increase the block_size for the next time.
  Increment_policy::increase_size(*this);
}

namespace internal {

  // **********************************************************************
  // specialization if we use index: in this case, an iterator use one more
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
    typedef typename boost::mpl::if_c< Const, const value_type*,
                                       value_type*>::type pointer;
    typedef typename boost::mpl::if_c< Const, const value_type&,
                                       value_type&>::type reference;
    typedef std::bidirectional_iterator_tag           iterator_category;

    typedef typename boost::mpl::if_c< Const, const DSC*, DSC*>::type
    cc_pointer;

    // the initialization with NULL is required by our Handle concept.
    CC_iterator_with_index() : m_ptr_to_cc(NULL),
      m_index(0)
    {}

    // Either a harmless copy-ctor,
    // or a conversion from iterator to const_iterator.
    CC_iterator_with_index (const iterator &it) : m_ptr_to_cc(it.m_ptr_to_cc),
      m_index(it.m_index)
    {}

    // Same for assignment operator (otherwise MipsPro warns)
    CC_iterator_with_index & operator= (const iterator &it)
    {
      m_ptr_to_cc = it.m_ptr_to_cc;
      m_index = it.m_index;
      return *this;
    }

    // Construction from NULL
    CC_iterator_with_index (Nullptr_t CGAL_assertion_code(n)) :
      m_ptr_to_cc(NULL),
      m_index(0)
    { CGAL_assertion (n == NULL); }

    operator size_type() const
    { return m_index; }

    size_type get_current() const
    { return m_index; }

  protected:
    void set_current(size_type dh)
    { m_index =  dh; }

  protected:

    // Only Compact_container should access these constructors.
    //template<class,class,class>
    friend class Compact_container_with_index<value_type, typename DSC::Al,
                                              typename DSC::Incr_policy, typename DSC::size_type >;

    //template<class,class,class>
    friend class Compact_container_with_index_2<value_type, typename DSC::Al,
                                                typename DSC::Incr_policy, typename DSC::size_type>;
    template<class,class,class>
    friend class Compact_container_with_index_3;

    cc_pointer m_ptr_to_cc;
    size_type m_index;

    // For begin()
    CC_iterator_with_index(cc_pointer ptr, int, int) : m_ptr_to_cc(ptr),
      m_index(0)
    {
      if(m_ptr_to_cc->type(m_index) != DSC::USED)
      { increment(); }
    }

    // Construction from raw pointer and for end().
    CC_iterator_with_index(cc_pointer ptr, size_type index) : m_ptr_to_cc(ptr),
      m_index(index)
    {}

    // NB : in case empty container, begin == end == NULL.
    void increment()
    {
      // It's either pointing to end(), or valid.
      CGAL_assertion_msg(m_ptr_to_cc != NULL,
         "Incrementing a singular iterator or an empty container iterator ?");
      CGAL_assertion_msg(m_index < m_ptr_to_cc->capacity_,
         "Incrementing end() ?");

      // If it's not end(), then it's valid, we can do ++.
      do
      {
        ++m_index;
      }
      while ( m_index < m_ptr_to_cc->capacity_ &&
              (m_ptr_to_cc->type(m_index) != DSC::USED) );
    }

    void decrement()
    {
      // It's either pointing to end(), or valid.
      CGAL_assertion_msg(m_ptr_to_cc != NULL,
         "Decrementing a singular iterator or an empty container iterator ?");
      CGAL_assertion_msg(m_index>0, "Decrementing begin() ?");

      // If it's not begin(), then it's valid, we can do --.
      do
      {
        --m_index;
      }
      while ( m_ptr_to_cc->type(m_index) != DSC::USED);
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
    return rhs.m_ptr_to_cc == lhs.m_ptr_to_cc &&
      rhs.m_index == lhs.m_index;
  }

  template < class DSC, bool Const1, bool Const2 >
  inline
  bool operator!=(const CC_iterator_with_index<DSC, Const1> &rhs,
                  const CC_iterator_with_index<DSC, Const2> &lhs)
  {
    return rhs.m_ptr_to_cc != lhs.m_ptr_to_cc ||
      rhs.m_index != lhs.m_index;
  }

  // Comparisons with NULL are part of CGAL's Handle concept...
  /*  template < class DSC, bool Const >
  inline
  bool operator==(const CC_iterator_with_index<DSC, Const> &rhs,
                  Nullptr_t CGAL_assertion_code(n))
  {
    CGAL_assertion( n == NULL);
    return rhs.m_index == 0;
  }

  template < class DSC, bool Const >
  inline
  bool operator!=(const CC_iterator_with_index<DSC, Const> &rhs,
		  Nullptr_t CGAL_assertion_code(n))
  {
    CGAL_assertion( n == NULL);
    return rhs.m_index != 0;
    }*/

  template<typename T, typename IT >
  class MyIndex
  {
  public:
    template < class, class, class, class >
    friend class Compact_container_with_index;
    template < class, class, class, class >
    friend class Compact_container_with_index_2;
    
    typedef MyIndex<T,IT> Self;
    typedef IT size_type;

    /// Constructor. Default construction creates a kind of "NULL" index.
    /// max/2 because the most significant bit must be equal to 0 (used).
    MyIndex(size_type idx=(std::numeric_limits<size_type>::max)()/2)
      : m_idx(idx)
    {}

    /// Get the underlying index
    operator size_t() const
    { return m_idx; }

    /// reset index to be NULL
    void reset()
    { m_idx = (std::numeric_limits<size_type>::max)()/2; }

    /// return whether the handle is valid
    bool is_valid() const
    { return m_idx != (std::numeric_limits<size_type>::max)()/2; }

    // /// are two indices equal?
    // bool operator==(const Self& rhs) const
    // { return m_idx == rhs.m_idx; }

    // /// are two handles different?
    // bool operator!=(const Self& rhs) const
    // { return m_idx != rhs.m_idx; }

    // /// Comparisons
    // bool operator<(const Self& rhs) const
    // { return m_idx < rhs.m_idx; }
    // bool operator>(const Self& rhs) const
    // { return m_idx > rhs.m_idx; }
    // bool operator<=(const Self& rhs) const
    // { return m_idx <= rhs.m_idx; }
    // bool operator>=(const Self& rhs) const
    // { return m_idx >= rhs.m_idx; }

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

} // namespace internal

} //namespace CGAL

#endif // CGAL_COMPACT_CONTAINER_WITH_INDEX_H

