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

#ifndef CGAL_COMPACT_CONTAINER_H
#define CGAL_COMPACT_CONTAINER_H

#include <CGAL/basic.h>
#include <CGAL/Default.h>

#include <iterator>
#include <algorithm>
#include <vector>
#include <cstring>

#include <CGAL/memory.h>
#include <CGAL/iterator.h>

#include <boost/mpl/if.hpp>

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

#define CGAL_INIT_COMPACT_CONTAINER_BLOCK_SIZE 14
#define CGAL_INCREMENT_COMPACT_CONTAINER_BLOCK_SIZE 16 
  
// The following base class can be used to easily add a squattable pointer
// to a class (maybe you loose a bit of compactness though).
// TODO : Shouldn't adding these bits be done automatically and transparently,
//        based on the traits class info ?
class Compact_container_base
{
  void * p;
public:
  Compact_container_base()
  : p(NULL) {}
  void *   for_compact_container() const { return p; }
  void * & for_compact_container()       { return p; }
};

// The traits class describes the way to access the pointer.
// It can be specialized.
template < class T >
struct Compact_container_traits {
  static void *   pointer(const T &t) { return t.for_compact_container(); }
  static void * & pointer(T &t)       { return t.for_compact_container(); }
};

namespace internal {
  template < class DSC, bool Const >
  class CC_iterator;
}

template < class T, class Allocator_ = Default >
class Compact_container
{
  typedef Allocator_                                Al;
  typedef typename Default::Get< Al, CGAL_ALLOCATOR(T) >::type Allocator;
  typedef Compact_container <T, Al>                 Self;
  typedef Compact_container_traits <T>              Traits;
public:
  typedef T                                         value_type;
  typedef Allocator                                 allocator_type;
  typedef typename Allocator::reference             reference;
  typedef typename Allocator::const_reference       const_reference;
  typedef typename Allocator::pointer               pointer;
  typedef typename Allocator::const_pointer         const_pointer;
  typedef typename Allocator::size_type             size_type;
  typedef typename Allocator::difference_type       difference_type;
  typedef internal::CC_iterator<Self, false>        iterator;
  typedef internal::CC_iterator<Self, true>         const_iterator;
  typedef std::reverse_iterator<iterator>           reverse_iterator;
  typedef std::reverse_iterator<const_iterator>     const_reverse_iterator;

  friend class internal::CC_iterator<Self, false>;
  friend class internal::CC_iterator<Self, true>;

  explicit Compact_container(const Allocator &a = Allocator())
  : alloc(a)
  {
    init();
  }

  template < class InputIterator >
  Compact_container(InputIterator first, InputIterator last,
                    const Allocator & a = Allocator())
  : alloc(a)
  {
    init();
    std::copy(first, last, CGAL::inserter(*this));
  }

  // The copy constructor and assignment operator preserve the iterator order
  Compact_container(const Compact_container &c)
  : alloc(c.get_allocator())
  {
    init();
    block_size = c.block_size;
    std::copy(c.begin(), c.end(), CGAL::inserter(*this));
  }

  Compact_container & operator=(const Compact_container &c)
  {
    if (&c != this) {
      Self tmp(c);
      swap(tmp);
    }
    return *this;
  }

  ~Compact_container()
  {
    clear();
  }

  void swap(Self &c)
  {
    std::swap(alloc, c.alloc);
    std::swap(capacity_, c.capacity_);
    std::swap(size_, c.size_);
    std::swap(block_size, c.block_size);
    std::swap(first_item, c.first_item);
    std::swap(last_item, c.last_item);
    std::swap(free_list, c.free_list);
    all_items.swap(c.all_items);
  }

  iterator begin() { return iterator(first_item, 0, 0); }
  iterator end()   { return iterator(last_item, 0); }

  const_iterator begin() const { return const_iterator(first_item, 0, 0); }
  const_iterator end()   const { return const_iterator(last_item, 0); }

  reverse_iterator rbegin() { return reverse_iterator(end()); }
  reverse_iterator rend()   { return reverse_iterator(begin()); }

  const_reverse_iterator
  rbegin() const { return const_reverse_iterator(end()); }
  const_reverse_iterator
  rend()   const { return const_reverse_iterator(begin()); }

  // Boost.Intrusive interface
  iterator iterator_to(reference value) const {
    return iterator(&value, 0);
  }
  const_iterator iterator_to(const_reference value) const {
    return const_iterator(&value, 0);
  }
  static iterator s_iterator_to(reference value) {
    return iterator(&value, 0);
  }
  static const_iterator s_iterator_to(const_reference value) {
    return const_iterator(&value, 0);
  }

  // Special insert methods that construct the objects in place
  // (just forward the arguments to the constructor, to optimize a copy).
#ifndef CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES
  template < typename... Args >
  iterator
  emplace(const Args&... args)
  {
    if (free_list == NULL)
      allocate_new_block();

    pointer ret = free_list;
    free_list = clean_pointee(ret);
    new (ret) value_type(args...);
    CGAL_assertion(type(ret) == USED);
    ++size_;
    return iterator(ret, 0);
  }
#else
  // inserts a default constructed item.
  iterator emplace()
  {
    if (free_list == NULL)
      allocate_new_block();

    pointer ret = free_list;
    free_list = clean_pointee(ret);
    new (ret) value_type();
    CGAL_assertion(type(ret) == USED);
    ++size_;
    return iterator(ret, 0);
  }

  template < typename T1 >
  iterator
  emplace(const T1 &t1)
  {
    if (free_list == NULL)
      allocate_new_block();

    pointer ret = free_list;
    free_list = clean_pointee(ret);
    new (ret) value_type(t1);
    CGAL_assertion(type(ret) == USED);
    ++size_;
    return iterator(ret, 0);
  }

  template < typename T1, typename T2 >
  iterator
  emplace(const T1 &t1, const T2 &t2)
  {
    if (free_list == NULL)
      allocate_new_block();

    pointer ret = free_list;
    free_list = clean_pointee(ret);
    new (ret) value_type(t1, t2);
    CGAL_assertion(type(ret) == USED);
    ++size_;
    return iterator(ret, 0);
  }

  template < typename T1, typename T2, typename T3 >
  iterator
  emplace(const T1 &t1, const T2 &t2, const T3 &t3)
  {
    if (free_list == NULL)
      allocate_new_block();

    pointer ret = free_list;
    free_list = clean_pointee(ret);
    new (ret) value_type(t1, t2, t3);
    CGAL_assertion(type(ret) == USED);
    ++size_;
    return iterator(ret, 0);
  }

  template < typename T1, typename T2, typename T3, typename T4 >
  iterator
  emplace(const T1 &t1, const T2 &t2, const T3 &t3, const T4 &t4)
  {
    if (free_list == NULL)
      allocate_new_block();

    pointer ret = free_list;
    free_list = clean_pointee(ret);
    new (ret) value_type(t1, t2, t3, t4);
    CGAL_assertion(type(ret) == USED);
    ++size_;
    return iterator(ret, 0);
  }

  template < typename T1, typename T2, typename T3, typename T4, typename T5 >
  iterator
  emplace(const T1 &t1, const T2 &t2, const T3 &t3, const T4 &t4,
	  const T5 &t5)
  {
    if (free_list == NULL)
      allocate_new_block();

    pointer ret = free_list;
    free_list = clean_pointee(ret);
    new (ret) value_type(t1, t2, t3, t4, t5);
    CGAL_assertion(type(ret) == USED);
    ++size_;
    return iterator(ret, 0);
  }

  template < typename T1, typename T2, typename T3, typename T4,
             typename T5, typename T6 >
  iterator
  emplace(const T1 &t1, const T2 &t2, const T3 &t3, const T4 &t4,
          const T5 &t5, const T6 &t6)
  {
    if (free_list == NULL)
      allocate_new_block();

    pointer ret = free_list;
    free_list = clean_pointee(ret);
    new (ret) value_type(t1, t2, t3, t4, t5, t6);
    CGAL_assertion(type(ret) == USED);
    ++size_;
    return iterator(ret, 0);
  }

  template < typename T1, typename T2, typename T3, typename T4,
             typename T5, typename T6, typename T7 >
  iterator
  emplace(const T1 &t1, const T2 &t2, const T3 &t3, const T4 &t4,
          const T5 &t5, const T6 &t6, const T7 &t7)
  {
    if (free_list == NULL)
      allocate_new_block();

    pointer ret = free_list;
    free_list = clean_pointee(ret);
    new (ret) value_type(t1, t2, t3, t4, t5, t6, t7);
    CGAL_assertion(type(ret) == USED);
    ++size_;
    return iterator(ret, 0);
  }

  template < typename T1, typename T2, typename T3, typename T4,
             typename T5, typename T6, typename T7, typename T8 >
  iterator
  emplace(const T1 &t1, const T2 &t2, const T3 &t3, const T4 &t4,
          const T5 &t5, const T6 &t6, const T7 &t7, const T8 &t8)
  {
    if (free_list == NULL)
      allocate_new_block();

    pointer ret = free_list;
    free_list = clean_pointee(ret);
    new (ret) value_type(t1, t2, t3, t4, t5, t6, t7, t8);
    CGAL_assertion(type(ret) == USED);
    ++size_;
    return iterator(ret, 0);
  }
#endif // CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES

  iterator insert(const T &t)
  {
    if (free_list == NULL)
      allocate_new_block();

    pointer ret = free_list;
    free_list = clean_pointee(ret);
    alloc.construct(ret, t);
    CGAL_assertion(type(ret) == USED);
    ++size_;
    return iterator(ret, 0);
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

  void erase(iterator x)
  {
    CGAL_precondition(type(&*x) == USED);
    alloc.destroy(&*x);
#ifndef CGAL_NO_ASSERTIONS
    std::memset(&*x, 0, sizeof(T));
#endif
    put_on_free_list(&*x);
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

    for (typename All_items::const_iterator it = all_items.begin(), itend = all_items.end();
         it != itend; ++it) {
      const_pointer p = it->first;
      size_type s = it->second;

      // Are we in the address range of this block (excluding first and last
      // elements) ?
      if (c <= p || (p+s-1) <= c)
        continue;

      CGAL_assertion_msg( (c-p)+p == c, "wrong alignment of iterator");

      return type(c) == USED;
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
    size_type tmp = block_size;
    block_size = (std::max)( n - capacity_, block_size );
    allocate_new_block();
    block_size = tmp+CGAL_INCREMENT_COMPACT_CONTAINER_BLOCK_SIZE;
  }
  
private:

  void allocate_new_block();

  void put_on_free_list(pointer x)
  {
    set_type(x, free_list, FREE);
    free_list = x;
  }

  // Definition of the bit squatting :
  // =================================
  // ptr is composed of a pointer part and the last 2 bits.
  // Here is the meaning of each of the 8 cases.
  //
  //                          value of the last 2 bits as "Type"
  // pointer part     0              1                2              3
  //         NULL     user elt       unused           free_list end  start/end
  //      != NULL     user elt       block boundary   free elt       unused
  //
  // meaning of ptr : user stuff     next/prev block  free_list      unused

  enum Type { USED = 0, BLOCK_BOUNDARY = 1, FREE = 2, START_END = 3 };

  // The bit squatting is implemented by casting pointers to (char *), then
  // subtracting to NULL, doing bit manipulations on the resulting integer,
  // and converting back.

  static char * clean_pointer(char * p)
  {
    return ((p - (char *) NULL) & ~ (std::ptrdiff_t) START_END) + (char *) NULL;
  }

  // Returns the pointee, cleaned up from the squatted bits.
  static pointer clean_pointee(const_pointer ptr)
  {
    return (pointer) clean_pointer((char *) Traits::pointer(*ptr));
  }

  // Get the type of the pointee.
  static Type type(const_pointer ptr)
  {
    char * p = (char *) Traits::pointer(*ptr);
    return (Type) (p - clean_pointer(p));
  }

  // Sets the pointer part and the type of the pointee.
  static void set_type(pointer ptr, void * p, Type t)
  {
    // This out of range compare is always true and causes lots of
    // unnecessary warnings.
    // CGAL_precondition(0 <= t && t < 4); 
    Traits::pointer(*ptr) = (void *) ((clean_pointer((char *) p)) + (int) t);
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
    block_size = CGAL_INIT_COMPACT_CONTAINER_BLOCK_SIZE;
    capacity_  = 0;
    size_      = 0;
    free_list  = NULL;
    first_item = NULL;
    last_item  = NULL;
    all_items  = All_items();
  }

  allocator_type   alloc;
  size_type        capacity_;
  size_type        size_;
  size_type        block_size;
  pointer          free_list;
  pointer          first_item;
  pointer          last_item;
  All_items        all_items;
};

template < class T, class Allocator >
void Compact_container<T, Allocator>::merge(Self &d)
{
  CGAL_precondition(&d != this);

  // Allocators must be "compatible" :
  CGAL_precondition(get_allocator() == d.get_allocator());

  // Concatenate the free_lists.
  if (free_list == NULL) {
    free_list = d.free_list;
  } else if (d.free_list != NULL) {
    pointer p = free_list;
    while (clean_pointee(p) != NULL)
      p = clean_pointee(p);
    set_type(p, d.free_list, FREE);
  }
  // Concatenate the blocks.
  if (last_item == NULL) { // empty...
    first_item = d.first_item;
    last_item  = d.last_item;
  } else if (d.last_item != NULL) {
    set_type(last_item, d.first_item, BLOCK_BOUNDARY);
    set_type(d.first_item, last_item, BLOCK_BOUNDARY);
    last_item = d.last_item;
  }
  // Add the sizes.
  size_ += d.size_;
  // Add the capacities.
  capacity_ += d.capacity_;
  // It seems reasonnable to take the max of the block sizes.
  block_size = (std::max)(block_size, d.block_size);
  // Clear d.
  d.init();
}

template < class T, class Allocator >
void Compact_container<T, Allocator>::clear()
{
  for (typename All_items::iterator it = all_items.begin(), itend = all_items.end();
       it != itend; ++it) {
    pointer p = it->first;
    size_type s = it->second;
    for (pointer pp = p + 1; pp != p + s - 1; ++pp) {
      if (type(pp) == USED)
        alloc.destroy(pp);
    }
    alloc.deallocate(p, s);
  }
  init();
}

template < class T, class Allocator >
void Compact_container<T, Allocator>::allocate_new_block()
{
  pointer new_block = alloc.allocate(block_size + 2);
  all_items.push_back(std::make_pair(new_block, block_size + 2));
  capacity_ += block_size;
  // We don't touch the first and the last one.
  // We mark them free in reverse order, so that the insertion order
  // will correspond to the iterator order...
  for (size_type i = block_size; i >= 1; --i)
    put_on_free_list(new_block + i);
  // We insert this new block at the end.
  if (last_item == NULL) // First time
  {
      first_item = new_block;
      last_item  = new_block + block_size + 1;
      set_type(first_item, NULL, START_END);
  }
  else
  {
      set_type(last_item, new_block, BLOCK_BOUNDARY);
      set_type(new_block, last_item, BLOCK_BOUNDARY);
      last_item = new_block + block_size + 1;
  }
  set_type(last_item, NULL, START_END);
  // Increase the block_size for the next time.
  block_size += CGAL_INCREMENT_COMPACT_CONTAINER_BLOCK_SIZE;
}

template < class T, class Allocator >
inline
bool operator==(const Compact_container<T, Allocator> &lhs,
                const Compact_container<T, Allocator> &rhs)
{
  return lhs.size() == rhs.size() &&
    std::equal(lhs.begin(), lhs.end(), rhs.begin());
}

template < class T, class Allocator >
inline
bool operator!=(const Compact_container<T, Allocator> &lhs,
                const Compact_container<T, Allocator> &rhs)
{
  return ! (lhs == rhs);
}

template < class T, class Allocator >
inline
bool operator< (const Compact_container<T, Allocator> &lhs,
                const Compact_container<T, Allocator> &rhs)
{
  return std::lexicographical_compare(lhs.begin(), lhs.end(),
                                      rhs.begin(), rhs.end());
}

template < class T, class Allocator >
inline
bool operator> (const Compact_container<T, Allocator> &lhs,
                const Compact_container<T, Allocator> &rhs)
{
  return rhs < lhs;
}

template < class T, class Allocator >
inline
bool operator<=(const Compact_container<T, Allocator> &lhs,
                const Compact_container<T, Allocator> &rhs)
{
  return ! (lhs > rhs);
}

template < class T, class Allocator >
inline
bool operator>=(const Compact_container<T, Allocator> &lhs,
                const Compact_container<T, Allocator> &rhs)
{
  return ! (lhs < rhs);
}

namespace internal {

  template < class DSC, bool Const >
  class CC_iterator
  {
    typedef typename DSC::iterator                    iterator;
    typedef CC_iterator<DSC, Const>                   Self;
  public:
    typedef typename DSC::value_type                  value_type;
    typedef typename DSC::size_type                   size_type;
    typedef typename DSC::difference_type             difference_type;
    typedef typename boost::mpl::if_c< Const, const value_type*,
                                       value_type*>::type pointer;
    typedef typename boost::mpl::if_c< Const, const value_type&,
                                       value_type&>::type reference;
    typedef std::bidirectional_iterator_tag           iterator_category;

    // the initialization with NULL is required by our Handle concept.
    CC_iterator()
    {
      m_ptr.p = NULL;
    }

    // Either a harmless copy-ctor,
    // or a conversion from iterator to const_iterator.
    CC_iterator (const iterator &it)
    {
      m_ptr.p = &(*it);
    }

    // Same for assignment operator (otherwise MipsPro warns)
    CC_iterator & operator= (const iterator &it)
    {
      m_ptr.p = &(*it);
      return *this;
    }

    // Construction from NULL
    CC_iterator (Nullptr_t CGAL_assertion_code(n))
    {
      CGAL_assertion (n == NULL);
      m_ptr.p = NULL;
    }

  private:

    union {
      pointer      p;
      void        *vp;
    } m_ptr;

    // Only Compact_container should access these constructors.
    friend class Compact_container<value_type, typename DSC::Al>;

    // For begin()
    CC_iterator(pointer ptr, int, int)
    {
      m_ptr.p = ptr;
      if (m_ptr.p == NULL) // empty container.
        return;

      ++(m_ptr.p); // if not empty, p = start
      if (DSC::type(m_ptr.p) == DSC::FREE)
        increment();
    }

    // Construction from raw pointer and for end().
    CC_iterator(pointer ptr, int)
    {
      m_ptr.p = ptr;
    }

    // NB : in case empty container, begin == end == NULL.
    void increment()
    {
      // It's either pointing to end(), or valid.
      CGAL_assertion_msg(m_ptr.p != NULL,
	 "Incrementing a singular iterator or an empty container iterator ?");
      CGAL_assertion_msg(DSC::type(m_ptr.p) != DSC::START_END,
	 "Incrementing end() ?");

      // If it's not end(), then it's valid, we can do ++.
      do {
        ++(m_ptr.p);
        if (DSC::type(m_ptr.p) == DSC::USED ||
            DSC::type(m_ptr.p) == DSC::START_END)
          return;

        if (DSC::type(m_ptr.p) == DSC::BLOCK_BOUNDARY)
          m_ptr.p = DSC::clean_pointee(m_ptr.p);
      } while (true);
    }

    void decrement()
    {
      // It's either pointing to end(), or valid.
      CGAL_assertion_msg(m_ptr.p != NULL,
	 "Decrementing a singular iterator or an empty container iterator ?");
      CGAL_assertion_msg(DSC::type(m_ptr.p - 1) != DSC::START_END,
	 "Decrementing begin() ?");

      // If it's not begin(), then it's valid, we can do --.
      do {
        --m_ptr.p;
        if (DSC::type(m_ptr.p) == DSC::USED ||
            DSC::type(m_ptr.p) == DSC::START_END)
          return;

        if (DSC::type(m_ptr.p) == DSC::BLOCK_BOUNDARY)
          m_ptr.p = DSC::clean_pointee(m_ptr.p);
      } while (true);
    }

  public:

    Self & operator++()
    {
      CGAL_assertion_msg(m_ptr.p != NULL,
	 "Incrementing a singular iterator or an empty container iterator ?");
      CGAL_assertion_msg(DSC::type(m_ptr.p) == DSC::USED,
                         "Incrementing an invalid iterator.");
      increment();
      return *this;
    }

    Self & operator--()
    {
      CGAL_assertion_msg(m_ptr.p != NULL,
	 "Decrementing a singular iterator or an empty container iterator ?");
      CGAL_assertion_msg(DSC::type(m_ptr.p) == DSC::USED
		      || DSC::type(m_ptr.p) == DSC::START_END,
                         "Decrementing an invalid iterator.");
      decrement();
      return *this;
    }

    Self operator++(int) { Self tmp(*this); ++(*this); return tmp; }
    Self operator--(int) { Self tmp(*this); --(*this); return tmp; }

    reference operator*() const { return *(m_ptr.p); }

    pointer   operator->() const { return (m_ptr.p); }

    // For std::less...
    bool operator<(const CC_iterator& other) const
    {
      return (m_ptr.p < other.m_ptr.p);
    }

    bool operator>(const CC_iterator& other) const
    {
      return (m_ptr.p > other.m_ptr.p);
    }

    bool operator<=(const CC_iterator& other) const
    {
      return (m_ptr.p <= other.m_ptr.p);
    }

    bool operator>=(const CC_iterator& other) const
    {
      return (m_ptr.p >= other.m_ptr.p);
    }

    // Can itself be used for bit-squatting.
    void *   for_compact_container() const { return (m_ptr.vp); }
    void * & for_compact_container()       { return (m_ptr.vp); }
  };

  template < class DSC, bool Const1, bool Const2 >
  inline
  bool operator==(const CC_iterator<DSC, Const1> &rhs,
                  const CC_iterator<DSC, Const2> &lhs)
  {
    return &*rhs == &*lhs;
  }

  template < class DSC, bool Const1, bool Const2 >
  inline
  bool operator!=(const CC_iterator<DSC, Const1> &rhs,
                  const CC_iterator<DSC, Const2> &lhs)
  {
    return &*rhs != &*lhs;
  }

  // Comparisons with NULL are part of CGAL's Handle concept...
  template < class DSC, bool Const >
  inline
  bool operator==(const CC_iterator<DSC, Const> &rhs,
                  Nullptr_t CGAL_assertion_code(n))
  {
    CGAL_assertion( n == NULL);
    return &*rhs == NULL;
  }

  template < class DSC, bool Const >
  inline
  bool operator!=(const CC_iterator<DSC, Const> &rhs,
		  Nullptr_t CGAL_assertion_code(n))
  {
    CGAL_assertion( n == NULL);
    return &*rhs != NULL;
  }

} // namespace internal

} //namespace CGAL

#endif // CGAL_COMPACT_CONTAINER_H
