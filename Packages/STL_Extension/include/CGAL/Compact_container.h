// ============================================================================
//
// Copyright (c) 2003 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision: $
// release_date  : $CGAL_Date: $
//
// file          : Compact_container.h
// chapter       : $CGAL_Chapter: STL Extensions for CGAL $
// package       : $CGAL_Package: STL_Extension $
// source        : stl_extension.fw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>
//                 Lutz Kettner <kettner@mpi-sb.mpg.de>
//                 Sylvain Pion <Sylvain.Pion@mpi-sb.mpg.de>
//
// maintainer    : Michael Hoffmann <hoffmann@inf.ethz.ch>
// coordinator   : ETH
//
// A compact container.
// ============================================================================

#ifndef CGAL_COMPACT_CONTAINER_H
#define CGAL_COMPACT_CONTAINER_H 1
#include <CGAL/basic.h>

#include <iterator>
#include <algorithm>

#include <CGAL/memory.h>
#include <CGAL/iterator.h>

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
// - CC_iterator<>'s default ctor initializing to NULL ?
// - For the block_size : increment it linearly, so that for N elements,
//   we get O(sqrt(N)) blocks of O(sqrt(N)) elements each ?
//   This way we will obtain the asymptotic memory complexity we claim.
//   However, in order to allocate large blocks, we need to be able to
//   specify the block size with a member function.
// - Add .reserve() and .resize() (and proper copy of capacity_).
// - Submit Insert_iterator (and reference it in the doc).
// - Write a test and an example program, documentation in STL_extension.
// - Add preconditions in input that real pointers need to have clean bits.
//   Also for the allocated memory alignment, and sizeof().
// - Do a benchmark before/after.
// - Check the end result with Valgrind.
// - The bit squatting mecanism will be reused for the conflict flag, maybe
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

CGAL_BEGIN_NAMESPACE

// The following base class can be used to easily add a squattable pointer
// to a class (maybe you loose a bit of compactness though).
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

namespace CGALi {
  template < class DSC, class Ptr, class Ref >
  class CC_iterator;
}

template < class T, class Allocator = CGAL_ALLOCATOR(T) >
class Compact_container
{
  typedef Compact_container <T, Allocator>          Self;
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
  typedef CGALi::CC_iterator<Self, pointer, reference>
    iterator;
  typedef CGALi::CC_iterator<Self, const_pointer, const_reference>
    const_iterator;
  typedef std::reverse_iterator<iterator>           reverse_iterator;
  typedef std::reverse_iterator<const_iterator>     const_reverse_iterator;

  friend class CGALi::CC_iterator<Self, pointer, reference>;
  friend class CGALi::CC_iterator<Self, const_pointer, const_reference>;

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
  }

  iterator begin() { return iterator(first_item); }
  iterator end()   { return iterator(last_item, 0); }

#if 0
  const_iterator begin() const { return const_iterator(first_item); }
  const_iterator end()   const { return const_iterator(last_item, 0); }
#else
  // We have to cheat, because Triangulation_[23] is not const-correct.
  // The cheating should be moved to the TDS_[23] level, I guess.
  // TODO : I can commit the changes to TDS_[23] already now...
  iterator begin() const { return const_cast<Self&>(*this).begin(); }
  iterator end()   const { return const_cast<Self&>(*this).end(); }
#endif

  reverse_iterator rbegin() { return reverse_iterator(end()); }
  reverse_iterator rend()   { return reverse_iterator(begin()); }

  const_reverse_iterator
  rbegin() const { return const_reverse_iterator(end()); }
  const_reverse_iterator
  rend()   const { return const_reverse_iterator(begin()); }

  iterator insert(const T &t)
  {
    if (free_list == NULL)
      allocate_new_block();

    pointer ret = free_list;
    free_list = clean_pointee(ret);
    alloc.construct(ret, t);
    Traits::pointer(*ret) = NULL; // FIXME : assertion instead ?
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
    put_on_free_list(&*x);
    --size_;
  }

  void erase(iterator first, iterator last)
  {
    for (; first != last; ++first)
      erase(first);
  }

  void clear();

  // Merge the content of d into *this.  d gets cleared.
  // The complexity is O(size(free list = capacity-size)).
  void merge(Self &d);

  // Historical Cruft for TDS :
  // But compared to insert(), we used to avoid a copy
  // => benchmark before using it (however, TDS can now use non default ctors)
  // pointer get_new_element()
  // {
  // return &*insert(T());
  // }

  // Historical cruft for TDS.
  void release_element(pointer x)
  {
    erase(iterator(x, 0));
  }

  // Historical cruft for TDS.  It should go in TDS.
  bool is_element(const_pointer x) const
  {
    for (const_iterator it = begin(); it != end(); ++it)
      if (&*it == x)
        return true;
    return false;
  }

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

  void reserve(size_type n); // TODO

  // void resize(size_type sz, T c = T()); // TODO  makes sense ???

  bool empty() const
  {
    return size_ == 0;
  }

  allocator_type get_allocator() const
  {
    return alloc;
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
  //                          value of the last 2 bits
  // pointer part     0              1                2              3
  //         NULL     user elt       unused           free_list end  start/end
  //      != NULL     user elt       block boundary   free elt       unused
  //
  // meaning of ptr : user stuff     next/prev block  free_list      unused

  enum Type { USED = 0, BLOCK_BOUNDARY = 1, FREE = 2, START_END = 3 };

  // Using a union is clean and should avoid aliasing problems.
  union menion {
    void *         p;
    unsigned int   t:2;

    menion(void * ptr)
    : p(ptr) {}

    menion(void * ptr, Type type)
    : p(ptr)
    {
      CGAL_precondition(0 <= type && type < 4);
      t = type;
    }
  };

  // Returns the pointee, cleaned from the squatted bits (the last 2 bits).
  static pointer clean_pointee(const_pointer ptr)
  {
    return (pointer) menion(Traits::pointer(*ptr), USED).p;
  }

  // Get the type of the pointee.
  static Type type(const_pointer ptr)
  {
    return (Type) menion(Traits::pointer(*ptr)).t;
  }

  // Sets the pointer part and the type of the pointee.
  static void set_type(pointer ptr, void * p, Type t)
  {
    Traits::pointer(*ptr) = menion(p, t).p;
  }

  void init()
  {
    block_size = 14;
    capacity_  = 0;
    size_      = 0;
    free_list  = NULL;
    first_item = NULL;
    last_item  = NULL;
  }

  allocator_type   alloc;
  size_type        capacity_;
  size_type        size_;
  size_type        block_size;
  pointer          free_list;
  pointer          first_item;
  pointer          last_item;
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
  block_size = std::max(block_size, d.block_size);
  // Clear d.
  d.init();
}

template < class T, class Allocator >
void Compact_container<T, Allocator>::clear()
{
  // erase(begin(), end()); // nicer, but doesn't free memory.
  pointer p = first_item;
  while (p != NULL) { // catches the empty container case.
    ++p;
    if (type(p) == USED)
      alloc.destroy(p); // destroy used elements
    else if (type(p) == BLOCK_BOUNDARY ||
             type(p) == START_END) {
      const_pointer end = p;
      p = clean_pointee(p);
      // p becomes NULL if end of block
      alloc.deallocate(first_item, end - first_item + 1);
      capacity_ -= end - first_item -1;
      first_item = p; // keep pointer to begining of current block.
    }
  };
  CGAL_assertion(capacity_==0);
  init();
}

template < class T, class Allocator >
void Compact_container<T, Allocator>::allocate_new_block()
{
  pointer new_block = alloc.allocate(block_size + 2);
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
      set_type(last_item, NULL, START_END);
    }
  else
    {
      set_type(last_item, new_block, BLOCK_BOUNDARY);
      set_type(new_block, last_item, BLOCK_BOUNDARY);
      last_item = new_block + block_size + 1;
      set_type(last_item, NULL, START_END);
    }
    // Increase the block_size for the next time.
    block_size += 16;
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

namespace CGALi {

  template < class DSC, class Ptr, class Ref >
  class CC_iterator
  {
    typedef typename DSC::iterator                    iterator;
    typedef CC_iterator<DSC, Ptr, Ref>                Iterator;
  public:
    typedef Ptr                                       pointer;
    typedef Ref                                       reference;
    typedef typename DSC::value_type                  value_type;
    typedef typename DSC::size_type                   size_type;
    typedef typename DSC::difference_type             difference_type;
    typedef std::bidirectional_iterator_tag           iterator_category;

    CC_iterator() {}

    // Either a harmless copy-ctor,
    // or a conversion from iterator to const_iterator.
    CC_iterator(const iterator &it)
    : p(&*it) {}

  private:

    pointer p;

    // Only Compact_container should access these constructors.
    friend class Compact_container<value_type, typename DSC::allocator_type>;

  public:
    // the following should be private and explicit, but TDS needs it now...

    // For begin()
    CC_iterator(pointer ptr)
    : p(ptr)
    {
      if (p == NULL) // empty container.
        return;
      ++p; // if not empty, p = start
      if (DSC::type(p) == DSC::FREE)
        increment();
    }

  private:
    // Construction from raw pointer and for end().
    CC_iterator(pointer ptr, int)
    : p(ptr) {}

    // NB : in case empty container, begin == end == NULL.
    void increment()
    {
      // It's either pointing to end(), or valid.
      CGAL_assertion_msg(p != NULL,
                         "Doing ++ on empty container iterator ?");
      CGAL_assertion_msg(DSC::type(p) != DSC::START_END,
                         "Doing ++ on end() ?");
      // If it's not end(), then it's valid, we can do ++.
      do {
        ++p;
        if (DSC::type(p) == DSC::USED ||
            DSC::type(p) == DSC::START_END)
          return;
        if (DSC::type(p) == DSC::BLOCK_BOUNDARY)
          p = DSC::clean_pointee(p);
      }
      while (true);
    }

    void decrement()
    {
      // It's either pointing to end(), or valid.
      CGAL_assertion_msg(p != NULL,
                         "Doing -- on empty container iterator ?");
      CGAL_assertion_msg(DSC::type(p-1) != DSC::START_END,
                         "Doing -- on begin() ?");
      // If it's not begin(), then it's valid, we can do --.
      do {
        --p;
        if (DSC::type(p) == DSC::USED ||
            DSC::type(p) == DSC::START_END)
          return;
        if (DSC::type(p) == DSC::BLOCK_BOUNDARY)
          p = DSC::clean_pointee(p);
      }
      while (true);
    }

  public:

    Iterator & operator++() { increment(); return *this; }
    Iterator & operator--() { decrement(); return *this; }

    Iterator operator++(int) { Iterator tmp(*this); ++(*this); return tmp; }
    Iterator operator--(int) { Iterator tmp(*this); --(*this); return tmp; }

    reference operator*() const { return *p; }
    pointer   operator->() const { return p; }

    // Can itself be used to bit-squatting.
    void *   for_compact_container() const { return (void *) p; }
    void * & for_compact_container()       { return (void * &) p; }
  };

  // Ptr/Ref/Val could be deduced from DSC...
  template < class DSC, class Ptr, class Ref >
  inline
  bool operator==(const CC_iterator<DSC, Ptr, Ref> &rhs,
                  const CC_iterator<DSC, Ptr, Ref> &lhs)
  {
    return &*rhs == &*lhs;
  }

  template < class DSC, class Val >
  inline
  bool operator==(const CC_iterator<DSC, Val*, Val&> &rhs,
                  const CC_iterator<DSC, const Val*, const Val&> &lhs)
  {
    return &*rhs == &*lhs;
  }

  template < class DSC, class Val >
  inline
  bool operator==(const CC_iterator<DSC, const Val*, const Val&> &rhs,
                  const CC_iterator<DSC, Val*, Val&> &lhs)
  {
    return &*rhs == &*lhs;
  }

  template < class DSC, class Ptr, class Ref >
  inline
  bool operator!=(const CC_iterator<DSC, Ptr, Ref> &rhs,
                  const CC_iterator<DSC, Ptr, Ref> &lhs)
  {
    return &*rhs != &*lhs;
  }

  template < class DSC, class Val >
  inline
  bool operator!=(const CC_iterator<DSC, Val*, Val&> &rhs,
                  const CC_iterator<DSC, const Val*, const Val&> &lhs)
  {
    return &*rhs != &*lhs;
  }

  template < class DSC, class Val >
  inline
  bool operator!=(const CC_iterator<DSC, const Val*, const Val&> &rhs,
                  const CC_iterator<DSC, Val*, Val&> &lhs)
  {
    return &*rhs != &*lhs;
  }


  // The following comparison operator is here so that the iterators
  // can be put directly in set/map (i.e. std::less<> just works).
  // But maybe it's confusing ?

  template < class DSC, class Ptr, class Ref >
  inline
  bool operator<(const CC_iterator<DSC, Ptr, Ref> &rhs,
                 const CC_iterator<DSC, Ptr, Ref> &lhs)
  {
    return &*rhs < &*lhs;
  }

} // namespace CGALi

CGAL_END_NAMESPACE
#endif // CGAL_COMPACT_CONTAINER_H //
// EOF //
