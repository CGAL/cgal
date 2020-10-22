// Copyright (c) 2003,2004,2007-2010  INRIA Sophia-Antipolis (France).
// Copyright (c) 2014  GeometryFactory Sarl (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Sylvain Pion

#ifndef CGAL_COMPACT_CONTAINER_H
#define CGAL_COMPACT_CONTAINER_H

#include <CGAL/disable_warnings.h>

#include <CGAL/config.h>
#include <CGAL/Default.h>

#include <cmath>
#include <cstddef>
#include <iterator>
#include <algorithm>
#include <vector>
#include <cstring>
#include <functional>
#include <atomic>

#include <CGAL/memory.h>
#include <CGAL/iterator.h>
#include <CGAL/CC_safe_handle.h>
#include <CGAL/Time_stamper.h>
#include <CGAL/Has_member.h>

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

template<unsigned int first_block_size_, unsigned int block_size_increment>
struct Addition_size_policy
{
  static const unsigned int first_block_size = first_block_size_;

  template<typename Compact_container>
  static void increase_size(Compact_container& cc)
  { cc.block_size += block_size_increment; }

  template<typename Compact_container>
  static void get_index_and_block(typename Compact_container::size_type i,
                                  typename Compact_container::size_type& index,
                                  typename Compact_container::size_type& block)
  {
    typedef typename Compact_container::size_type ST;
    const ST TWO_M_N = 2*first_block_size_ - block_size_increment;
    ST delta = TWO_M_N*TWO_M_N + 8*block_size_increment*i;
    block= (static_cast<ST>(std::sqrt(static_cast<double>(delta))) - TWO_M_N)
      / (2*block_size_increment);

    if ( block==0 )
    { index = i + 1; }
    else
    {
      const typename Compact_container::size_type first_element_in_block =
        block*(first_block_size_+ (block_size_increment*(block - 1))/2);

      index=i - first_element_in_block + 1;
    }
  }
};

template<unsigned int k>
struct Constant_size_policy
{
  static const unsigned int first_block_size = k;

  template<typename Compact_container>
  static void increase_size(Compact_container& /*cc*/)
  {}

  template<typename Compact_container>
  static void get_index_and_block(typename Compact_container::size_type i,
                                  typename Compact_container::size_type& index,
                                  typename Compact_container::size_type& block)
  {
    block=i/k;
    index=(i%k)+1;
  }
};

// The following base class can be used to easily add a squattable pointer
// to a class (maybe you lose a bit of compactness though).
// TODO : Shouldn't adding these bits be done automatically and transparently,
//        based on the traits class info ?
class Compact_container_base
{
  void * p;
public:
  Compact_container_base()
  : p(nullptr) {}
  void *   for_compact_container() const { return p; }
  void for_compact_container(void* ptr)  { p = ptr; }
};

// The traits class describes the way to access the pointer.
// It can be specialized.
template < class T >
struct Compact_container_traits {
  static void *   pointer(const T &t)    { return t.for_compact_container(); }
  static void set_pointer(T &t, void* p) { t.for_compact_container(p); }
};

namespace internal {
  template < class DSC, bool Const >
  class CC_iterator;

  CGAL_GENERATE_MEMBER_DETECTOR(increment_erase_counter);

  // A basic "no erase counter" strategy
  template <bool Has_erase_counter_tag>
  class Erase_counter_strategy {
  public:
    // Do nothing
    template <typename Element>
    static unsigned int erase_counter(const Element &) { return 0; }
    template <typename Element>
    static void set_erase_counter(Element &, unsigned int) {}
    template <typename Element>
    static void increment_erase_counter(Element &) {}
  };


  // A strategy managing an internal counter
  template <>
  class Erase_counter_strategy<true>
  {
  public:
    template <typename Element>
    static unsigned int erase_counter(const Element &e)
    {
      return e.erase_counter();
    }

    template <typename Element>
    static void set_erase_counter(Element &e, unsigned int c)
    {
      e.set_erase_counter(c);
    }

    template <typename Element>
    static void increment_erase_counter(Element &e)
    {
      e.increment_erase_counter();
    }
  };
}

template < class T,
           class Allocator_ = Default,
           class Increment_policy_ = Default,
           class TimeStamper_ = Default >
class Compact_container
{
  typedef Allocator_                                Al;
  typedef typename Default::Get< Al, CGAL_ALLOCATOR(T) >::type Allocator;
  typedef Increment_policy_                         Ip;
  typedef typename Default::Get< Ip,
            Addition_size_policy<CGAL_INIT_COMPACT_CONTAINER_BLOCK_SIZE,
                             CGAL_INCREMENT_COMPACT_CONTAINER_BLOCK_SIZE>
          >::type                                   Increment_policy;
  typedef TimeStamper_                              Ts;
  typedef Compact_container <T, Al, Ip, Ts>         Self;
  typedef Compact_container_traits <T>              Traits;
public:
  typedef typename Default::Get< TimeStamper_,
                                 CGAL::Time_stamper_impl<T> >::type
                                                    Time_stamper;
  typedef Time_stamper                              Time_stamper_impl; // backward-compatibility

  typedef T                                         value_type;
  typedef Allocator                                 allocator_type;

  typedef value_type&                               reference;
  typedef const value_type&                         const_reference;

  typedef typename std::allocator_traits<Allocator>::pointer               pointer;
  typedef typename std::allocator_traits<Allocator>::const_pointer         const_pointer;
  typedef typename std::allocator_traits<Allocator>::size_type             size_type;
  typedef typename std::allocator_traits<Allocator>::difference_type       difference_type;

  typedef internal::CC_iterator<Self, false>        iterator;
  typedef internal::CC_iterator<Self, true>         const_iterator;
  typedef std::reverse_iterator<iterator>           reverse_iterator;
  typedef std::reverse_iterator<const_iterator>     const_reverse_iterator;

  friend class internal::CC_iterator<Self, false>;
  friend class internal::CC_iterator<Self, true>;

  template<unsigned int first_block_size_, unsigned int block_size_increment>
    friend struct Addition_size_policy;
  template<unsigned int k> friend struct Constant_size_policy;

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
    time_stamp = c.time_stamp.load();
    std::copy(c.begin(), c.end(), CGAL::inserter(*this));
  }

  Compact_container(Compact_container&& c) noexcept
  : alloc(c.get_allocator())
  {
    c.swap(*this);
  }

  Compact_container & operator=(const Compact_container &c)
  {
    if (&c != this) {
      Self tmp(c);
      swap(tmp);
    }
    return *this;
  }

  Compact_container & operator=(Compact_container&& c) noexcept
  {
    Self tmp(std::move(c));
    tmp.swap(*this);
    return *this;
  }

  ~Compact_container()
  {
    clear();
  }

  bool is_used(const_iterator ptr) const
  {
    return (type(&*ptr)==USED);
  }

  bool is_used(size_type i) const
  {
    typename Self::size_type block_number, index_in_block;
    Increment_policy::template get_index_and_block<Self>(i,
                                                         index_in_block,
                                                         block_number);
    return (type(&all_items[block_number].first[index_in_block])
                 == USED);
  }

  const T& operator[] (size_type i) const
  {
    CGAL_assertion( is_used(i) );

    typename Self::size_type block_number, index_in_block;
    Increment_policy::template get_index_and_block<Self>(i,
                                                         index_in_block,
                                                         block_number);
    return all_items[block_number].first[index_in_block];
  }

  T& operator[] (size_type i)
  {
    CGAL_assertion( is_used(i) );

    typename Self::size_type block_number, index_in_block;
    Increment_policy::template get_index_and_block<Self>(i,
                                                         index_in_block,
                                                         block_number);
    return all_items[block_number].first[index_in_block];
  }

  friend void swap(Compact_container& a, Compact_container b) {
    a.swap(b);
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
  template < typename... Args >
  iterator
  emplace(const Args&... args)
  {
    if (free_list == nullptr)
      allocate_new_block();

    pointer ret = free_list;
    free_list = clean_pointee(ret);
    new (ret) value_type(args...);
    CGAL_assertion(type(ret) == USED);
    ++size_;
    Time_stamper::set_time_stamp(ret, time_stamp);
    return iterator(ret, 0);
  }

  iterator insert(const T &t)
  {
    if (free_list == nullptr)
      allocate_new_block();

    pointer ret = free_list;
    free_list = clean_pointee(ret);
    std::allocator_traits<allocator_type>::construct(alloc, ret, t);
    CGAL_assertion(type(ret) == USED);
    ++size_;
    Time_stamper::set_time_stamp(ret, time_stamp);
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
    typedef internal::Erase_counter_strategy<
      internal::has_increment_erase_counter<T>::value> EraseCounterStrategy;

    CGAL_precondition(type(&*x) == USED);
    EraseCounterStrategy::increment_erase_counter(*x);
    std::allocator_traits<allocator_type>::destroy(alloc, &*x);
/*#ifndef CGAL_NO_ASSERTIONS
    std::memset(&*x, 0, sizeof(T));
#endif*/
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

  // Returns the index of the iterator "cit", i.e. the number n so that
  // operator[](n)==*cit.
  // Complexity : O(#blocks) = O(sqrt(capacity())).
  // This function is mostly useful for purposes of efficient debugging at
  // higher levels.
  size_type index(const_iterator cit) const
  {
    // We use the block structure to provide an efficient version :
    // we check if the address is in the range of each block.

    assert(cit != end());

    const_pointer c = &*cit;
    size_type res=0;

    for (typename All_items::const_iterator it = all_items.begin(), itend = all_items.end();
         it != itend; ++it) {
      const_pointer p = it->first;
      size_type s = it->second;

      // Are we in the address range of this block (excluding first and last
      // elements) ?
      if ( p<c && c<(p+s-1) )
      {
        CGAL_assertion_msg( (c-p)+p == c, "wrong alignment of iterator");
        return res+(c-p-1);
      }

      res += s-2;
    }

    return (size_type)-1; // cit does not belong to this compact container
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

    size_type lastblock = all_items.size();

    while ( capacity_<n )
    { // Pb because the order of free list is no more the order of
      // allocate_new_block();
      pointer new_block = alloc.allocate(block_size + 2);
      all_items.push_back(std::make_pair(new_block, block_size + 2));
      capacity_ += block_size;
      // We insert this new block at the end.
      if (last_item == nullptr) // First time
      {
        first_item = new_block;
        last_item  = new_block + block_size + 1;
        set_type(first_item, nullptr, START_END);
      }
      else
      {
        set_type(last_item, new_block, BLOCK_BOUNDARY);
        set_type(new_block, last_item, BLOCK_BOUNDARY);
        last_item = new_block + block_size + 1;
      }
      set_type(last_item, nullptr, START_END);
      // Increase the block_size for the next time.
      Increment_policy::increase_size(*this);
    }

    // Now we put all the new elements on freelist, starting from the last block
    // inserted and mark them free in reverse order, so that the insertion order
    // will correspond to the iterator order...
    // We don't touch the first and the last one.
    size_type curblock=all_items.size();
    do
    {
      --curblock; // We are sure we have at least create a new block
      pointer new_block = all_items[curblock].first;
      for (size_type i = all_items[curblock].second-2; i >= 1; --i)
        put_on_free_list(new_block + i);
    }
    while ( curblock>lastblock );
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
  //         nullptr     user elt       unused           free_list end  start/end
  //      != nullptr     user elt       block boundary   free elt       unused
  //
  // meaning of ptr : user stuff     next/prev block  free_list      unused

  enum Type { USED = 0, BLOCK_BOUNDARY = 1, FREE = 2, START_END = 3 };

  // The bit squatting is implemented by casting pointers to (char *), then
  // subtracting to nullptr, doing bit manipulations on the resulting integer,
  // and converting back.

  static char * clean_pointer(char * p)
  {
    return reinterpret_cast<char*>(reinterpret_cast<std::ptrdiff_t>(p) &
                                   ~ (std::ptrdiff_t) START_END);
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
    return (Type) (reinterpret_cast<std::ptrdiff_t>(p) -
                   reinterpret_cast<std::ptrdiff_t>(clean_pointer(p)));
  }

  // Sets the pointer part and the type of the pointee.
  static void set_type(pointer ptr, void * p, Type t)
  {
    // This out of range compare is always true and causes lots of
    // unnecessary warnings.
    // CGAL_precondition(0 <= t && t < 4);
    Traits::set_pointer(*ptr, reinterpret_cast<void *>
      (reinterpret_cast<std::ptrdiff_t>(clean_pointer((char *) p)) + (int) t));
  }

public:
  // @return true iff pts is on the beginning or on the end of its block.
  static bool is_begin_or_end(const_pointer ptr)
  { return type(ptr)==START_END; }

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

    // non-atomic swap of time_stamp:
    c.time_stamp = time_stamp.exchange(c.time_stamp.load());
  }
private:
  // We store a vector of pointers to all allocated blocks and their sizes.
  // Knowing all pointers, we don't have to walk to the end of a block to reach
  // the pointer to the next block.
  // Knowing the sizes allows to deallocate() without having to compute the size
  // by walking through the block till its end.
  // This opens up the possibility for the compiler to optimize the clear()
  // function considerably when has_trivial_destructor<T>.
  using All_items = std::vector<std::pair<pointer, size_type> >;

  using time_stamp_t = std::atomic<std::size_t>;

  void init()
  {
    block_size = Increment_policy::first_block_size;
    capacity_  = 0;
    size_      = 0;
    free_list  = nullptr;
    first_item = nullptr;
    last_item  = nullptr;
    all_items  = All_items();
    time_stamp = 0;
  }

  allocator_type   alloc;
  size_type        capacity_   = 0;
  size_type        size_       = 0;
  size_type        block_size  = Increment_policy::first_block_size;
  pointer          free_list   = nullptr;
  pointer          first_item  = nullptr;
  pointer          last_item   = nullptr;
  All_items        all_items   = {};
  time_stamp_t     time_stamp  = {};
};

template < class T, class Allocator, class Increment_policy, class TimeStamper >
void Compact_container<T, Allocator, Increment_policy, TimeStamper>::merge(Self &d)
{
  CGAL_precondition(&d != this);

  // Allocators must be "compatible" :
  CGAL_precondition(get_allocator() == d.get_allocator());

  // Concatenate the free_lists.
  if (free_list == nullptr) {
    free_list = d.free_list;
  } else if (d.free_list != nullptr) {
    pointer p = free_list;
    while (clean_pointee(p) != nullptr)
      p = clean_pointee(p);
    set_type(p, d.free_list, FREE);
  }
  // Concatenate the blocks.
  if (last_item == nullptr) { // empty...
    first_item = d.first_item;
    last_item  = d.last_item;
  } else if (d.last_item != nullptr) {
    set_type(last_item, d.first_item, BLOCK_BOUNDARY);
    set_type(d.first_item, last_item, BLOCK_BOUNDARY);
    last_item = d.last_item;
  }
  all_items.insert(all_items.end(), d.all_items.begin(), d.all_items.end());
  // Add the sizes.
  size_ += d.size_;
  // Add the capacities.
  capacity_ += d.capacity_;
  // It seems reasonnable to take the max of the block sizes.
  block_size = (std::max)(block_size, d.block_size);
  // Clear d.
  d.init();
}

template < class T, class Allocator, class Increment_policy, class TimeStamper >
void Compact_container<T, Allocator, Increment_policy, TimeStamper>::clear()
{
  for (typename All_items::iterator it = all_items.begin(), itend = all_items.end();
       it != itend; ++it) {
    pointer p = it->first;
    size_type s = it->second;
    for (pointer pp = p + 1; pp != p + s - 1; ++pp) {
      if (type(pp) == USED)
      {
        std::allocator_traits<allocator_type>::destroy(alloc, pp);
        set_type(pp, nullptr, FREE);
      }
    }
    alloc.deallocate(p, s);
  }
  init();
}

template < class T, class Allocator, class Increment_policy, class TimeStamper >
void Compact_container<T, Allocator, Increment_policy, TimeStamper>::allocate_new_block()
{
  typedef internal::Erase_counter_strategy<
    internal::has_increment_erase_counter<T>::value> EraseCounterStrategy;

  pointer new_block = alloc.allocate(block_size + 2);
  all_items.push_back(std::make_pair(new_block, block_size + 2));
  capacity_ += block_size;
  // We don't touch the first and the last one.
  // We mark them free in reverse order, so that the insertion order
  // will correspond to the iterator order...
  for (size_type i = block_size; i >= 1; --i)
  {
    EraseCounterStrategy::set_erase_counter(*(new_block + i), 0);
    Time_stamper::initialize_time_stamp(new_block + i);
    put_on_free_list(new_block + i);
  }
  // We insert this new block at the end.
  if (last_item == nullptr) // First time
  {
      first_item = new_block;
      last_item  = new_block + block_size + 1;
      set_type(first_item, nullptr, START_END);
  }
  else
  {
      set_type(last_item, new_block, BLOCK_BOUNDARY);
      set_type(new_block, last_item, BLOCK_BOUNDARY);
      last_item = new_block + block_size + 1;
  }
  set_type(last_item, nullptr, START_END);
  // Increase the block_size for the next time.
  Increment_policy::increase_size(*this);
}

template < class T, class Allocator, class Increment_policy, class TimeStamper >
inline
bool operator==(const Compact_container<T, Allocator, Increment_policy, TimeStamper> &lhs,
                const Compact_container<T, Allocator, Increment_policy, TimeStamper> &rhs)
{
  return lhs.size() == rhs.size() &&
    std::equal(lhs.begin(), lhs.end(), rhs.begin());
}

template < class T, class Allocator, class Increment_policy, class TimeStamper >
inline
bool operator!=(const Compact_container<T, Allocator, Increment_policy, TimeStamper> &lhs,
                const Compact_container<T, Allocator, Increment_policy, TimeStamper> &rhs)
{
  return ! (lhs == rhs);
}

template < class T, class Allocator, class Increment_policy, class TimeStamper >
inline
bool operator< (const Compact_container<T, Allocator, Increment_policy, TimeStamper> &lhs,
                const Compact_container<T, Allocator, Increment_policy, TimeStamper> &rhs)
{
  return std::lexicographical_compare(lhs.begin(), lhs.end(),
                                      rhs.begin(), rhs.end());
}

template < class T, class Allocator, class Increment_policy, class TimeStamper >
inline
bool operator> (const Compact_container<T, Allocator, Increment_policy, TimeStamper> &lhs,
                const Compact_container<T, Allocator, Increment_policy, TimeStamper> &rhs)
{
  return rhs < lhs;
}

template < class T, class Allocator, class Increment_policy, class TimeStamper >
inline
bool operator<=(const Compact_container<T, Allocator, Increment_policy, TimeStamper> &lhs,
                const Compact_container<T, Allocator, Increment_policy, TimeStamper> &rhs)
{
  return ! (lhs > rhs);
}

template < class T, class Allocator, class Increment_policy, class TimeStamper >
inline
bool operator>=(const Compact_container<T, Allocator, Increment_policy, TimeStamper> &lhs,
                const Compact_container<T, Allocator, Increment_policy, TimeStamper> &rhs)
{
  return ! (lhs < rhs);
}

// forward-declare Concurrent_compact_container, for CC_iterator
template < class T, class Allocator_ >
class Concurrent_compact_container;

namespace internal {

  template < class DSC, bool Const >
  class CC_iterator
  {
    typedef CC_iterator<DSC, Const>                   Self;
  public:
    typedef DSC                                       CC;
    typedef typename DSC::value_type                  value_type;
    typedef typename DSC::size_type                   size_type;
    typedef typename DSC::difference_type             difference_type;
    typedef typename boost::mpl::if_c< Const, const value_type*,
                                       value_type*>::type pointer;
    typedef typename boost::mpl::if_c< Const, const value_type&,
                                       value_type&>::type reference;
    typedef std::bidirectional_iterator_tag           iterator_category;

    // the initialization with nullptr is required by our Handle concept.
    CC_iterator()
#ifdef CGAL_COMPACT_CONTAINER_DEBUG_TIME_STAMP
      : ts(0)
#endif
    {
      m_ptr = nullptr;
    }

    // Converting constructor from mutable to constant iterator
    template <bool OtherConst>
    CC_iterator(const CC_iterator<
                typename std::enable_if<(!OtherConst && Const), DSC>::type,
                OtherConst> &const_it)
#ifdef CGAL_COMPACT_CONTAINER_DEBUG_TIME_STAMP
        : ts(Time_stamper::time_stamp(const_it.operator->()))
#endif
    {
      m_ptr = const_it.operator->();
    }

    // Assignment operator from mutable to constant iterator
    template <bool OtherConst>
    CC_iterator & operator= (const CC_iterator<
                typename std::enable_if<(!OtherConst && Const), DSC>::type,
                OtherConst> &const_it)
    {
      m_ptr = const_it.operator->();
#ifdef CGAL_COMPACT_CONTAINER_DEBUG_TIME_STAMP
      ts = Time_stamper::time_stamp(const_it.operator->());
#endif
      return *this;
    }

    // Construction from nullptr
    CC_iterator (std::nullptr_t /*CGAL_assertion_code(n)*/)
#ifdef CGAL_COMPACT_CONTAINER_DEBUG_TIME_STAMP
      : ts(0)
#endif
    {
      //CGAL_assertion (n == nullptr);
      m_ptr = nullptr;
    }

  private:

    typedef typename DSC::Time_stamper           Time_stamper;
#ifdef CGAL_COMPACT_CONTAINER_DEBUG_TIME_STAMP
    std::size_t ts;
#endif
    pointer m_ptr;

    // Only Compact_container and Concurrent_compact_container should
    // access these constructors.
    template <typename T, typename Al, typename Ip, typename Ts>
    friend class CGAL::Compact_container;

    friend class CGAL::Concurrent_compact_container<value_type,
                                                    typename DSC::Al>;

    // For begin()
    CC_iterator(pointer ptr, int, int)
#ifdef CGAL_COMPACT_CONTAINER_DEBUG_TIME_STAMP
      : ts(0)
#endif
    {
      m_ptr = ptr;
      if (m_ptr == nullptr) // empty container.
        return;

      ++(m_ptr); // if not empty, p = start
      if (DSC::type(m_ptr) == DSC::FREE)
        increment();
#ifdef CGAL_COMPACT_CONTAINER_DEBUG_TIME_STAMP
      else
        ts = Time_stamper::time_stamp(m_ptr);
#endif // CGAL_COMPACT_CONTAINER_DEBUG_TIME_STAMP
    }

    // Construction from raw pointer and for end().
    CC_iterator(pointer ptr, int)
#ifdef CGAL_COMPACT_CONTAINER_DEBUG_TIME_STAMP
      : ts(0)
#endif
    {
      m_ptr = ptr;
#ifdef CGAL_COMPACT_CONTAINER_DEBUG_TIME_STAMP
      if(ptr != nullptr){
        ts = Time_stamper::time_stamp(m_ptr);
      }
#endif // end CGAL_COMPACT_CONTAINER_DEBUG_TIME_STAMP
    }

    // NB : in case empty container, begin == end == nullptr.
    void increment()
    {
      // It's either pointing to end(), or valid.
      CGAL_assertion_msg(m_ptr != nullptr,
         "Incrementing a singular iterator or an empty container iterator ?");
      CGAL_assertion_msg(DSC::type(m_ptr) != DSC::START_END,
         "Incrementing end() ?");

      // If it's not end(), then it's valid, we can do ++.
      do {
        ++(m_ptr);
        if (DSC::type(m_ptr) == DSC::USED ||
            DSC::type(m_ptr) == DSC::START_END)
        {
#ifdef CGAL_COMPACT_CONTAINER_DEBUG_TIME_STAMP
          ts = Time_stamper::time_stamp(m_ptr);
#endif
          return;
        }
        if (DSC::type(m_ptr) == DSC::BLOCK_BOUNDARY)
          m_ptr = DSC::clean_pointee(m_ptr);
      } while (true);
    }

    void decrement()
    {
      // It's either pointing to end(), or valid.
      CGAL_assertion_msg(m_ptr != nullptr,
         "Decrementing a singular iterator or an empty container iterator ?");
      CGAL_assertion_msg(DSC::type(m_ptr - 1) != DSC::START_END,
         "Decrementing begin() ?");

      // If it's not begin(), then it's valid, we can do --.
      do {
        --m_ptr;
        if (DSC::type(m_ptr) == DSC::USED ||
            DSC::type(m_ptr) == DSC::START_END)
        {
#ifdef CGAL_COMPACT_CONTAINER_DEBUG_TIME_STAMP
          ts = Time_stamper::time_stamp(m_ptr);
#endif
          return;
        }

        if (DSC::type(m_ptr) == DSC::BLOCK_BOUNDARY)
          m_ptr = DSC::clean_pointee(m_ptr);
      } while (true);
    }

  public:

    Self & operator++()
    {
      CGAL_assertion_msg(m_ptr != nullptr,
         "Incrementing a singular iterator or an empty container iterator ?");
      /* CGAL_assertion_msg(DSC::type(m_ptr) == DSC::USED,
         "Incrementing an invalid iterator."); */
      increment();
      return *this;
    }

    Self & operator--()
    {
      CGAL_assertion_msg(m_ptr != nullptr,
         "Decrementing a singular iterator or an empty container iterator ?");
      /*CGAL_assertion_msg(DSC::type(m_ptr) == DSC::USED
                      || DSC::type(m_ptr) == DSC::START_END,
                      "Decrementing an invalid iterator.");*/
      decrement();
      return *this;
    }

    Self operator++(int) { Self tmp(*this); ++(*this); return tmp; }
    Self operator--(int) { Self tmp(*this); --(*this); return tmp; }

#ifdef CGAL_COMPACT_CONTAINER_DEBUG_TIME_STAMP
    bool is_time_stamp_valid() const
    {
      return (ts == 0) || (ts == Time_stamper::time_stamp(m_ptr));
    }
#endif // CGAL_COMPACT_CONTAINER_DEBUG_TIME_STAMP

    reference operator*() const { return *(m_ptr); }

    pointer   operator->() const { return (m_ptr); }

    // For std::less...
    bool operator<(const CC_iterator& other) const
    {
#ifdef CGAL_COMPACT_CONTAINER_DEBUG_TIME_STAMP
      assert( is_time_stamp_valid() );
#endif
      return Time_stamper::less(m_ptr, other.m_ptr);
    }

    bool operator>(const CC_iterator& other) const
    {
#ifdef CGAL_COMPACT_CONTAINER_DEBUG_TIME_STAMP
      assert( is_time_stamp_valid() );
#endif
      return Time_stamper::less(other.m_ptr, m_ptr);
    }

    bool operator<=(const CC_iterator& other) const
    {
#ifdef CGAL_COMPACT_CONTAINER_DEBUG_TIME_STAMP
      assert( is_time_stamp_valid() );
#endif
      return Time_stamper::less(m_ptr, other.m_ptr)
          || (*this == other);
    }

    bool operator>=(const CC_iterator& other) const
    {
#ifdef CGAL_COMPACT_CONTAINER_DEBUG_TIME_STAMP
      assert( is_time_stamp_valid() );
#endif
      return Time_stamper::less(other.m_ptr, m_ptr)
          || (*this == other);
    }

    // Can itself be used for bit-squatting.
    void * for_compact_container() const { return m_ptr; }
    void for_compact_container(void* p) { m_ptr = static_cast<pointer>(p); }
  };

  template < class DSC, bool Const1, bool Const2 >
  inline
  bool operator==(const CC_iterator<DSC, Const1> &rhs,
                  const CC_iterator<DSC, Const2> &lhs)
  {
    return rhs.operator->() == lhs.operator->();
  }

  template < class DSC, bool Const1, bool Const2 >
  inline
  bool operator!=(const CC_iterator<DSC, Const1> &rhs,
                  const CC_iterator<DSC, Const2> &lhs)
  {
    return rhs.operator->() != lhs.operator->();
  }

  // Comparisons with nullptr are part of CGAL's Handle concept...
  template < class DSC, bool Const >
  inline
  bool operator==(const CC_iterator<DSC, Const> &rhs,
                  std::nullptr_t /*CGAL_assertion_code(n)*/)
  {
    //CGAL_assertion( n == nullptr);
    return rhs.operator->() == nullptr;
  }

  template < class DSC, bool Const >
  inline
  bool operator!=(const CC_iterator<DSC, Const> &rhs,
                  std::nullptr_t /*CGAL_assertion_code(n)*/)
  {
    //CGAL_assertion( n == nullptr);
    return rhs.operator->() != nullptr;
  }

  template <class DSC, bool Const>
  std::size_t hash_value(const CC_iterator<DSC, Const>&  i)
  {
    typedef Time_stamper_impl<typename DSC::value_type> Stamper;
    return Stamper::hash_value(i.operator->());
  }

namespace handle {
  // supply a specialization for Hash_functor

  // forward declare base template
  template <class H> struct Hash_functor;

  template<class DSC, bool Const>
  struct Hash_functor<CC_iterator<DSC, Const> >{
    std::size_t
    operator()(const CC_iterator<DSC, Const>& i)
    {
      return hash_value(i);
    }
  };
} // namespace handle

} // namespace internal

} //namespace CGAL

namespace std {

#ifndef CGAL_CFG_NO_STD_HASH

  template < class DSC, bool Const >
  struct hash<CGAL::internal::CC_iterator<DSC, Const> >
    : public CGAL::cpp98::unary_function<CGAL::internal::CC_iterator<DSC, Const>, std::size_t> {

    std::size_t operator()(const CGAL::internal::CC_iterator<DSC, Const>& i) const
    {
      return reinterpret_cast<std::size_t>(&*i) / sizeof(typename DSC::value_type);
    }
  };
#endif // CGAL_CFG_NO_STD_HASH


} // namespace std

#include <CGAL/enable_warnings.h>

#endif // CGAL_COMPACT_CONTAINER_H
