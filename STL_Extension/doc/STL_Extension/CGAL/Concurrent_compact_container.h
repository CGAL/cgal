
namespace CGAL {

/// \ingroup PkgStlExtension

/*!
\ingroup CompactContainer

The traits class `Concurrent_compact_container_traits` provides 
the way to access the internal pointer required for `T` to be 
used in a `Concurrent_compact_container<T, Allocator>`. Note that this 
pointer needs to be accessible even when the object is not constructed, 
which means it has to reside in the same memory place as `T`. 

You can specialize this class for your own type `T` 
if the default template is not suitable. 

You can also use `Compact_container_base` as base class for your own 
types `T` to make them usable with the default `Concurrent_compact_container`. 

\cgalHeading{Parameters}

`T` is any type providing the following member functions: 
`void * t.for_compact_container() const;` 
`void *& t.for_compact_container();`.
*/
template< typename T >
class Concurrent_compact_container_traits {
public:

/// \name Operations 
/// @{ 
  /*! 
  Returns the pointer held by `t`. 
  The template version defines this function as: `return t.for_compact_container(); 
  */ 
  static void * pointer(const T &t); 

/// @} 

/// \name Operations 
/// @{ 
  /*! 
  Returns a reference to the pointer held by `t`. 
  The template version defines this function as: `return t.for_compact_container();` 

  */ 
  static void * & pointer(T &t); 

/// @} 

}; /* end Concurrent_compact_container_traits */


/*!
\ingroup CompactContainer

An object of the class `Concurrent_compact_container` 
is a container of objects of type `T`, which allows to call
`insert` and `erase` operations concurrently.
Other operations are not concurrency-safe.
For example, one should not parse the container while others are modifying it.
It matches all the 
standard requirements for reversible containers, except that 
the complexity of its iterator increment and decrement operations 
is not always guaranteed to be amortized constant time. 

This container is not a standard <I>sequence</I> nor <I>associative</I> container, 
which means the elements are stored in no particular order, and it is not 
possible to specify a particular place in the iterator sequence where to 
insert new objects. However, all dereferenceable iterators are 
still valid after calls to `insert()` and `erase()`, except those 
that have been erased (it behaves similarly to `std::list`). 

The main feature of this container is that it is very memory efficient: 
its memory size is `N*sizeof(T)+o(N)`, where `N` is the maximum size 
that the container has had in its past history, its `capacity()` 
(the memory of erased elements is not deallocated until destruction of the 
container or a call to `clear()`). This container has been developed in 
order to store large graph-like data structures like the triangulation and 
the halfedge data structures. 

It supports bidirectional iterators and allows a constant time amortized 
`insert()` operation. You cannot specify where to insert new objects 
(i.e.\ you don't know where they will end up in the iterator sequence, 
although `insert()` returns an iterator pointing to the newly inserted 
object). You can erase any element with a constant time complexity. 

Summary of the differences with `std::list`: it is more compact in 
memory since it doesn't store two additional pointers for the iterator needs. 
It doesn't deallocate elements until the destruction or `clear()` of the 
container. The iterator does not have constant amortized time complexity for 
the increment and decrement operations in all cases, only when not too many 
elements have not been freed (i.e.\ when the `size()` is close to the 
`capacity()`). Iterating from `begin()` to `end()` takes 
`O(capacity())` time, not `size()`. In the case where the container 
has a small `size()` compared to its `capacity()`, we advise to 
\"defragment the memory\" by copying the container if the iterator performance 
is needed. 

The iterators themselves can be used as `T`, they provide the necessary 
functions to be used by `Compact_container_traits<T>`. Moreover, they 
also provide a default constructor value which is not singular: it is 
copyable, comparable, and guaranteed to be unique under comparison 
(like `NULL` for pointers). This makes them suitable for use in 
geometric graphs like handles to vertices in triangulations. 

In addition, in a way inspired from the Boost.Intrusive containers, it is 
possible to construct iterators from references to values in containers 
using the `iterator_to` and `s_iterator_to` functions. 

The objects stored in the `Concurrent_compact_container` can optionally store an 
"erase counter". If it exists, i.e.\ if the object is a model of the
`ObjectWithEraseCounter` concept, each time an object is erased from the 
container, the erase counter of the object will be incremented.
For example, this erase counter can be exploited using the `CC_safe_handle` 
helper class, so that one can know if a handle is still pointing to the same
element.
Note that this is meaningful only because the 
`CGAL::Concurrent_compact_container` doesn't 
deallocate elements until the destruction or clear() of the container.

\cgalHeading{Parameters}

The parameter `T` is required to have a copy constructor and an 
assignment operator. It also needs to provide access to an internal 
pointer via `Compact_container_traits<T>`. 

The equality test and the relational order require the operators 
`==` and `<` for `T` respectively. 

The parameter `Allocator` has to match the standard allocator 
requirements, with value type `T`. This parameter has the default 
value `CGAL_ALLOCATOR(T)`.

*/

template < class T, class Allocator >
class Concurrent_compact_container
{
public:
/// \name Types 
/// @{ 
  typedef unspecified_type value_type;
  typedef unspecified_type allocator_type;
  typedef unspecified_type reference;
  typedef unspecified_type const_reference;
  typedef unspecified_type pointer;
  typedef unspecified_type const_pointer;
  typedef unspecified_type size_type;
  typedef unspecified_type difference_type;
  typedef unspecified_type iterator;
  typedef unspecified_type const_iterator;
  typedef unspecified_type reverse_iterator;
  typedef unspecified_type const_reverse_iterator;
/// @} 


/// \name Creation 
/// @{ 
/*! 
introduces an empty container `ccc`, eventually specifying a particular 
allocator `a` as well. 
*/ 
  explicit Concurrent_compact_container(const Allocator &a = Allocator());

/*! 
a container with copies from the range [`first,last`), eventually 
specifying a particular allocator. 
*/ 
  template < class InputIterator >
  Concurrent_compact_container(InputIterator first, InputIterator last,
                    const Allocator & a = Allocator());

/*! 
copy constructor. Each item in `ccc2` is copied. The allocator 
is copied. The iterator order is preserved. 
*/ 
  // The copy constructor and assignment operator preserve the iterator order
  Concurrent_compact_container(const Concurrent_compact_container &ccc2);
  
/*! 
assignment. Each item in `ccc2` is copied. The allocator is copied. 
Each item in `ccc` is deleted. The iterator order is preserved. 
*/ 
  Concurrent_compact_container & operator=(const Concurrent_compact_container &ccc2);

/*! 
swaps the contents of `ccc` and `ccc2` in constant time 
complexity. No exception is thrown. 
*/ 
  void swap(Self &ccc2);
  
/// @} 

/// \name Access Member Functions 
/// @{ 

  /*!
  returns true if the element `pos` is used (i.e.\ valid).
  */
  bool is_used(const_iterator pos) const;

  /// returns a mutable iterator referring to the first element in `ccc`.
  iterator begin();
  /// returns a constant iterator referring to the first element in `ccc`.
  const_iterator begin() const;
  /// returns a mutable iterator which is the past-end-value of `ccc`.
  iterator end();
  /// returns a constant iterator which is the past-end-value of `ccc`. 
  const_iterator end();

  /// returns a mutable reverse iterator referring to the reverse beginning in `ccc`. 
  reverse_iterator rbegin();
  /// returns a constant reverse iterator referring to the reverse beginning in `ccc`.
  const_reverse_iterator rbegin() const;
  /// returns a mutable reverse iterator which is the reverse past-end-value of `ccc`.
  reverse_iterator rend();
  /// returns a constant reverse iterator which is the reverse past-end-value of `ccc`. 
  const_reverse_iterator rend() const;

  /// returns an iterator which points to `value`.
  iterator iterator_to(reference value) const;
  /// returns a constant iterator which points to `value`.
  const_iterator iterator_to(const_reference value) const;
  /// returns an iterator which points to `value`.
  static iterator s_iterator_to(reference value);
  /// returns a constant iterator which points to `value`.
  static const_iterator s_iterator_to(const_reference value);

  /// returns `true` iff `ccc` is empty. 
  bool empty() const;
  /// returns the number of items in `ccc`. 
  /// Note: do not call this function while others are inserting/erasing elements
  size_type size() const;
  /// returns the maximum possible size of the container `ccc`.
  /// This is the allocator's max_size value
  size_type max_size() const; 
  /// returns the total number of elements that `ccc` can hold without requiring reallocation. 
  size_type capacity() const;
  /// returns the allocator
  Allocator get_allocator() const; 

/// @} 


/// \name Insertion 
/// @{ 
  /*! 
  constructs an object of type `T` with the constructor that takes 
  `t1` as argument, inserts it in `ccc`, and returns the iterator pointing 
  to it. Overloads of this member function are defined that take additional 
  arguments, up to 9. 
  */ 
  template < class T1 > 
  iterator emplace(const T1& t1); 
  
  /*! 
  inserts a copy of `t` in `ccc` and returns the iterator pointing 
  to it. 
  */ 
  iterator insert(const T &t);

  /// inserts the range [`first, last`) in `ccc`.
  template < class InputIterator >
  void insert(InputIterator first, InputIterator last);

  /*! 
  erases all the elements of `ccc`, then inserts the range 
  [`first, last`) in `ccc`. 
  */
  template < class InputIterator >
  void assign(InputIterator first, InputIterator last);
/// @} 


/// \name Removal 
/// @{
  /// removes the item pointed by `pos` from `ccc`.
  void erase(iterator x);
  /// removes the items from the range [`first, last`) from `ccc`.
  void erase(iterator first, iterator last);
  /*! 
  all items in `ccc` are deleted, and the memory is deallocated. 
  After this call, `ccc` is in the same state as if just default 
  constructed. 
  */ 
  void clear();
/// @} 

/// \name Ownership testing 
/// The following functions are mostly helpful for efficient debugging, since 
/// their complexity is \f$ O(\sqrt{\mathrm{c.capacity()}})\f$. 
/// @{ 
  /// returns whether `pos` is in the range `[ccc.begin(),  ccc.end()]` (`ccc.end()` included).
  bool owns(const_iterator pos); 
  /// returns whether `pos` is in the range `[ccc.begin(), ccc`.end())` (`ccc.end()` excluded). 
  bool owns_dereferencable(const_iterator pos); 
  
/// @} 

/// \name Merging 
/// @{ 
/*! 
adds the items of `ccc2` to the end of `ccc` and `ccc2` becomes empty. 
The time complexity is O(`ccc`.`capacity()`-`ccc`.`size()`). 
\pre `ccc2` must not be the same as `ccc`, and the allocators of `ccc` and `ccc2` must be compatible: `ccc.get_allocator() == ccc2.get_allocator()`. 
*/ 
void merge(Concurrent_compact_container<T, Allocator> &ccc2); 

/// @}
  
/// \name Comparison Operations 
/// @{ 
  /*! 
  test for equality: Two containers are equal, iff they have the 
  same size and if their corresponding elements are equal. 
  */ 
  bool operator==(const Concurrent_compact_container<T, Allocator> &ccc2) const; 
  /// test for inequality: returns `!(ccc == ccc2)`. 
  bool operator!=(const Concurrent_compact_container<T, Allocator> &ccc2) const;
  /// compares in lexicographical order. 
  bool operator<(const Concurrent_compact_container<T, Allocator> &ccc2) const; 
  /// returns `ccc2 < ccc`.
  bool operator>(const Concurrent_compact_container<T, Allocator> &ccc2) const;
  /// returns `!(ccc > ccc2)`.
  bool operator<=(const Concurrent_compact_container<T, Allocator> &ccc2) const;
  /// returns `!(ccc < ccc2)`.
  bool operator>=(const Concurrent_compact_container<T, Allocator> &ccc2) const;
/// @} 

}; /* end Concurrent_compact_container */
} /* end namespace CGAL */
