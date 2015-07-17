
namespace CGAL {

/// \defgroup CompactContainer Compact Container
/// \ingroup PkgStlExtension

/*!
\ingroup CompactContainer



The class `Compact_container_base` can be used as a base class for 
your own type `T`, so that `T` can be used directly within 
`Compact_container<T, Allocator>`. This class stores a `void *` 
pointer only for this purpose, so it may not be the most memory efficient 
way to achieve this goal. The other ways are to provide in `T` the 
necessary member functions so that the template 
`Compact_container_traits<T>` works, or to specialize it for the 
particular type `T` that you want to use. 

*/
class Compact_container_base {
public:

/// \name Operations 
/// @{ 
/*!
Returns the pointer necessary for `Compact_container_traits<T>`. 
*/ 
void * for_compact_container() const; 

/*!
Returns a reference to the pointer necessary for 
`Compact_container_traits<T>`. 
*/ 
void * & for_compact_container(); 

/// @} 



}; /* end Compact_container_base */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup CompactContainer

An object of the class `Compact_container` 
is a container of objects of type `T`. 

This container matches all the 
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
"defragment the memory" by copying the container if the iterator performance 
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

The objects stored in the `Compact_container` can optionally store an 
"erase counter". If it exists, i.e.\ if the object is a model of the
`ObjectWithEraseCounter` concept, each time an object is erased from the 
container, the erase counter of the object will be incremented.
For example, this erase counter can be exploited using the `CC_safe_handle` 
helper class, so that one can know if a handle is still pointing to the same
element.
Note that this is meaningful only because the 
`CGAL::Compact_container` doesn't 
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
template< typename T, typename Allocator >
class Compact_container {
public:


/// \name Types 
/// @{ 
/*!

*/ 
typedef unspecified_type value_type; 
/// @} 


/// \name Types 
/// @{ 
/*!

*/ 
typedef unspecified_type reference; 
/// @} 


/// \name Types 
/// @{ 
/*!

*/ 
typedef unspecified_type const_reference; 
/// @} 


/// \name Types 
/// @{ 
/*!

*/ 
typedef unspecified_type pointer; 
/// @} 


/// \name Types 
/// @{ 
/*!

*/ 
typedef unspecified_type const_pointer; 
/// @} 


/// \name Types 
/// @{ 
/*!

*/ 
typedef unspecified_type size_type; 
/// @} 


/// \name Types 
/// @{ 
/*!

*/ 
typedef unspecified_type difference_type; 
/// @} 


/// \name Types 
/// @{ 
/*!

*/ 
typedef unspecified_type iterator; 
/// @} 


/// \name Types 
/// @{ 
/*!

*/ 
typedef unspecified_type const_iterator; 
/// @} 


/// \name Types 
/// @{ 
/*!

*/ 
typedef unspecified_type reverse_iterator; 
/// @} 


/// \name Types 
/// @{ 
/*!

*/ 
typedef unspecified_type const_reverse_iterator; 
/// @} 


/// \name Types 
/// @{ 
/*!

*/ 
typedef unspecified_type allocator_type; 
/// @} 

/// \name Creation 
/// @{ 
/*!
introduces an empty container `cc`, eventually specifying a particular 
allocator `a` as well. 
*/ 
explicit Compact_container(const Allocator &a = Allocator()); 



/// @} 


/// \name Creation 
/// @{ 
/*!
a container with copies from the range [`first,last`), eventually 
specifying a particular allocator. 
*/ 
template <class InputIterator> Compact_container( 
InputIterator first, InputIterator last, 
const Allocator &a = Allocator()); 



/// @} 


/// \name Creation 
/// @{ 
/*!
copy constructor. Each item in `cc2` is copied. The allocator 
is copied. The iterator order is preserved. 
*/ 
Compact_container(const Compact_container<T, Allocator> &cc2); 



/// @} 


/// \name Creation 
/// @{ 
/*!
assignment. Each item in `cc2` is copied. The allocator is copied. 
Each item in `c` is deleted. The iterator order is preserved. 
*/ 
Compact_container<T, Allocator> & operator=(const 
Compact_container<T, Allocator> &cc2); 



/// @} 


/// \name Creation 
/// @{ 
/*!
swaps the contents of `cc` and `cc2` in constant time 
complexity. No exception is thrown. 
*/ 
void swap(Compact_container<T, Allocator> &cc2); 



/// @} 


/// \name Creation 
/// @{ 
/*!
if `value` is less than or equal to `capacity()`, this call 
has no effect. Otherwise, it is a request for allocation of 
additional memory so that then `capacity()` is greater than or 
equal to value. `size()` is unchanged. 
*/ 
void reserve(size_type value); 



/// @} 


/// \name Access Member Functions 
/// @{ 
/*!
returns a mutable iterator referring to the first element in `cc`. 
*/ 
iterator begin(); 



/// @} 


/// \name Access Member Functions 
/// @{ 
/*!
returns a constant iterator referring to the first element in `cc`. 
*/ 
const_iterator begin() const; 



/// @} 


/// \name Access Member Functions 
/// @{ 
/*!
returns a mutable iterator which is the past-end-value of `cc`. 
*/ 
iterator end(); 



/// @} 


/// \name Access Member Functions 
/// @{ 
/*!
returns a constant iterator which is the past-end-value of `cc`. 
*/ 
const_iterator end() const; 



/// @} 


/// \name Access Member Functions 
/// @{ 
/*!

*/ 
reverse_iterator rbegin(); 



/// @} 


/// \name Access Member Functions 
/// @{ 
/*!

*/ 
const_reverse_iterator rbegin() const; 



/// @} 


/// \name Access Member Functions 
/// @{ 
/*!

*/ 
reverse_iterator rend(); 



/// @} 


/// \name Access Member Functions 
/// @{ 
/*!

*/ 
const_reverse_iterator rend() const; 



/// @} 


/// \name Access Member Functions 
/// @{ 
/*!
returns an iterator which points to `value`. 
*/ 
iterator iterator_to(reference value) const; 



/// @} 


/// \name Access Member Functions 
/// @{ 
/*!
returns an iterator which points to `value`. 
*/ 
const_iterator iterator_to(const_reference value) const; 



/// @} 


/// \name Access Member Functions 
/// @{ 
/*!
returns an iterator which points to `value`; 
*/ 
static iterator s_iterator_to(reference value); 



/// @} 


/// \name Access Member Functions 
/// @{ 
/*!
returns an iterator which points to `value`; 
*/ 
static const_iterator s_iterator_to(const_reference value); 



/// @} 


/// \name Access Member Functions 
/// @{ 
/*!
returns `true` iff `cc` is empty. 
*/ 
bool empty() const; 



/// @} 


/// \name Access Member Functions 
/// @{ 
/*!
returns the number of items in `cc`. 
*/ 
size_type size() const; 



/// @} 


/// \name Access Member Functions 
/// @{ 
/*!
returns the maximum possible size of the container `cc`.
This is the allocator's max_size value.
*/ 
size_type max_size() const; 



/// @} 


/// \name Access Member Functions 
/// @{ 
/*!
returns the total number of elements that `cc` can hold without requiring 
reallocation. 
*/ 
size_type capacity() const; 


/// \name Access Member Functions 
/// @{ 

/*!
returns true if the element `pos` is used (i.e.\ valid).
*/
bool is_used(const_iterator pos) const;

/*!
returns true if the element at position `i` in the container is used
(i.e.\ valid).

\pre \f$ 0 \leq \f$ `i` \f$ < \f$ `capacity()`
*/

bool is_used(size_type i) const;

/// @} 

/// \name Access Member Functions 
/// @{ 
/*!
returns the element at pos `i` in the container. 

\pre `is_used(i) == true` and \f$ 0 \leq \f$ `i` \f$ < \f$ `capacity()`
*/ 

const T& operator[] (size_type i) const;

/// @}


/// \name Access Member Functions 
/// @{ 
/*!
returns the element at pos `i` in the container. 

\pre `is_used(i) == true` and \f$ 0 \leq \f$ `i` \f$ < \f$ `capacity()`
*/ 

T& operator[] (size_type i);

/// @} 


/// \name Access Member Functions 
/// @{ 
/*!
returns the allocator. 
*/ 
Allocator get_allocator() const; 



/// @} 


/// \name Insertion 
/// @{ 
/*!
inserts a copy of `t` in `cc` and returns the iterator pointing 
to it. 
*/ 
iterator insert(const T& t); 



/// @} 


/// \name Insertion 
/// @{ 
/*!
inserts the range [`first, last`) in `cc`. 
*/ 
template <class InputIterator> 
void insert(InputIterator first, InputIterator last); 



/// @} 


/// \name Insertion 
/// @{ 
/*!
erases all the elements of `cc`, then inserts the range 
[`first, last`) in `cc`. 
*/ 
template <class InputIterator> 
void assign(InputIterator first, InputIterator last); 



/// @} 


/// \name Insertion 
/// @{ 
/*!
constructs an object of type `T` with the constructor that takes 
`t1` as argument, inserts it in `cc`, and returns the iterator pointing 
to it. Overloads of this member function are defined that take additional 
arguments, up to 9. 
*/ 
template < class T1 > 
iterator emplace(const T1& t1); 



/// @} 


/// \name Removal 
/// @{ 
/*!
removes the item pointed by `pos` from `cc`. 
*/ 
void erase(iterator pos); 



/// @} 


/// \name Removal 
/// @{ 
/*!
removes the items from the range [`first, last`) from `cc`. 
*/ 
void erase(iterator first, iterator last); 



/// @} 


/// \name Removal 
/// @{ 
/*!
all items in `cc` are deleted, and the memory is deallocated. 
After this call, `cc` is in the same state as if just default 
constructed. 
*/ 
void clear(); 



/// @} 


/// \name Ownership testing 
/// The following functions are mostly helpful for efficient debugging, since 
/// their complexity is \f$ O(\sqrt{\mathrm{c.capacity()}})\f$. 
/// @{ 

/*!
 * returns whether `pos` is in the range `[cc.begin(),  cc.end()]` (`cc.end()` included).
 */ 
bool owns(const_iterator pos); 

/*!
 * returns whether `pos` is in the range `[cc.begin(), cc`.end())` (`cc.end()` excluded). 
 */ 
bool owns_dereferencable(const_iterator pos); 

/// @} 


/// \name Merging 
/// @{ 
/*!
adds the items of `cc2` to the end of `cc` and `cc2` becomes empty. 
The time complexity is O(`cc`.`capacity()`-`cc`.`size()`). 
\pre `cc2` must not be the same as `cc`, and the allocators of `cc` and `cc2` must be compatible: `cc.get_allocator() == cc2.get_allocator()`. 
*/ 
void merge(Compact_container<T, Allocator> &cc); 



/// @} 


/// \name Comparison Operations 
/// @{ 
/*!
test for equality: Two containers are equal, iff they have the 
same size and if their corresponding elements are equal. 
*/ 
bool operator==(const Compact_container<T, Allocator> &cc) const; 



/// @} 


/// \name Comparison Operations 
/// @{ 
/*!
test for inequality: returns `!(c == cc)`. 
*/ 
bool operator!=(const Compact_container<T, Allocator> &cc) const; 



/// @} 


/// \name Comparison Operations 
/// @{ 
/*!
compares in lexicographical order. 
*/ 
bool operator<(const Compact_container<T, Allocator> &cc2) const; 



/// @} 


/// \name Comparison Operations 
/// @{ 
/*!
returns `cc2 <cc`. 
*/ 
bool operator>(const Compact_container<T, Allocator> &cc2) const; 



/// @} 


/// \name Comparison Operations 
/// @{ 
/*!
  returns `!(cc > cc2)`. 
*/ 
bool operator<=(const Compact_container<T, Allocator> &cc2) const; 



/// @} 


/// \name Comparison Operations 
/// @{ 
/*!
returns `!(cc < cc2)`. 
*/ 
bool operator>=(const Compact_container<T, Allocator> &cc2) const; 



/// @} 



}; /* end Compact_container */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup CompactContainer



The traits class `Compact_container_traits` provides 
the way to access the internal pointer required for `T` to be 
used in a `Compact_container<T, Allocator>`. Note that this 
pointer needs to be accessible even when the object is not constructed, 
which means it has to reside in the same memory place as `T`. 

You can specialize this class for your own type `T` 
if the default template is not suitable. 

You can also use `Compact_container_base` as base class for your own 
types `T` to make them usable with the default `Compact_container_traits`. 



\cgalHeading{Parameters}

`T` is any type providing the following member functions: 

`void * t.for_compact_container() const;` 

`void *& t.for_compact_container();`. 


*/
template< typename T >
class Compact_container_traits {
public:

/// \name Operations 
/// @{ 
/*!
Returns the pointer held by `t`. 
The template version defines this function as: `return t.for_compact_container(); `

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



}; /* end Compact_container_traits */

/*!
returns a hash value for the pointee of `i`. 
\relates Compact_container
*/ 
  template <class T, class A>
  std::size_t hash_value(const Compact_container<T,A>::iterator i);

} /* end namespace CGAL */
