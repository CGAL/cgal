/// \defgroup inplacelist Doubly-Connected List Managing Items in Place
/// \ingroup PkgStlExtension


namespace CGAL {

/*!
\ingroup inplacelist



The node base classes provides pointers to build 
linked lists. The class `In_place_sl_list_base<T>` provides 
a pointer `next_link` for a single linked list. The class 
`In_place_list_base<T>` provides an additional pointer 
`prev_link` for doubly linked lists. These names conform to 
the default parameters used in the template argument lists of the 
container classes. The pointers are public members. 
*/
template< typename T >
class In_place_list_base {
public:

/// \name Variables 
/// @{ 
/*! 
forward pointer 
*/ 
T* next_link; 



/// @} 


/// \name Variables 
/// @{ 
/*! 
backward pointer 
*/ 
T* prev_link; 



/// @} 



}; /* end In_place_list_base */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup inplacelist



An object of the class `In_place_list` 
represents a sequence of items of type `T` that supports 
bidirectional iterators and allows constant time insert and erase 
operations anywhere within the sequence. The functionality is 
similar to the `list<T>` in the \stl. 

The `In_place_list` manages the items in place, i.e., inserted 
items are not copied. Two pointers of type `T*` are expected 
to be reserved in `T` for the list management. The base class 
`In_place_list_base<T>` can be used to obtain such pointers. 

The `In_place_list` does not copy element items during 
insertion (unless otherwise stated for a function). On removal of an 
item or destruction of the list the items are not deleted by 
default. The second template parameter `bool` is set to 
`false` in this case. If the `In_place_list` should 
take the responsibility for the stored objects the `bool` 
parameter could be set to `true`, in which case the list 
will delete removed items and will delete all remaining items on 
destruction. In any case, the `destroy()` member function 
deletes all items. Note that these two possible versions of 
`In_place_list` are not assignable to each other to avoid 
confusions between the different storage responsibilities. 

### Parameters ###

The full class name is `In_place_list<T, bool managed = 
false, class Alloc = CGAL_ALLOCATOR(T)>`. 

The parameter `T` is supposed to have a default constructor, 
a copy constructor and an assignment operator. The copy constructor 
and the assignment may copy the pointers in `T` for the list 
management, but they do not have to. The equality test and the 
relational order require the operators `==` and `<` 
for \f$ T\f$ respectively. These operators must not compare the pointers 
in `T`. 

### Example ###

\cgalexample{STL_Extension/in_place_list_prog.cpp} 


*/
template< typename T, bool >
class In_place_list {
public:


/// \name Types 
/// @{ 
/*! 

*/ 
typedef Hidden_type iterator; 
/// @} 


/// \name Types 
/// @{ 
/*! 

*/ 
typedef Hidden_type const_iterator; 
/// @} 


/// \name Types 
/// @{ 
/*! 

*/ 
typedef Hidden_type value_type; 
/// @} 


/// \name Types 
/// @{ 
/*! 

*/ 
typedef Hidden_type reference; 
/// @} 


/// \name Types 
/// @{ 
/*! 

*/ 
typedef Hidden_type const_reference; 
/// @} 


/// \name Types 
/// @{ 
/*! 

*/ 
typedef Hidden_type size_type; 
/// @} 


/// \name Types 
/// @{ 
/*! 

*/ 
typedef Hidden_type difference_type; 
/// @} 


/// \name Types 
/// @{ 
/*! 

*/ 
typedef Hidden_type reverse_iterator; 
/// @} 


/// \name Types 
/// @{ 
/*! 

*/ 
typedef Hidden_type const_reverse_iterator; 
/// @} 


/// \name Types 
/// @{ 
/*! 

*/ 
typedef Hidden_type allocator_type; 
/// @} 

/// \name Creation 
/// @{ 
/*! 
introduces an empty list. 
*/ 
In_place_list(); 



/// @} 


/// \name Creation 
/// @{ 
/*! 
copy constructor. 
Each item in `l1` is copied. 
*/ 
In_place_list(const list<T> &l1); 



/// @} 


/// \name Creation 
/// @{ 
/*! 
introduces a list with \f$ n\f$ items, all initialized with copies 
of \f$ t\f$. 
*/ 
In_place_list(size_type n, const T& t = T() 
); 



/// @} 


/// \name Creation 
/// @{ 
/*! 
a list with copies from 
the range [`first,last`). 
*/ 
template <class InputIterator> In_place_list( 
InputIterator first, InputIterator last); 



/// @} 


/// \name Creation 
/// @{ 
/*! 
non-member-template version. 
*/ 
In_place_list( const T* first, const T* 
last); 



/// @} 


/// \name Creation 
/// @{ 
/*! 
assignment. Each item in `l1` 
is copied. Each item in `l` is deleted if the `bool` 
parameter is `true`. 
*/ 
In_place_list<T,bool> & operator=(const 
In_place_list<T,bool> &l1); 



/// @} 


/// \name Creation 
/// @{ 
/*! 
swaps the 
contents of `l` with `l1`. 
*/ 
void swap( const In_place_list<T,bool> &l1); 



/// @} 


/// \name Creation 
/// @{ 
/*! 
all items in `l` are deleted 
regardless of the `bool` parameter. 
*/ 
void destroy(); 



/// @} 


/// \name Comparison Operations 
/// @{ 
/*! 
test for equality: Two lists are equal, iff they have the 
same size and if their corresponding elements are equal. 
*/ 
bool operator==(const In_place_list<T,bool> &l1) 
const; 



/// @} 


/// \name Comparison Operations 
/// @{ 
/*! 
compares in lexicographical order. 
*/ 
bool operator<(const In_place_list<T,bool> &l1) 
const; 



/// @} 


/// \name Access Member Functions 
/// @{ 
/*! 
returns a mutable iterator referring to the first 
element in `l`. 
*/ 
iterator 
begin(); 



/// @} 


/// \name Access Member Functions 
/// @{ 
/*! 
returns a constant 
iterator referring to the first element in `l`. 
*/ 
const_iterator begin() const; 



/// @} 


/// \name Access Member Functions 
/// @{ 
/*! 
returns a mutable iterator which 
is the past-end-value of `l`. 
*/ 
iterator end(); 



/// @} 


/// \name Access Member Functions 
/// @{ 
/*! 
returns a constant 
iterator which is the past-end-value of `l`. 
*/ 
const_iterator end() const; 



/// @} 


/// \name Access Member Functions 
/// @{ 
/*! 
returns `true` if `l` is 
empty. 
*/ 
bool empty() const; 



/// @} 


/// \name Access Member Functions 
/// @{ 
/*! 
returns the number of 
items in list `l`. 
*/ 
size_type size() const; 



/// @} 


/// \name Access Member Functions 
/// @{ 
/*! 
returns the maximum 
possible size of the list `l`. 
*/ 
size_type max_size() const; 



/// @} 


/// \name Access Member Functions 
/// @{ 
/*! 
returns the first item in list `l`. 
*/ 
T& front(); 



/// @} 


/// \name Access Member Functions 
/// @{ 
/*! 
returns the last item in list `l`. 
*/ 
T& back(); 



/// @} 


/// \name Access Member Functions 
/// @{ 
/*! 
returns the allocator. 
*/ 
allocator_type get_allocator() const; 



/// @} 


/// \name Insertion 
/// @{ 
/*! 
inserts an item in front of 
list `l`. 
*/ 
void push_front( T&); 



/// @} 


/// \name Insertion 
/// @{ 
/*! 
inserts an item at the back 
of list `l`. 
*/ 
void push_back( T&); 



/// @} 


/// \name Insertion 
/// @{ 
/*! 

*/ 
iterator insert(iterator pos, T& t); 



/// @} 


/// \name Insertion 
/// @{ 
/*! 
inserts `t` 
in front of `pos`. The return value points to the 
inserted item. 
*/ 
iterator insert(T* pos, T& t); 



/// @} 


/// \name Insertion 
/// @{ 
/*! 

*/ 
void insert(iterator pos, size_type n, const T& t = T()); 



/// @} 


/// \name Insertion 
/// @{ 
/*! 
inserts \f$ n\f$ copies of `t` in front of 
`pos`. 
*/ 
void insert(T* pos, size_type n, const T& t = 
T()); 



/// @} 


/// \name Insertion 
/// @{ 
/*! 

*/ 
template <class InputIterator> void insert(iterator pos, 
InputIterator first, InputIterator last); 



/// @} 


/// \name Insertion 
/// @{ 
/*! 
inserts the range 
[`first, last`) in front of iterator `pos`. 
*/ 
template <class InputIterator> void insert(T* pos, 
InputIterator first, InputIterator last); 



/// @} 


/// \name Removal 
/// @{ 
/*! 
removes the first item from 
list `l`. 
*/ 
void pop_front(); 



/// @} 


/// \name Removal 
/// @{ 
/*! 
removes the last item from 
list `l`. 
*/ 
void pop_back(); 



/// @} 


/// \name Removal 
/// @{ 
/*! 
removes the item from 
list `l`, where `pos` refers to. 
*/ 
void erase(iterator pos); 



/// @} 


/// \name Removal 
/// @{ 
/*! 
removes the item from 
list `l`, where `pos` refers to. 
*/ 
void erase(T* pos); 



/// @} 


/// \name Removal 
/// @{ 
/*! 

*/ 
void erase(iterator first, iterator last); 



/// @} 


/// \name Removal 
/// @{ 
/*! 
removes the items 
in the range [`first, last`) from `l`. 
*/ 
void erase(T* first, T* last); 



/// @} 


/// \name Special List Operations 
/// @{ 
/*! 

*/ 
void splice(iterator pos, In_place_list<T,bool>& x); 



/// @} 


/// \name Special List Operations 
/// @{ 
/*! 
inserts the list \f$ x\f$ before position `pos` and \f$ x\f$ 
becomes empty. It takes constant time. \pre `&` `l``!= &x`. 
*/ 
void splice(T* pos, In_place_list<T,bool>& 
x); 



/// @} 


/// \name Special List Operations 
/// @{ 
/*! 

*/ 
void splice(iterator pos, In_place_list<T,bool>& x, 
iterator i); 



/// @} 


/// \name Special List Operations 
/// @{ 
/*! 
inserts an element pointed to by \f$ i\f$ from list \f$ x\f$ before 
position `pos` and removes the element from \f$ x\f$. It takes 
constant time. \f$ i\f$ is a valid dereferenceable iterator of \f$ x\f$. 
The result is unchanged if `pos == i` or `pos == 
++i`. 
*/ 
void splice(T* pos, In_place_list<T,bool>& x, T* 
i); 



/// @} 


/// \name Special List Operations 
/// @{ 
/*! 

*/ 
void splice(iterator pos, In_place_list<T,bool>& x, 
iterator first, iterator last); 



/// @} 


/// \name Special List Operations 
/// @{ 
/*! 
inserts elements in the range [`first, 
last`) before position `pos` and removes the elements 
from \f$ x\f$. It takes constant time if `&x == &``l`; 
otherwise, it takes linear time. [`first, last`) is a 
valid range in \f$ x\f$. \pre `pos` is not in the range [`first, last`). 
*/ 
void splice(T* pos, In_place_list<T,bool>& x, T* 
first, T* last); 



/// @} 


/// \name Special List Operations 
/// @{ 
/*! 
erases all elements \f$ e\f$ in 
the list `l` for which `e == value`. It is stable. 
\pre a suitable `operator==` for the type \f$ T\f$. 
*/ 
void remove(const T& value); 



/// @} 


/// \name Special List Operations 
/// @{ 
/*! 
erases all but the first element from 
every consecutive group of equal elements in the list `l`. 
\pre a suitable `operator==` for the type \f$ T\f$. 
*/ 
void unique(); 



/// @} 


/// \name Special List Operations 
/// @{ 
/*! 
merges the list \f$ x\f$ 
into the list `l` and \f$ x\f$ becomes empty. It is stable. 
\pre Both lists are increasingly sorted. A suitable `operator<` for the type \f$ T\f$. 
*/ 
void merge(In_place_list<T,bool>& x); 



/// @} 


/// \name Special List Operations 
/// @{ 
/*! 
reverses the order of the elements in 
`l` in linear time. 
*/ 
void reverse(); 



/// @} 


/// \name Special List Operations 
/// @{ 
/*! 
sorts the list `l` according to the 
`operator<` in time \f$ O(n \logn)\f$ where `n = 
size()`. It is stable. \pre a suitable `operator<` for the type \f$ T\f$. 
*/ 
void sort(); 



/// @} 



}; /* end In_place_list */
} /* end namespace CGAL */
