
namespace CGAL {

/*!
\ingroup PkgProfilingTools

An instance `P` of the data type `Union_find<T,A>` is a 
partition of values of type `T` into disjoint sets. The template 
parameter `A` has to be a model of the allocator concept as defined 
in the C++ standard. It has a default argument `CGAL_ALLOCATOR(T)`. 

### Implementation ###

`Union_find<T,A>` is implemented with union by rank and path 
compression. The running time for \f$ m\f$ set operations on \f$ n\f$ elements 
is \f$ O(n \alpha(m,n))\f$ where \f$ \alpha(m,n)\f$ is the extremely slow growing 
inverse of Ackermann's function. 

*/
template< typename T, typename A >
class Union_find {
public:

/// \name Types 
/// There are also constant versions `const_handle` and `const_iterator`.
/// @{

/*! 
values stored in items (equal to `T`). 
*/ 
typedef Hidden_type value_type; 

/*! 
handle to values. 
*/ 
typedef Hidden_type handle; 

/*! 
iterator over values. 
*/ 
typedef Hidden_type iterator; 

/*! 
allocator. 
*/ 
typedef Hidden_type allocator; 

/// @} 

/// \name Creation 
/// @{

/*! 
creates an instance `P` of type 
`Union_find<T,A>` and initializes it to the empty partition. 
*/ 
Union_find<T,A>(); 

/// @} 

/// \name Operations 
/// @{

/*! 
the allocator of `P`. 
*/ 
allocator get_allocator() ; 

/*! 
returns the number of disjoint 
sets of `P`. 
*/ 
std::size_t number_of_sets() ; 

/*! 
returns the number of values of `P`. 
*/ 
std::size_t size() ; 

/*! 
returns the memory consumed by `P`. 
*/ 
std::size_t bytes() ; 

/*! 
returns the size of the set 
containing \f$ p\f$. 
*/ 
std::size_t size( const_handle p) ; 

/*! 
reinitializes `P` to an empty partition. 
*/ 
void clear(); 

/*! 
creates a new singleton set 
containing `x` and returns a handle to it. 
*/ 
handle make_set(const T& x); 

/*! 
same as `make_set(x)`. 
*/ 
handle push_back(const T& x) ; 

/*! 
insert 
the range of values referenced by `[first,beyond)`. 
\requires value type of `Forward_iterator` is `T`. 
*/ 
template <class Forward_iterator> void 
insert(Forward_iterator first, Forward_iterator beyond) ; 

/*! 

*/ 
handle find(handle p) ; 

/*! 
returns a 
canonical handle of the set that contains `p`, i.e., 
`P.same_set(p,q)` iff `P.find(p)` and `P.find(q)` 
return the same handle. 
\pre `p` is a handle in `P`. 
*/ 
const_handle find( const_handle p) ; 

/*! 
unites the sets of 
partition `P` containing \f$ p\f$ and \f$ q\f$. \pre \f$ p\f$ and \f$ q\f$ are in `P`. 
*/ 
void unify_sets( handle p, handle q); 

/*! 
returns 
true iff \f$ p\f$ and \f$ q\f$ belong to the same set of `P`. 
\pre \f$ p\f$ and \f$ q\f$ are in `P`. 
*/ 
bool same_set( const_handle p, const_handle q) ; 

/*! 
returns an iterator pointing to the 
first value of `P`. 
*/ 
iterator begin() ; 

/*! 
returns an iterator pointing beyond the 
last value of `P`. 
*/ 
iterator end() ; 

/// @}

}; /* end Union_find */
} /* end namespace CGAL */
