
namespace CGAL {

/*!
\ingroup MiscellanyRef

An instance `P` of the data type `Union_find<T,A>` is a
partition of values of type `T` into disjoint sets. The template
parameter `A` has to be a model of the allocator concept as defined
in the C++ standard. It has a default argument `CGAL_ALLOCATOR(T)`.

\cgalHeading{Implementation}

`Union_find<T,A>` is implemented with union by rank and path
compression. The running time for \f$ m\f$ set operations on \f$ n\f$ elements
is \cgalBigO{n \alpha(m,n)} where \f$ \alpha(m,n)\f$ is the extremely slow growing
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
typedef unspecified_type value_type;

/*!
handle to values.
*/
typedef unspecified_type handle;

/*!
iterator over values.
*/
typedef unspecified_type iterator;

/*!
allocator.
*/
typedef unspecified_type allocator;

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
the allocator of the partition.
*/
allocator get_allocator() const;

/*!
returns the number of disjoint sets of the partition.
*/
std::size_t number_of_sets() const;

/*!
returns the number of values of the partition.
*/
std::size_t size() const;

/*!
returns the memory consumed by the partition.
*/
std::size_t bytes() const;

/*!
returns the size of the set
containing `h`.
*/
std::size_t size( const_handle h) const;

/*!
reinitializes to an empty partition.
*/
void clear();

/*!
creates a new singleton set
containing `x` and returns a handle to it.
*/
handle make_set(const T& x);

/*!
same as `make_set()`.
*/
handle push_back(const T& x) ;

/*!
inserts the range of values referenced by `[first,beyond)`.
\tparam Forward_iterator must be a forward iterator with value type `T`.
*/
template <class Forward_iterator> void
insert(Forward_iterator first, Forward_iterator beyond) ;

/*!

*/
handle find(handle h) const;

/*!
returns a
canonical handle of the set that contains `h`, i.e.,
`P.same_set(h,h2)` iff `P.find(h) == P.find(h2)`.
\pre `h` is a handle in `P`.
*/
const_handle find( const_handle p) const;

/*!
unites the sets of
partition `P` containing `h1` and `h2`.
\pre `h1` and `h2` are in `P`.
*/
void unify_sets( handle h1, handle h2);

/*!
returns
true iff `h1` and `h2` belong to the same set of `P`.
\pre `h1` and `h2` are in `P`.
*/
bool same_set( const_handle h1, const_handle h2) const;

/*!
returns an iterator pointing to the
first value of the partition.
*/
iterator begin() const;

/*!
returns an iterator pointing beyond the
last value of the partition.
*/
iterator end() const;

/// @}

}; /* end Union_find */
} /* end namespace CGAL */
