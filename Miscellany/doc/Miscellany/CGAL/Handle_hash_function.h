
namespace CGAL {

/*!
\ingroup MiscellanyRef

The class `Handle_hash_function` is a model for the `UniqueHashFunction`
concept. It is applicable for all key types with pointer-like
functionality, such as handles, iterators, and circulators.
Specifically, for a `key` value the expression `&*key` must
return a unique address.

\cgalModels `UniqueHashFunction`

\sa `CGAL::Unique_hash_map<Key,Data,UniqueHashFunction>`

\cgalHeading{Implementation}

Plain type cast of `&*key` to `std::size_t` and devided
by the size of the `std::iterator_traits<Handle>::%value_type` to
avoid correlations with the internal table size, which is a power of
two.

*/

struct Handle_hash_function {
public:

/// \name Creation
/// @{

/*!
%Default constructor.
*/
Handle_hash_function();

/// @}

/// \name Operations
/// @{

/*!

Returns unique hash value for any `Handle`
type for which `&*key` gives a unique address.

The type `std::iterator_traits<Handle>::%value_type` has to be defined
(which it is already for pointers, handles, iterators, and
circulators).
*/
template <class Handle>
std::size_t operator()( const Handle& key);

/// @}

}; /* end Handle_hash_function */
} /* end namespace CGAL */
