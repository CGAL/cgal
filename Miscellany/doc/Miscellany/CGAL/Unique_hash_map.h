
namespace CGAL {

/*!
\ingroup PkgProfilingTools

An instance of the class template `Unique_hash_map` is an 
injective mapping from the set of keys of type `Key` to the set of 
variables of type `Data`. New keys can be inserted at any time, 
however keys cannot be individually deleted. 

An object `hash` of the type `UniqueHashFunction` returns a 
unique integer index `hash(key)` of type `std::size_t` for all 
objects \f$ key\f$ stored in `map`. The template parameter has as default 
the `Handle_hash_function` that hashes all types of pointers, handles, 
iterators, and circulators. 

The parameter `Allocator` has to match the standard allocator 
requirements, with value type `Data`. This parameter has the default 
value `CGAL_ALLOCATOR(Data)`. 

All variables are initialized to `default_data`, a value 
of type `Data` specified in the definition of `map`. 

\sa `UniqueHashFunction` 
\sa `CGAL::Handle_hash_function` 

\cgalHeading{Implementation}

`Unique_hash_map` is implemented via a chained hashing scheme. Access 
operations `map``[i]` take expected time \f$ O(1)\f$. The `table_size` 
parameter passed to chained hashing can be used to avoid unnecessary 
rehashing when set to the number of expected elements in the map. 
The design is derived from the \stl `hash_map` and the \leda type 
`map`. Its specialization on insertion only and unique hash values 
allow for a more time- and space-efficient implementation, see also 
\cgalCite{mn-lpcgc-00}, Chapter 5. This implementation makes also use 
of sentinels that lead to defined keys that have not been inserted. 

*/
template< typename Key, typename Data, typename UniqueHashFunction, typename Allocator>
class Unique_hash_map {
public:

/// \name Types 
/// In compliance with stl, the types `key_type`, `data_type`, and `hasher` are defined as well.
/// @{

/*!
the `Key` type. 
*/ 
typedef unspecified_type Key; 

/*!
the `Data` type. 
*/ 
typedef unspecified_type Data; 

/*!
the unique hash function type. 
*/ 
typedef unspecified_type Hash_function; 

/// @} 

/// \name Creation 
/// @{

/*!

creates an injective function from `Key` to the set of unused 
variables of type `Data`, sets `default_data` to `default`, 
passes the `table_size` as argument to the internal implementation, 
and initializes the hash function with `fct`. 
*/ 
Unique_hash_map<Key,Data,UniqueHashFunction>( 
const Data& default = Data(), 
std::size_t table_size = 1, 
const Hash_function& fct = Hash_function()); 

/*!

creates an injective function from `Key` to the set of unused 
variables of type `Data`, sets `default_data` to `default`, 
passes the `table_size` as argument to the internal implementation, 
initializes the hash function with `fct`, and inserts all keys 
from the range `[first1,beyond1)`. The data variable for each 
inserted `key` is initialized with the corresponding value from 
the range `[first2, first2 + (beyond1-first1))`. 
\pre The increment operator must be defined for values 
of type `Key` and for values of type `Data`. `beyond1` 
must be reachable from `first1` using increments. 
*/ 
Unique_hash_map<Key,Data,UniqueHashFunction>( 
Key first1, Key beyond1, Data first2, 
const Data& default = Data(), 
std::size_t table_size = 1, 
const Hash_function& fct = Hash_function()); 

/// @} 

/// \name Operations 
/// @{

/*!
the current `default_value`. 
*/ 
Data default_value() const; 

/*!
the current hash function. 
*/ 
Hash_function hash_function() const; 

/*!
returns true if \f$ key\f$ is 
defined in `*this`. Note that there can be keys defined that have not 
been inserted explicitly. Their variables are initialized to 
`default_value`. 
*/ 
bool is_defined( const Key& key) const; 

/*!

resets `*this` to the injective function from `Key` to the 
set of unused variables of type `Data`. The `default_data` 
remains unchanged. 
*/ 
void clear(); 

/*!

resets `*this` to the injective function from `Key` to the 
set of unused variables of type `Data` and sets `default_data` 
to `default`. 
*/ 
void clear(const Data& default); 

/*!

returns a reference to the variable `map``(key)`. If `key` 
has not been inserted into `map` before, `key` is inserted and 
initialized with `default_value`. 
*/ 
Data& operator[](const Key& key); 

/*!

returns a const reference to the variable `*this``(key)`. If `key` 
has not been inserted into `*this` before, a const reference to the 
`default_value` is returned. However, `key` is not inserted 
into `*this`. 
*/ 
const Data& operator[](const Key& key) const; 

/*!

inserts all keys from the range `[first1,beyond1)`. 
The data variable for each inserted `key` is initilized with the 
corresponding value from the range `[first2, first2 + 
(beyond1-first1))`. Returns `first2 + (beyond1-first1)`. 
\pre The increment operator must be defined for values 
of type `Key` and for values of type `Data`. `beyond1` 
must be reachable from `first1` using increments. 
*/ 
Data insert( Key first1, Key beyond1, Data first2); 

/// @}

}; /* end Unique_hash_map */
} /* end namespace CGAL */
