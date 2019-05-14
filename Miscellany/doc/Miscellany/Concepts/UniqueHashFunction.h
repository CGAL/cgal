
/*!
\ingroup PkgProfilingToolsConcepts
\cgalConcept

`UniqueHashFunction` is a concept for a hash function with unique hash values. 
An instance `hash` for a model of the `UniqueHashFunction` concept is a 
function object. It maps objects of its domain type `Key` to 
the integral image type `std::size_t`. The image values have to 
be unique for all keys in the domain type `Key`. 

\cgalHasModel `CGAL::Handle_hash_function` 

\sa `CGAL::Unique_hash_map<Key,Data,UniqueHashFunction>` 

*/

class UniqueHashFunction {
public:

/// \name Types 
/// @{

/*!
type of the hash value. 
*/ 
typedef std::size_t result_type; 

/// @} 

/// \name Creation 
/// @{

/*!

copy constructor. 
*/ 
UniqueHashFunction( const UniqueHashFunction& hash2); 

/*!
assignment. 
*/ 
UniqueHashFunction& operator=(const UniqueHashFunction& hash2); 

/// @} 

/// \name Operations 
/// @{

/*!

returns a unique hash value for the `key` value. 
*/ 
std::size_t operator()( const Key& key); 

/// @}

}; /* end UniqueHashFunction */

