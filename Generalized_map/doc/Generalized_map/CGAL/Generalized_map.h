
namespace CGAL {

/*!
\ingroup PkgGeneralizedMapsClasses

The class `Generalized_map` represents a <I>d</I>D generalized map.

Two versions exist: one where darts and non void attributes are stored in memory using `Compact_container`, using `Alloc` as allocator, and use handles as descriptors; a second one where darts and non void attributes are stored in an internal `std::vector` like data-structure, and use indices as descriptors. The choice between the two versions is done through the item class.

\cgalModels{GeneralizedMap}

\tparam d the dimension of the map.

\tparam Items a model of the `GenericMapItems` concept. Equal to `Generic_map_min_items` by default.

\tparam Alloc has to match the standard allocator requirements. The `rebind` mechanism  `Alloc` will be used to create appropriate allocators internally with value type `Dart`. Equal to `CGAL_ALLOCATOR(int)` from `<CGAL/memory.h>` by default.

\cgalHeading{Complexity}

The complexity of \link GeneralizedMap::sew `sew`\endlink and \link GeneralizedMap::unsew `unsew`\endlink is in <I>O</I>( \f$ | \f$ <I>S</I> \f$ | \f$ \f$ \times \f$ \f$ | \f$ <I>c</I> \f$ | \f$ ), <I>S</I> being the set of darts of the orbit \f$ \langle{}\f$\f$ \alpha_0\f$,\f$ \ldots\f$,\f$ \alpha_{i-2}\f$,\f$ \alpha_{i+2}\f$,\f$ \ldots\f$,\f$ \alpha_d\f$\f$ \rangle{}\f$ for the considered dart, and <I>c</I> the biggest <I>j</I>-cell merged or split during the sew such that <I>j</I>-attributes are non void.
The complexity of \link GeneralizedMap::is_sewable `is_sewable`\endlink is in <I>O</I>( \f$ | \f$ <I>S</I> \f$ | \f$ ).

The complexity of \link GenericMap::set_attribute `set_attribute`\endlink is in <I>O</I>( \f$ | \f$ <I>c</I> \f$ | \f$ ), <I>c</I> being the <I>i</I>-cell containing the considered dart.

The complexity of \link GenericMap::is_without_boundary(unsigned int i) const `is_without_boundary(i)`\endlink is in <I>O</I>( \f$ | \f$ <I>D</I> \f$ | \f$ ), <I>D</I> being the set of darts of the generalized map, and the complexity of \link GenericMap::is_without_boundary() const `is_without_boundary()`\endlink is in <I>O</I>( \f$ | \f$ <I>D</I> \f$ | \f$ \f$ \times \f$ <I>d</I> ).

The complexity of \link GenericMap::unmark_all `unmark_all`\endlink and \link GenericMap::free_mark `free_mark`\endlink is in <I>O</I>( 1 ) if all the darts of the generalized map have the same mark, and in <I>O</I>( \f$ | \f$ <I>D</I> \f$ | \f$ ) otherwise.

The complexity of  \link GeneralizedMap::is_valid `is_valid`\endlink is in <I>O</I>( \f$ | \f$ <I>D</I> \f$ | \f$ \f$ \times \f$ <I>d</I>\f$ ^2\f$ ).

The complexity of \link GenericMap::clear `clear`\endlink is in <I>O</I>( \f$ | \f$ <I>D</I> \f$ | \f$ \f$ \times \f$ <I>d</I> ).

Other methods have all a constant time complexity.

\sa `GenericMapItems`

*/
template< unsigned int d, typename Items, typename Alloc >
class Generalized_map {
public:

/// \name Constants
/// @{

/*!
The dimension of the generalized map.
*/
static const unsigned int dimension = d;

/// @}

/// \name Types
/// @{

/*!

*/
typedef Generalized_map<d,Items,Alloc> Self;

/*!
Information associated with darts. Equal to `void` if `Dart_info` is not defined in the items class.
*/
typedef Items::Dart_wrapper<Self>::Dart_info Dart_info;

/*!
The tuple of cell attributes. Equal to `std::tuple<>` if `Attributes` is not defined in the items class.
*/
typedef Items::Dart_wrapper<Self>::Attributes Attributes;

/// @}

}; /* end Generalized_map */
} /* end namespace CGAL */
