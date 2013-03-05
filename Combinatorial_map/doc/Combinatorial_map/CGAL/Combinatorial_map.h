
namespace CGAL {

/*!
\ingroup PkgCombinatorialMapsClasses

The class `Combinatorial_map` represents a <I>d</I>D combinatorial map.

Darts and non void attributes are stored in memory using
`Compact_container`, using `Alloc` as allocator.

\cgalModels `CombinatorialMap`


\tparam d is an integer for the dimension of the map.

\tparam Items must be a model of the `CombinatorialMapItems` concept.

\tparam Alloc has to match the standard allocator requirements. The `rebind` mechanism  `Alloc` will be used to create appropriate allocators internally with value type `Dart`.

There are two default template arguments:
`Combinatorial_map_min_items<d>` for `Items` and
`CGAL_ALLOCATOR(int)` from the `<CGAL/memory.h>` header file for `Alloc`.

\cgalHeading{Complexity}


The complexity of `sew` and `unsew` is in <I>O</I>(\f$ |\f$<I>S</I>\f$ |\f$\f$ \times\f$\f$ |\f$<I>c</I>\f$ |\f$), <I>S</I>
being the set of darts of the orbit
\f$ \langle{}\f$\f$ \beta_1\f$,\f$ \ldots\f$,\f$ \beta_{i-2}\f$,\f$ \beta_{i+2}\f$,\f$ \ldots\f$,\f$ \beta_d\f$\f$ \rangle{}\f$ for the
considered dart, and <I>c</I> the biggest <I>j</I>-cell merged or
split during the sew such that <I>j</I>-attributes are non void.
The complexity of `is_sewable` is in <I>O</I>(\f$ |\f$<I>S</I>\f$ |\f$).

The complexity of `set_attribute` is in <I>O</I>(\f$ |\f$<I>c</I>\f$ |\f$), <I>c</I> being the
<I>i</I>-cell containing the considered dart.

The complexity of `is_without_boundary(unsigned int i)` is in <I>O</I>(\f$ |\f$<I>D</I>\f$ |\f$),
<I>D</I> being the set of darts of the combinatorial map, and the complexity of
`is_without_boundary()` is in <I>O</I>(\f$ |\f$<I>D</I>\f$ |\f$\f$ \times\f$<I>d</I>).

The complexity of `unmark_all` and `free_mark` is in <I>O</I>(1) if
all the darts of the combinatorial map have the same mark, and in
<I>O</I>(\f$ |\f$<I>D</I>\f$ |\f$) otherwise.

The complexity of `is_valid` is in
<I>O</I>(\f$ |\f$<I>D</I>\f$ |\f$\f$ \times\f$<I>d</I>\f$ ^2\f$).

The complexity of `clear` is in <I>O</I>(\f$ |\f$<I>D</I>\f$ |\f$\f$ \times\f$<I>d</I>).

Other methods have all a constant time complexity.

\sa `CombinatorialMapItems`
\sa `Dart`

*/
template< unsigned int d, typename Items, typename Alloc >
class Combinatorial_map {
public:

/// \name Types
/// @{

/*!

*/
typedef Combinatorial_map<d,Items,Alloc> Self;

/*!

*/
typedef Items::Dart_wrapper<Self>::Dart Dart;

/// @}

}; /* end Combinatorial_map */
} /* end namespace CGAL */
