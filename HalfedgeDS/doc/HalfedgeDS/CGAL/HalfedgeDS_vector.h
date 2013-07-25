namespace CGAL {

/*!
\ingroup PkgHDS_HDS

The class `HalfedgeDS_vector` is a model for the `HalfedgeDS` concept. 
`HalfedgeDS_vector` is a vector-based representation with random 
access iterators that does not support removal. 

\cgalModels `HalfedgeDS<Traits,Items,Alloc>`

\sa `CGAL::HalfedgeDS_default` 
\sa `CGAL::HalfedgeDS_list` 
\sa `HalfedgeDSItems` 
\sa `CGAL::Polyhedron_3<Traits>` 
\sa `CGAL::HalfedgeDS_items_decorator<HDS>` 
\sa `CGAL::HalfedgeDS_decorator<HDS>` 
\sa `CGAL::HalfedgeDS_const_decorator<HDS>` 

\cgalHeading{Implementation}

`HalfedgeDS_vector` uses internally the \stl `std::vector` container 
class. We require that we can create a `std::vector::iterator` 
from a pointer. If this will not be true any longer for any major \stl distribution we might switch to an internal implementation of a vector. 

The capacity is restricted to the reserved size. Allocations 
are not possible beyond the capacity without calling reserve again. 
All handles and iterators are invalidated upon a reserve call that 
increases the capacity. 

`CGAL_ALLOCATOR(int)` is used as default argument for the 
`Alloc` template parameter. 

*/
template< typename Traits, typename HalfedgeDSItems, typename Alloc >
class HalfedgeDS_vector {
public:

/// \name Types 
/// @{

/*!

*/ 
typedef random_access_iterator_tag iterator_category; 

/*!

*/ 
typedef CGAL::Tag_false Supports_removal; 

/// @}

}; /* end HalfedgeDS_vector */
} /* end namespace CGAL */
