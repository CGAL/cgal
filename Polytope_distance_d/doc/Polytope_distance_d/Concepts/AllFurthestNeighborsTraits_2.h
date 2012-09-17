
/*!
\ingroup PkgOptimalDistancesConcepts
\cgalconcept

The concept `AllFurthestNeighborsTraits_2` defines types and operations 
needed to compute all furthest neighbors for the vertices of a 
convex polygon using the function `all_furthest_neighbors_2`. 

\hasModel `CGAL::Cartesian<FieldNumberType>`
\hasModel `CGAL::Homogeneous<RingNumberType>`
\hasModel `CGAL::Simple_cartesian<FieldNumberType>`
\hasModel `CGAL::Simple_homogeneous<RingNumberType>`

\sa `CGAL::all_furthest_neighbors_2` 

### Notes ###

<UL> 
<LI>`AllFurthestNeighborsTraits_2``::Less_xy_2` and 
`AllFurthestNeighborsTraits_2``::Orientation_2` are used for (expensive) 
precondition checking only. Therefore, they need not to be 
specified, in case that precondition checking is disabled. 
</UL> 

*/

class AllFurthestNeighborsTraits_2 {
public:

/// \name Types 
/// @{

/*! 
model for `FieldNumberType`. 
*/ 
typedef Hidden_type FT; 

/*! 
model for 
`Kernel::Point_2`. 
*/ 
typedef Hidden_type Point_2; 

/*! 
model for 
`Kernel::Compute_squared_distance_2`. 
*/ 
typedef Hidden_type Compute_squared_distance_2; 

/*! 
model for 
`Kernel::Less_xy_2`. 
*/ 
typedef Hidden_type Less_xy_2; 

/*! 
model for 
`Kernel::Orientation_2`. 
*/ 
typedef Hidden_type Orientation_2; 

/// @} 

/// \name Operations 
/// The following member functions return function objects of the types listed above.
/// @{

/*! 

*/ 
Compute_squared_distance_2 
compute_squared_distance_2_object(); 

/*! 

*/ 
Less_xy_2 less_xy_2_object(); 

/*! 

*/ 
Orientation_2 orientation_2_object(); 

/// @}

}; /* end AllFurthestNeighborsTraits_2 */

