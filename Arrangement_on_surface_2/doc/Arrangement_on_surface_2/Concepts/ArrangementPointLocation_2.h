
/*!
\ingroup PkgArrangement2Concepts
\cgalconcept

A model of the `ArrangementPointLocation_2` concept can be attached to an `Arrangement_2` 
instance and answer point-location queries on this arrangement. Namely, given 
a `Arrangement_2::Point_2` object, representing a point in the plane, 
it returns the arrangement cell containing it. In the general case, the 
query point is contained inside an arrangement face, but in degenerate 
situations it may lie on an edge or coincide with an arrangement vertex. 

\hasModel Arr_naive_point_location<Arrangement> 
\hasModel Arr_walk_along_a_line_point_location<Arrangement> 
\hasModel Arr_trapezoid_ric_point_location<Arrangement> 
\hasModel Arr_landmarks_point_location<Arrangement,Generator> 

*/

class ArrangementPointLocation_2 {
public:

/// \name Types 
/// @{

/*! 
the associated arrangement type. 
*/ 
typedef Hidden_type Arrangement_2; 

/*! 
equivalent to `Arrangement_2::Point_2`. 
*/ 
typedef Hidden_type Point_2; 

/// @} 

/// \name Creation 
/// @{

/*! 
default constructor. 
*/ 
ArrangementPointLocation_2(); 

/*! 
constructs a point-location object `pl` attached to the given 
arrangement `arr`. 
*/ 
ArrangementPointLocation_2 (const Arrangement_2& arr); 

/// @} 

/// \name Query Functions 
/// @{

/*! 
locates the arrangement cell that contains the query point `q` 
and returns a handle for this cell. 
The function returns an `Object` instance that wraps either of the 
following types: 
<UL> 
<LI>`Arrangement_2::Face_const_handle`, in case `q` is 
contained inside an arrangement face; 
<LI>`Arrangement_2::Halfedge_const_handle`, in case `q` lies 
on an arrangement edge; 
<LI>`Arrangement_2::Vertex_const_handle`, in case `q` coincides 
with an arrangement vertex. 
</UL> 
\pre `pl` is attached to a valid arrangement instance. 
*/ 
Object locate (const Point_2& q) const; 

/// @} 

/// \name Operations 
/// @{

/*! 
attaches `pl` to the given arrangement `arr`. 
*/ 
void attach (const Arrangement_2& arr); 

/*! 
detaches `pl` from the arrangement it is currently attached to. 
*/ 
void detach (); 

/// @}

}; /* end ArrangementPointLocation_2 */

