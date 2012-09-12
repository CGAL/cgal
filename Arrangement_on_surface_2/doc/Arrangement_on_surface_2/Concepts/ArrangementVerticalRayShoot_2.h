
/*!
\ingroup PkgArrangement2Concepts
\cgalconcept

A model of the `ArrangementVerticalRayShoot_2` concept can be attached to an `Arrangement_2` 
instance and answer vertical ray-shooting queries on this arrangement. 
Namely, given a `Arrangement_2::Point_2` object, representing a point in 
the plane, it returns the arrangement feature (edge or vertex) that lies 
strictly above it (or below it). By "strictly" we mean that if the 
query point lies on an arrangement edge (or on an arrangement vertex) this 
edge will <I>not</I> be the query result, but the feature lying above or 
below it. (An exception to this rule is the degenerate situation where the 
query point lies in the interior of a vertical edge.) Note that it may happen 
that the query point lies above the upper envelope (or below the lower 
envelope) of the arrangement, so that the vertical ray emanating from it 
may go to infinity without hitting any arrangement feature on its way. In this 
case the unbounded face is returned. 

\hasModel Arr_naive_point_location<Arrangement> 
\hasModel Arr_walk_along_a_line_point_location<Arrangement> 
\hasModel Arr_trapezoid_ric_point_location<Arrangement> 
\hasModel Arr_landmarks_point_location<Arrangement,Generator> 

*/

class ArrangementVerticalRayShoot_2 {
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
ArrangementVerticalRayShoot_2(); 

/*! 
constructs a ray-shooting object `rs` attached to the given 
arrangement `arr`. 
*/ 
ArrangementVerticalRayShoot_2 (const Arrangement_2& arr); 

/// @} 

/// \name Query Functions 
/// @{

/*! 
locates the arrangement feature that is first hit by an upward-directed 
vertical ray emanating from the query point `q`, 
and returns a handle for this feature. 
The function returns an `Object` instance that is a wrapper for 
one of the following types: 
<UL> 
<LI>`Arrangement_2::Halfedge_const_handle`, in case the vertical 
ray hits an arrangement edge; 
<LI>`Arrangement_2::Vertex_const_handle`, in case the vertical 
ray hits an arrangement vertex. 
<LI>`Arrangement_2::Face_const_handle` for the unbounded arrangement 
face, in case `q` lies above the upper envelope of the 
arrangement. 
</UL> 
\pre `rs` is attached to a valid arrangement instance. 
*/ 
Object ray_shoot_up (const Point_2& q) const; 

/*! 
locates the arrangement feature that is first hit by a downward-directed 
vertical ray emanating from the query point `q`, 
and returns a handle for this feature. 
The function returns an `Object` instance that is a wrapper for 
one of the following types: 
<UL> 
<LI>`Arrangement_2::Halfedge_const_handle`, in case the vertical 
ray hits an arrangement edge; 
<LI>`Arrangement_2::Vertex_const_handle`, in case the vertical 
ray hits an arrangement vertex. 
<LI>`Arrangement_2::Face_const_handle` for the unbounded arrangement 
face, in case `q` lies below the lower envelope of the 
arrangement. 
</UL> 
\pre `rs` is attached to a valid arrangement instance. 
*/ 
Object ray_shoot_down (const Point_2& q) const; 

/// @} 

/// \name Operations 
/// @{

/*! 
attaches `rs` to the given arrangement `arr`. 
*/ 
void attach (const Arrangement_2& arr); 

/*! 
detaches `rs` from the arrangement it is currently attached to. 
*/ 
void detach (); 

/// @}

}; /* end ArrangementVerticalRayShoot_2 */

