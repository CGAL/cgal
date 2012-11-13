namespace CGAL {

/*!
\ingroup PkgIntervalSkipList

The class `Level_interval` represents intervals for the minimum and 
maximum value of the `z`-coordinate of a face of a triangulation. 

### Requirements ###

The `value_type` of `FaceHandle` must be `Face`, which must have a 
nested type `Vertex`, which must have a nested type `Point`, 
whose `Kernel_traits<Point>::Kernel` must have a nested type `FT`. 
These requirements are fulfilled, if one uses a <span class="textsc">Cgal</span> triangulation 
and a <span class="textsc">Cgal</span> kernel. 

\cgalModels ::Interval 

*/
template< typename FaceHandle >
class Level_interval {
public:

/// \name Types 
/// @{

/*! 
The type of the \f$ z\f$-coordinate of points stored in vertices of faces. 
*/ 
typedef FT Value; 

/// @} 

/// \name Creation 
/// @{

/*! 
%Default constructor. 
*/ 
Level_interval(); 

/*! 
Constructs the interval with smallest and largest `z` coordinate of the points 
stored in the vertices of the face `fh` points to. 
*/ 
Level_interval(FaceHandle fh); 

/// @} 

/// \name Operations 
/// @{

/*! 
Returns the face handle. 
*/ 
FaceHandle face_handle(); 

/// @}

}; /* end Level_interval */

/*! 
Inserts the interval `i` into the stream `os`. 
\pre The output operator for `*Face_handle` is defined. 
\relates Level_interval 
*/ 
template <typename FaceHandle>
ostream& operator<<(ostream& os, 
const Level_interval<FaceHandle>& i); 

} /* end namespace CGAL */
