
/*!
\ingroup PkgBoundingVolumesConcepts
\cgalConcept

The concept `RectangularPCenterTraits_2` defines types and operations 
needed to compute rectilinear \f$ p\f$-centers of a planar point set 
using the function `CGAL::rectangular_p_center_2`. 

\cgalHasModel `CGAL::Rectangular_p_center_default_traits_2<K>` 

\sa `CGAL::rectangular_p_center_2()` 

*/

class RectangularPCenterTraits_2 {
public:

/// \name Types 
/// @{

/*!
model for `FieldNumberType`. 
*/ 
typedef unspecified_type FT; 

/*!
model for 
`Kernel::Point_2`. 
*/ 
typedef unspecified_type Point_2; 

/*!
model for 
`Kernel::Iso_rectangle_2`. 
*/ 
typedef unspecified_type Iso_rectangle_2; 

/*!
model for 
`Kernel::Less_x_2`. 
*/ 
typedef unspecified_type Less_x_2; 

/*!
model for 
`Kernel::Less_y_2`. 
*/ 
typedef unspecified_type Less_y_2; 

/*!
model for 
`Kernel::Construct_vertex_2`. 
*/ 
typedef unspecified_type Construct_vertex_2; 

/*!
model for 
`Kernel::Construct_iso_rectangle_2`. 
*/ 
typedef unspecified_type Construct_iso_rectangle_2; 

/*!
adaptable binary function 
class: `Point_2` \f$ \times\f$ `Point_2` \f$ \rightarrow\f$ 
`FT` returns the signed distance of two points' 
\f$ x\f$-coordinates. 
*/ 
typedef unspecified_type Signed_x_distance_2; 

/*!
adaptable binary function 
class: `Point_2` \f$ \times\f$ `Point_2` \f$ \rightarrow\f$ 
`FT` returns the signed distance of two points' 
\f$ y\f$-coordinates. 
*/ 
typedef unspecified_type Signed_y_distance_2; 

/*!
adaptable binary function 
class: `Point_2` \f$ \times\f$ `Point_2` \f$ \rightarrow\f$ 
`FT` returns the \f$ ||\cdot||_{\infty}\f$ distance of two 
points. 
*/ 
typedef unspecified_type Infinity_distance_2; 

/*!
adaptable binary 
function class: `Point_2` \f$ \times\f$ `Point_2` 
\f$ \rightarrow\f$ `FT` returns the signed \f$ ||\cdot||_{\infty}\f$ 
distance of two points. 
*/ 
typedef unspecified_type Signed_infinity_distance_2; 

/*!

3-argument function class: `Point_2` \f$ \times\f$ `Point_2` 
\f$ \times\f$ `FT` \f$ \rightarrow\f$ `Point_2`. For arguments 
\f$ (p,\,q,\,r)\f$ it returns the lower-left corner of the iso-oriented 
square with sidelength \f$ r\f$ and upper-right corner at the 
intersection of the vertical line through \f$ p\f$ and the horizontal 
line through \f$ q\f$. 
*/ 
typedef unspecified_type Construct_point_2_below_left_implicit_point_2; 

/*!

3-argument function class: `Point_2` \f$ \times\f$ `Point_2` 
\f$ \times\f$ `FT` \f$ \rightarrow\f$ `Point_2`. For arguments 
\f$ (p,\,q,\,r)\f$ it returns the lower-right corner of the 
iso-oriented square with sidelength \f$ r\f$ and upper-left corner at 
the intersection of the vertical line through \f$ p\f$ and the 
horizontal line through \f$ q\f$. 
*/ 
typedef unspecified_type Construct_point_2_below_right_implicit_point_2; 

/*!

3-argument function class: `Point_2` \f$ \times\f$ `Point_2` 
\f$ \times\f$ `FT` \f$ \rightarrow\f$ `Point_2`. For arguments 
\f$ (p,\,q,\,r)\f$ it returns the upper-right corner of the 
iso-oriented square with sidelength \f$ r\f$ and lower-left corner at 
the intersection of the vertical line through \f$ p\f$ and the 
horizontal line through \f$ q\f$. 
*/ 
typedef unspecified_type Construct_point_2_above_right_implicit_point_2; 

/*!

3-argument function class: `Point_2` \f$ \times\f$ `Point_2` 
\f$ \times\f$ `FT` \f$ \rightarrow\f$ `Point_2`. For arguments 
\f$ (p,\,q,\,r)\f$ it returns the upper-left corner of the iso-oriented 
square with sidelength \f$ r\f$ and lower-right corner at the 
intersection of the vertical line through \f$ p\f$ and the horizontal 
line through \f$ q\f$. 
*/ 
typedef unspecified_type Construct_point_2_above_left_implicit_point_2; 

/// @} 

/// \name Operations 
/// For every function class listed above there is a member function
/// to fetch the corresponding function object.
/// @{

/*!

*/ 
Inf_distance_2 inf_distance_2_object() const; 

/*!

*/ 
Signed_inf_distance_2 
signed_inf_distance_2_object() const; 

/*!

*/ 
Construct_vertex_2 
construct_vertex_2_object() const; 

/*!

*/ 
Construct_iso_rectangle_2 
construct_iso_rectangle_2_object() const; 

/*!

*/ 
Construct_iso_rectangle_2_below_left_point_2 
construct_iso_rectangle_2_below_left_point_2_object() const; 

/*!

*/ 
Construct_iso_rectangle_2_above_left_point_2 
construct_iso_rectangle_2_above_left_point_2_object() const; 

/*!

*/ 
Construct_iso_rectangle_2_below_right_point_2 
construct_iso_rectangle_2_below_right_point_2_object() const; 

/*!

*/ 
Construct_iso_rectangle_2_above_right_point_2 
construct_iso_rectangle_2_above_right_point_2_object() const; 

/// @}

}; /* end RectangularPCenterTraits_2 */

