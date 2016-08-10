
/*!
\ingroup PkgAABB_treeConcepts
\cgalConcept

The concept `AABBTraits` provides the geometric primitive types and methods for the class `CGAL::AABB_tree<AABBTraits>`.

\cgalHasModel `CGAL::AABB_traits<AABBGeomTraits,AABBPrimitive>`

\cgalRefines `SearchGeomTraits_3`

\sa `CGAL::AABB_traits<AABBGeomTraits,AABBPrimitive>`
\sa `CGAL::AABB_tree<AABBTraits>`
\sa `AABBPrimitive`

*/
class AABBTraits {
public:

/// \name Types 
/// @{

/*!
Value type of the `Squared_distance` functor. 
*/ 
typedef unspecified_type FT; 

/*!
Type of a 3D point. 
*/ 
typedef unspecified_type Point_3; 

/*!
Type of primitive. 
Must be a model of the concepts `AABBPrimitive` or `AABBPrimitiveWithSharedData`.
*/ 
typedef unspecified_type Primitive; 

/*!
Bounding box type. 
*/ 
typedef unspecified_type Bounding_box; 

  /*!
  */ 
  enum Axis { 
    X_AXIS,
    Y_AXIS,
    Z_AXIS
  };

/*!

*/ 
typedef std::pair<Point_3, Primitive::Id> Point_and_primitive_id; 

/*!
\deprecated 
This requirement is deprecated and is no longer needed.
*/ 
typedef std::pair<Object, Primitive::Id> Object_and_primitive_id;


/*!
A nested class template providing as a pair the intersection result of a `Query` object
and a `Primitive::Datum`, together with the `Primitive::Id` of the primitive intersected.
The type of the pair is  `%Intersection_and_primitive_id<Query>::%Type`.
*/

template <typename Query>
using Intersection_and_primitive_id = unspecified_type;

/// @}

/// \name Splitting
/// During the construction of the AABB tree, the primitives are
/// sorted according to some comparison functions related to the \f$x\f$,
/// \f$ y\f$ or \f$ z\f$ coordinate axis:
/// @{

/*!
A functor object to split a range of primitives into two sub-ranges along the X-axis. Provides the operator: 
`void operator()(InputIterator first, InputIterator beyond);` %Iterator type `InputIterator` must be a model of RandomAccessIterator and have `Primitive` as value type. The operator is used for determining the primitives assigned to the two children nodes of a given node, assuming that the goal is to split the X-dimension of the bounding box of the node. The primitives assigned to this node are passed as argument to the operator. It should modify the iterator range in such a way that its first half and its second half correspond to the two children nodes. 
*/ 
typedef unspecified_type Split_primitives_along_x_axis; 

/*!
A functor object to split a range of primitives into two sub-ranges along the Y-axis. See `Split_primitives_along_x_axis` for the detailed description. 
*/ 
typedef unspecified_type Split_primitives_along_y_axis; 

/*!
A functor object to split a range of primitives into two sub-ranges along the Z-axis. See `Split_primitives_along_x_axis` for the detailed description. 
*/ 
typedef unspecified_type Split_primitives_along_z_axis; 

/*!
A functor object to compute the bounding box of a set of primitives. Provides the operator: 
`Bounding_box operator()(Input_iterator begin, Input_iterator beyond);` %Iterator type `InputIterator` must have `Primitive` as value type. 
*/ 
typedef unspecified_type Compute_bbox; 


// remove as not used any where in the code:
// A functor object to specify the direction along which the bounding box should be split: 
// `Axis operator()(const Bounding_box& bbox);` which returns the
// direction used for splitting `bbox`. It is usually the axis aligned
// with the longest edge of `bbox`.
// typedef unspecified_type Splitting_direction; 

/// @}


/// \name Intersections
/// The following predicates are required for each type `Query` for
/// which the class `CGAL::AABB_tree<AABBTraits>` may receive an intersection
/// detection or computation query.
/// @{

/*!
A functor object to compute intersection predicates between the query and the nodes of the tree. Provides the operators: 
- `bool operator()(const Query & q, const Bounding_box & box);` which returns `true` iff the query intersects the bounding box
- `bool operator()(const Query & q, const Primitive & primitive);` which returns `true` iff the query intersects the primitive
*/ 
typedef unspecified_type Do_intersect; 

/*!
A functor object to compute the intersection of a query and a primitive. Provides the operator: 
`boost::optional<Intersection_and_primitive_id<Query>::%Type > operator()(const Query & q, const Primitive& primitive);` which returns the intersection as a pair composed of an object and a primitive id, iff the query intersects the primitive. 

\cgalHeading{Note on Backward Compatibility}
Before the release 4.3 of \cgal, the return type of this function used to be `boost::optional<Object_and_primitive_id>`.
*/ 
typedef unspecified_type Intersect; 

/// \name Distance Queries
/// The following predicates are required for each
/// type `Query` for which the class `CGAL::AABB_tree<AABBTraits>` may receive a
/// distance query.
/// @{

/*!
A functor object to compute distance comparisons between the query and the nodes of the tree. Provides the operators: 
- `bool operator()(const Query & query, const Bounding_box& box, const Point & closest);` which returns `true` iff the bounding box is closer to `query` than `closest` is
- `bool operator()(const Query & query, const Primitive & primitive, const Point & closest);` which returns `true` iff `primitive` is closer to the `query` than `closest` is
*/ 
typedef unspecified_type Compare_distance; 

/*!
A functor object to compute closest point from the query on a primitive. Provides the operator: 
`Point_3 operator()(const Query& query, const Primitive& primitive, const Point_3 & closest);` which returns the closest point to `query`, among `closest` and all points of the primitive. 
*/ 
typedef unspecified_type Closest_point; 

/*!
A functor object to compute the squared distance between two points. Provides the operator: 
`FT operator()(const Point& query, const Point_3 & p);` which returns the squared distance between `p` and `q`. 
*/ 
typedef unspecified_type Squared_distance; 

/// @} 

/// \name Operations 
/// @{

/*!
Returns the primitive splitting functor for the X axis. 
*/ 
Split_primitives_along_x_axis split_primitives_along_x_axis_object(); 

/*!
Returns the primitive splitting functor for the Y axis. 
*/ 
Split_primitives_along_y_axis split_primitives_along_y_axis_object(); 

/*!
Returns the primitive splitting functor for the Z axis. 
*/ 
Split_primitives_along_z_axis split_primitives_along_z_axis_object(); 

/*!
Returns the bounding box constructor. 
*/ 
Compute_bbox compute_bbox_object(); 

/*!
Returns the intersection detection functor. 
*/ 
Do_intersect do_intersect_object(); 

/*!
Returns the intersection constructor. 
*/ 
Intersect intersect_object(); 

/*!
Returns the distance comparison functor. 
*/ 
Compare_distance compare_distance_object(); 

/*!
Returns the closest point constructor. 
*/ 
Closest_point closest_point_object(); 

/*!
Returns the squared distance functor.
*/ 
Squared_distance squared_distance_object(); 

/// @}

/// \name Primitive with Shared Data
/// In addition, if `Primitive` is a model of the concept `AABBPrimitiveWithSharedData`,
/// the following functions are part of the concept:
/// @{

/*!
the signature of that function must be the same as the static function `Primitive::construct_shared_data`.
The type `Primitive` expects that the data constructed by a call to `Primitive::construct_shared_data(t...)`
is the one given back when accessing the reference point and the datum of a primitive.
*/
template <class ... T>
void set_shared_data(T ... t);
{}

/*!
Returns the shared data of the primitive constructed after a call to `set_shared_data`.
If no call to `set_shared_data` has been done, `Primitive::Shared_data()` is returned.
*/
const Primitive::Shared_data& shared_data() const;

/// @}


}; /* end AABBTraits */

