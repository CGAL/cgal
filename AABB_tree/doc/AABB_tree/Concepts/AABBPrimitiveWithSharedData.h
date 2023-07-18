/*!
\ingroup PkgAABBTreeConcepts
\cgalConcept

The concept `AABBPrimitiveWithSharedData` describes the requirements for the primitives
stored in the AABB tree data structure. The concept encapsulates a type for the input
datum (a geometric object) and an identifier (id) type through which those primitives
are referred to. The concept `AABBPrimitiveWithSharedData` also refines the concepts
`DefaultConstructible` and `Assignable`.
The concept is similar to `AABBPrimitive` except that some data stored outside
of the primitives are required to access the datum and the reference point.

\sa `CGAL::AABB_tree<AABBTraits>`
\sa `AABBPrimitive`

\cgalHeading{Example}

The `Primitive` type can be a wrapper around an integer that refers to the position
of an object in a vector. Assume for instance that the input objects are some triangles.
The `Datum` would be a `Triangle_3` and the `Id` a `std::size_t`. The shared data here is a
`std::vector<Triangle_3>`.
The method `datum(const Shared_data&)` then returns a triangle from the vector.

\cgalHasModelsBegin
\cgalHasModels{CGAL::AABB_primitive<Id,ObjectPropertyMap,PointPropertyMap,Tag_true,CacheDatum>}
\cgalHasModels{CGAL::AABB_halfedge_graph_segment_primitive<HalfedgeGraph,VertexPointPMap,Tag_true,CacheDatum>}
\cgalHasModels{CGAL::AABB_face_graph_triangle_primitive<FaceGraph,VertexPointPMap,Tag_true,CacheDatum>}
\cgalHasModelsEnd
*/

class AABBPrimitiveWithSharedData {
public:

/// \name Types
/// @{
/*!
3D point type.
*/
typedef unspecified_type Point;

/*!
Type of input datum.
*/
typedef unspecified_type Datum;

/*!
Point reference type returned by the function `point(const Shared_data&)`. It is convertible to the type `Point`.
 */
typedef unspecified_type Point_reference;

/*!
Datum reference type returned by the function `datum(const Shared_data&)`. It is convertible to the type `Datum`.
*/
typedef unspecified_type Datum_reference;

/*!
Type of identifiers through which the input objects are referred to. It must be a model of the concepts DefaultConstructible and Assignable.
*/
typedef unspecified_type Id;

/*!
Type of the data shared amongst primitives
*/
typedef unspecified_type Shared_data;

/// @}

/// \name Operations
/// @{
/*!
returns the datum (geometric object) represented by the primitive.
*/
Datum_reference datum(const Shared_data& data);

/*!
returns the corresponding identifier. This identifier is only used as a reference for the objects in the output of the `AABB_tree` methods.
*/
Id id();

/*!
returns a 3D point located on the geometric object represented by the primitive. This function is used to sort the primitives during the AABB tree construction as well as to construct the search KD-tree internal to the AABB tree used to accelerate distance queries.
*/
Point_reference reference_point(const Shared_data& data);

/*!
constructs the shared data of a primitive.
The parameter pack is such that there exists a constructor `template <class T1, class ... T> AABBPrimitiveWithSharedData (T1,T...)`.
*/
template <class ... T>
static Shared_data construct_shared_data(T ... t);
/// @}

};
