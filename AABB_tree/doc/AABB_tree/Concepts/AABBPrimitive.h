
/*!
\ingroup PkgAABB_treeConcepts
\cgalconcept

The concept `AABBPrimitive` describes the requirements for the primitives stored in the AABB tree data structure. The concept encapsulates a type for the input datum (a geometric object) and an identifier (id) type through which those primitives are referred to. The concept `AABBPrimitive` also refines the concepts DefaultConstructible and Assignable. 

\sa `AABB_tree<AT>` 

### Example ###

The `Primitive` type can be, e.g., a wrapper around a `Handle`. Assume for instance that the input objects are the triangle faces of a mesh stored as a `CGAL::Polyhedron`. The `Datum` would be a `Triangle_3` and the `Id` would be a polyhedron `Face_handle`. Method `datum()` can return either a `Triangle_3` constructed on the fly from the face handle or a `Triangle_3` stored internally. This provides a way for the user to trade memory for efficiency. 

\hasModel `CGAL::AABB_polyhedron_triangle_primitive<GeomTraits,Polyhedron>`
\hasModel `CGAL::AABB_polyhedron_segment_primitive<GeomTraits,Polyhedron>` 

*/

class AABBPrimitive {
public:

/// \name Types 
/// @{

/*! 
3D point type. 
*/ 
typedef Hidden_type Point; 

/*! 
Type of input datum. 
*/ 
typedef Hidden_type Datum; 

/*! 
Type of identifiers through which the input objects are referred to. It must be a model of the concepts DefaultConstructible and Assignable. 
*/ 
typedef Hidden_type Id; 

/// @} 

/// \name Operations 
/// @{

/*! 
Returns the datum (geometric object) represented by the primitive. 
*/ 
Datum datum(); 

/*! 
Returns the corresponding identifier. This identifier is only used as a reference for the objects in the output of the `AABB_tree` methods. 
*/ 
Id id(); 

/*! 
Returns a 3D point located on the geometric object represented by the primitive. This function is used to sort the primitives during the AABB tree construction as well as to construct the search KD-tree internal to the AABB tree used to accelerate distance queries. 
*/ 
Point reference_point(); 

/// @}

}; /* end AABBPrimitive */

