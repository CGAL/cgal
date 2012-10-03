


/*!
\ingroup PkgMesh2Concepts
\cgalconcept


The concept `ConformingDelaunayTriangulationTraits_2` refines the concept 
`ConstrainedDelaunayTriangulationTraits_2` by providing a numeric 
field type `FT`, a type `Vector_2` and several constructors on 
`Vector_2`, `Point_2`, and a predicate on angles. 
The field type has to be a model of 
the concept `::SqrtFieldNumberType`. This field type and the 
constructors are used by the conforming algorithm to compute Steiner 
points on constrained edges. 

\refines ::DelaunayTriangulationTraits_2 

\hasModel Any model of `Kernel` concept. In particular, all \cgal kernels. 


*/

class ConformingDelaunayTriangulationTraits_2 {
public:


/// \name Types 
/// @{

/*! 
The field type. It must be a model of 
`SqrtFieldNumberType`, that is must be a number type 
supporting the operations \f$ +\f$, \f$ -\f$, \f$ *\f$, \f$ /\f$, and \f$ \sqrt{\cdot}\f$. 
*/ 
typedef Hidden_type FT; 


/*! 
The vector type. 
*/ 
typedef Hidden_type Vector_2; 


/*! 
Constructor object. Must 
provide the operator `Vector_2 operator()(Point a, Point b)` 
that computes the vector \f$ b-a\f$. 
*/ 
typedef Hidden_type Construct_vector_2; 


/*! 
Constructor object. Must 
provide the operator `Vector_2 operator()(Vector_2 v, FT scale)` 
that computes the vector \f$ scale \cdot\mathbf{v}\f$. 
*/ 
typedef Hidden_type Construct_scaled_vector_2; 


/*! 
Constructor object. Must 
provide the operator `Point_2 operator()(Point_2 p, Vector_2 v)` 
that computes the point \f$ p + \mathbf{v}\f$. 
*/ 
typedef Hidden_type Construct_translated_point_2; 


/*! 
Constructor object. Must provide 
the operator `Point_2 operator()(Point_2 a, Point_2 b)` that 
computes the midpoint of the segment \f$ ab\f$. 
*/ 
typedef Hidden_type Construct_midpoint_2; 


/*! 
Constructor object. Must 
provide the operator `FT operator()(Point_2 a, Point_2 b)` that 
computes the squared distance between \f$ a\f$ and \f$ b\f$. 
*/ 
typedef Hidden_type Compute_squared_distance_2; 


/*! 
Predicate object. Must provide the operator 
`CGAL::Angle operator()(Point_2 p, Point_2 q, Point_2 r)` that 
returns OBTUSE, RIGHT or ACUTE depending on the angle formed by the three 
points p, q, r (q being the vertex of the angle). 
*/ 
typedef Hidden_type Angle_2; 

/// @} 


/// \name Access to predicate and constructor objects 
/// @{

/*! 

*/ 
Construct_vector_2 construct_vector_2_object(); 





/*! 

*/ 
Construct_scaled_vector_2 construct_scaled_vector_2_object(); 





/*! 

*/ 
Construct_translated_point_2 
construct_translated_point_2_object(); 





/*! 

*/ 
Construct_midpoint_2 construct_midpoint_2_object(); 





/*! 

*/ 
Compute_squared_distance_2 
compute_squared_distance_2_object(); 





/*! 

*/ 
Angle_2 angle_2_object(); 





/// @}

}; /* end ConformingDelaunayTriangulationTraits_2 */

