
/*!
\ingroup PkgTetrahedralRemeshingConcepts
\cgalConcept

\cgalRefines TriangulationTraits_3

The concept `RemeshingTriangulationTraits_3` is the first template parameter
of the class `Remeshing_triangulation_3`. It defines the geometric objects
(points, segments, triangles, and tetrahedra) forming the triangulation together with a few
geometric predicates and constructions on these objects.

\cgalHasModel All models of `Kernel`.

\sa `CGAL::Triangulation_3`
*/

class RemeshingTriangulationTraits_3 {
public:

/// \name Types
/// @{

/*!
A constructor object model of `ConstructCrossProductVector_3`
*/
typedef unspecified_type Construct_cross_product_vector_3;

/*!
A constructor object model of `ConstructVector_3`
*/
typedef unspecified_type Construct_vector_3;

/*!
A constructor object model of `ConstructScaledVector_3 `
*/
typedef unspecified_type Construct_scaled_vector_3;

/*!
A constructor object model of `ConstructSumOfVectors_3`
*/
typedef unspecified_type Construct_sum_of_vectors_3;

/*!
A constructor object model of `ConstructOppositeVector_3`
*/
typedef unspecified_type Construct_opposite_vector_3;

/*!
A constructor object model of `ComputeSquaredLength_3`
*/
typedef unspecified_type Compute_squared_length_3;

/*!
A constructor object model of `ConstructDividedVector_3`
*/
typedef unspecified_type Construct_divided_vector_3;

/*!
A constructor object model of `ConstructTranslatedPoint_3`
*/
typedef unspecified_type Construct_translated_point_3;

/*!
A constructor object model of `ConstructMidpoint_3`
*/
typedef unspecified_type Construct_midpoint_3;

/*!
A constructor obeject model of `ComputeApproximateDihedralAngle_3`
*/
typedef unspecified_type Compute_approximate_dihedral_angle_3;

/////*!
////A predicate object that must provide the function operator
////
////`Comparison_result operator()(Point_3 p, Point_3 q)`,
////
////which returns `EQUAL` if the two points are equal. Otherwise it must
////return a consistent order for any two points chosen in a same line.
////*/
////typedef unspecified_type Compare_xyz_3;


/// @}


/// \name Operations
/// The following functions give access to the predicate and construction objects:
/// @{

/*!
*/
Construct_cross_product_vector_3 construct_cross_product_vector_3_object();

/*!
*/
Construct_vector_3 construct_vector_3_object();

/*!
*/
Construct_scaled_vector_3 construct_scaled_vector_3_object();

/*!
*/
Construct_sum_of_vectors_3 construct_sum_of_vectors_3_object();

/*!
*/
Construct_opposite_vector_3 construct_opposite_vector_3_object();

/*!
*/
Compute_squared_length_3 compute_squared_length_3_object();

/*!
*/
Construct_divided_vector_3 construct_divided_vector_3_object();

/*!
*/
Construct_translated_point_3 construct_translated_point_3_object();

/*!
*/
Construct_midpoint_3         construct_midpoint_3_object();

/*!
*/
Compute_approximate_dihedral_angle_3 compute_approximate_dihedral_angle_3_object();

/// @}

}; /* end RemeshingTriangulationTraits_3 */

