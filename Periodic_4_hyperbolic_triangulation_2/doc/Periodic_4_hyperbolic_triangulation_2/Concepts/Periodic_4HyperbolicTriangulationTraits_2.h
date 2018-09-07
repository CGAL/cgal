// Copyright (c) 2016-2017 INRIA Nancy - Grand Est (France).
// All rights reserved.

/*!
\ingroup PkgPeriodic4HyperbolicTriangulation2Concepts
\cgalConcept

\cgalRefines HyperbolicDelaunayTriangulationTraits_2

The concept `Periodic_4HyperbolicTriangulationTraits_2` refines the concept 
`HyperbolicDelaunayTriangulationTraits_2`. It describes the set of requirements to be 
fulfilled by any class used to instantiate the first template parameter of the class 
`CGAL::Periodic_4_hyperbolic_triangulation_2`. In addition to the geometric 
types and the operations defined on them in `HyperbolicDelaunayTriangulationTraits_2`, 
it defines the hyperbolic translations that allow to encode the periodicity of the 
triangulation.

\cgalModifBegin
	The concept requires that two nested types are provided:
	<ul>
		<li> A field number type `FT` that must support exact computations with algebraic 
			 numbers, in particular with nested square roots;
		<li> A `Hyperbolic_translation` type, which satisfies the requirements described 
			 in the concept `HyperbolicOctagonTranslation`.
	</ul>
\cgalModifEnd

*/



class Periodic_4HyperbolicTriangulationTraits_2 {

public:


	/// \name Types
	/// @{

		/*!
			\cgalModifBegin
			Represents a hyperbolic circle, i.e., a circle contained in the unit disk.
			\cgalModifEnd
		*/
		typedef unspecified_type 				Circle_2;

		/*!
			Represents a hyperbolic translation. 
			Must be a model of the concept `HyperbolicOctagonTranslation`. 
		*/
		typedef unspecified_type 				Hyperbolic_translation;
	/// @}


	/// \name Predicate Types
	/// @{
		/*!
			Predicate type. Must provide the function operator

			`Bounded_side operator()(Point_2 p),`

			which returns the position of point `p` relative to the half-open
			regular hyperbolic octagon.
		*/
		typedef unspecified_type 				Side_of_original_octagon;
	/// @}


	/// \name Construction Types
	/// @{
		/*!
			Construction type. Must provide the function operator

			`Point_2 operator() ( const Point_2& pt, const HyperbolicOctagonTranslation& tr ) const,`

			which returns the image of the point `pt` under the action of the 
			hyperbolic translation `tr`.
		*/
		typedef unspecified_type		      	Construct_point_2;
	/// @}
	


	/// \name Operations
	/// The following functions give access to the predicate objects:
	/// @{
		Side_of_original_octagon 
		side_of_original_octagon_object() const;

	/// @}

	/// \name
	/// The following functions give access to the construction objects:
	/// @{
		Construct_point_2 
		construct_point_2_object() const;
	/// @}

}; 

