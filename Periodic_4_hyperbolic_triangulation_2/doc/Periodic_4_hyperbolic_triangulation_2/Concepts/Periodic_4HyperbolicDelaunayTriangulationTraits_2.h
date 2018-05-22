// Copyright (c) 2016-2017 INRIA Nancy - Grand Est (France).
// All rights reserved.

/*!
\ingroup PkgPeriodic4HyperbolicTriangulation2Concepts
\cgalConcept

The concept `Periodic_4HyperbolicDelaunayTriangulationTraits_2` refines the concept 
`HyperbolicDelaunayTriangulationTraits_2`. It describes the set of requirements to be 
fulfilled by any class used to instantiate the first template parameter of the class 
`CGAL::Periodic_4_hyperbolic_Delaunay_triangulation_2`. In addition to the geometric 
types and the operations defined on them in `HyperbolicDelaunayTriangulationTraits_2`, 
it defines the hyperbolic translations that allow to encode the periodicity of the 
triangulation.

\cgalHasModel CGAL::Periodic_4_hyperbolic_Delaunay_triangulation_traits_2
*/



class Periodic_4HyperbolicDelaunayTriangulationTraits_2 : public HyperbolicDelaunayTriangulationTraits_2 {

public:


	/// \name Types
	/// @{
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
	
		
	/// \name
	/// The following type is not really a construction, nor a predicate.
	/// @{
		/*!
			Type that permits to compute the hyperbolic diameter of a circle.
			Must provide the function operator

			`double operator()(Circle_2 c),`

			which returns a floating-point approximation of the hyperbolic diameter
			of the circle `c`.
			
			\todo Move to the triangulation class instead?
		*/
		typedef unspecified_type 				Compute_approximate_hyperbolic_diameter;
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

