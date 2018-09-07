// Copyright (c) 2016-2017 INRIA Nancy - Grand Est (France).
// All rights reserved.

/*!
\ingroup PkgPeriodic4HyperbolicTriangulation2Concepts
\cgalConcept

\cgalRefines Periodic_4HyperbolicTriangulationTraits_2

The concept `Periodic_4HyperbolicDelaunayTriangulationTraits_2` refines the concept 
`Periodic_4HyperbolicDelaunayTriangulationTraits_2`. It adds a requirement that needs
to be fulfilled by any class used to instantiate the first template parameter of the class 
`CGAL::Periodic_4_hyperbolic_Delaunay_triangulation_2`. The added requirement is 
essential for the removal of existing vertices in a Delaunay triangulation of the
Bolza surface.

\cgalModifBegin
	The concept requires that two nested types are provided:
	<ul>
		<li> A field number type `FT` that must support exact computations with algebraic 
			 numbers, in particular with nested square roots;
		<li> A `Hyperbolic_translation` type, which satisfies the requirements described 
			 in the concept `HyperbolicOctagonTranslation`.
	</ul>
\cgalModifEnd

\cgalHasModel CGAL::Periodic_4_hyperbolic_Delaunay_triangulation_traits_2
*/



class Periodic_4HyperbolicDelaunayTriangulationTraits_2 {

public:

	/// \name Computation Types
	/// @{
		/*!
			Type that permits to compute the hyperbolic diameter of a circle.
			Must provide the function operator

			`double operator()(Circle_2 c),`

			which returns a floating-point approximation of the hyperbolic diameter
			of the circle `c`.
		*/
		typedef unspecified_type 				Compute_approximate_hyperbolic_diameter;
	/// @}


	/// \name Operations
	/// \cgalModifBegin
	/// The following function gives access to the computation object:
	/// \cgalModifEnd
	/// @{

		Compute_approximate_hyperbolic_diameter	
		compute_approximate_hyperbolic_diameter_object() const;
	///	@}

}; 

