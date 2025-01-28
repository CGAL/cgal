// Copyright (c) 2016-2017 INRIA Nancy - Grand Est (France).
// All rights reserved.

/*!
\ingroup PkgPeriodic4HyperbolicTriangulation2Concepts
\cgalConcept

\cgalRefines{HyperbolicDelaunayTriangulationTraits_2}

The concept `Periodic_4HyperbolicTriangulationTraits_2` describes the set of requirements
to be fulfilled by any class used to instantiate the first template parameter of the class
`CGAL::Periodic_4_hyperbolic_triangulation_2`. In addition to the geometric types and the
operations defined on them in `HyperbolicDelaunayTriangulationTraits_2`, it defines the
hyperbolic translations that enable the encoding of the periodicity of the triangulation.

The concept requires that the field number type `FT` defined in the concept
`HyperbolicDelaunayTriangulationTraits_2` supports exact computations with algebraic numbers,
in particular with nested square roots.

*/



class Periodic_4HyperbolicTriangulationTraits_2 {

public:


        /// \name Types
        /// @{

                /*!
                        Represents a hyperbolic translation.
                */
                typedef CGAL::Hyperbolic_octagon_translation<FT>                 Hyperbolic_translation;
        /// @}


        /// \name Predicate Types
        /// @{
                /*!
                        Predicate type. Must provide the function operator

                        `Bounded_side operator()(Hyperbolic_point_2 p),`

                        which returns the position of point `p` relative to the half-open
                        regular hyperbolic octagon.
                */
                typedef unspecified_type                                 Side_of_original_octagon;
        /// @}


        /// \name Construction Types
        /// @{
                /*!
                        Construction type. Must provide the function operator

                        `Hyperbolic_point_2 operator() ( const Hyperbolic_point_2& pt, const Hyperbolic_translation& tr ) const,`

                        which returns the image of the point `pt` under the action of the
                        hyperbolic translation `tr`.
                */
                typedef unspecified_type                              Construct_hyperbolic_point_2;
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
                Construct_hyperbolic_point_2
                construct_hyperbolic_point_2_object() const;
        /// @}

};

