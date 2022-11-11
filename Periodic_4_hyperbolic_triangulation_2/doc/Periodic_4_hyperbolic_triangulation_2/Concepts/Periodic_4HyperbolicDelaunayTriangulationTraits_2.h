// Copyright (c) 2016-2017 INRIA Nancy - Grand Est (France).
// All rights reserved.

/*!
\ingroup PkgPeriodic4HyperbolicTriangulation2Concepts
\cgalConcept

\cgalRefines `Periodic_4HyperbolicTriangulationTraits_2`

The concept `Periodic_4HyperbolicDelaunayTriangulationTraits_2` adds a requirement
to `Periodic_4HyperbolicTriangulationTraits_2` that needs to be fulfilled
by any class used to instantiate the first template parameter of the class
`CGAL::Periodic_4_hyperbolic_Delaunay_triangulation_2`.

\cgalHasModel `CGAL::Periodic_4_hyperbolic_Delaunay_triangulation_traits_2`
*/



class Periodic_4HyperbolicDelaunayTriangulationTraits_2 {

public:

        /// \name Computation Types
        /// @{
                /*!
                        Must provide the function operator

                        `double operator()(Hyperbolic_point_2 p1, Hyperbolic_point_2 p2, Hyperbolic_point_2 p3),`

                        which returns a floating-point approximation of the hyperbolic diameter
                        of the circle defined by the points `p1, p2,` and `p3`.
                */
                typedef unspecified_type                                 Compute_approximate_hyperbolic_diameter;
        /// @}


        /// \name Operations
                /// @{

                Compute_approximate_hyperbolic_diameter
                compute_approximate_hyperbolic_diameter_object() const;
        ///        @}

};

