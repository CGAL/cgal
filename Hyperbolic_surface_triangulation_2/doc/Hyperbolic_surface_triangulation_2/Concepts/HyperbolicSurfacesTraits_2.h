// Copyright (c) 2024 INRIA Nancy - Grand Est (France). LIGM Marne-la-Vallée (France)
// All rights reserved.

/*!
\ingroup PkgHyperbolicSurfaceTriangulation2Concepts
\cgalConcept

The concept `HyperbolicSurfacesTraits_2` describes the set of requirements
to be fulfilled by any class instantiating the template parameter of the main classes of the package.

\cgalRefines{HyperbolicDelaunayTriangulationTraits_2}

\cgalHasModelsBegin
\cgalHasModels{CGAL::Hyperbolic_surfaces_traits_2}
\cgalHasModelsEnd
*/
class HyperbolicSurfacesTraits_2 {
public:
  /// \name Types
  /// @{
  /*!
          Represents a complex number over a field: must be a model of 'CGAL::ComplexWithoutSqrt'.
  */
  typedef unspecified_type               Complex;
  /// @}
};
