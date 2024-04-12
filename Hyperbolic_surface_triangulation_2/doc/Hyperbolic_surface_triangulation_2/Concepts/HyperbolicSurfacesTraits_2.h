// Copyright (c) 2024 INRIA Nancy - Grand Est (France). LIGM Marne-la-Vallée (France)
// All rights reserved.

/*!
\ingroup PkgHyperbolicSurfaceTriangulation2Concepts
\cgalConcept

\cgalRefines{HyperbolicDelaunayTriangulationTraits_2}

The concept `HyperbolicSurfacesTraits_2` describes the set of requirements
to be fulfilled by any class instantiating the first template parameter of the main classes of the package.
Notably, it defines the field type, the complex number type, and the point type used by thoses classes.

\cgalHasModelsBegin
\cgalHasModels{CGAL::Hyperbolic_surfaces_traits_2}
\cgalHasModelsEnd
*/

class HyperbolicSurfacesTraits_2 {
public:
  /// \name Types
  /// @{
  /*!
          Field number type.
  */
  typedef unspecified_type                          FT;

  /*!
          Represents a point in the Poincaré disk model.
  */
  typedef unspecified_type          Hyperbolic_point_2;

  /*!
          Represents a complex number over FT.
  */
  typedef unspecified_type               Complex;
  /// @}
};
