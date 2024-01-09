// Copyright (c) 2016-2017 INRIA Nancy - Grand Est (France).
// All rights reserved.

/*!
\ingroup PkgHyperbolicTriangulation2Concepts
\cgalConcept

The concept `HyperbolicTriangulationFaceBase_2` describes the requirements for an internal
data class that must be provided by the model of the `HyperbolicTriangulationFaceBase_2`.

\sa `HyperbolicTriangulationFaceBase_2`
*/
class HyperbolicFaceData {
public:
  /// sets the face as infinite and non-Delaunay hyperbolic
  void set_infinite();

  /// sets the face as Delaunay hyperbolic
  void set_Delaunay_hyperbolic();

  /// returns whether the face is Delaunay hyperbolic, or not.
  bool is_Delaunay_hyperbolic();

  /// sets the `i`-th edge of the face as non-hyperbolic
  void set_Delaunay_non_hyperbolic(int i);

  /// returns whether the `i`-th edge of the face is non-hyperbolic
  bool is_Delaunay_non_hyperbolic(int i);
};

/*!
\ingroup PkgHyperbolicTriangulation2Concepts
\cgalConcept

\cgalRefines{TriangulationFaceBase_2}

The concept `HyperbolicTriangulationFaceBase_2` describes the requirements for the base
face class of a hyperbolic triangulation data structure.

This concept provides an interface for the functionality needed in faces to compute
Delaunay triangulations in the hyperbolic plane. The function `tds_data()` is used
internally by the triangulation class during the insertion of points in the triangulation.

\cgalHasModelsBegin
\cgalHasModels{CGAL::Hyperbolic_triangulation_face_base_2}
\cgalHasModelsEnd

\sa `TriangulationDataStructure_2`
\sa `HyperbolicFaceData`
*/
class HyperbolicTriangulationFaceBase_2 {
public:

  /// \name Internal Access Functions
  /// \cgalAdvancedBegin
  /// These functions are used internally by the hyperbolic Delaunay triangulation.
  /// The user is not encouraged to use them directly as they may change in the future.
  /// \cgalAdvancedEnd
  /// @{

    /*!
      This function gives non-`const` access to a variable that is a model of `HyperbolicFaceData`.
      \cgalAdvancedFunction
    */
    HyperbolicFaceData& hyperbolic_data();

    /*!
      This function gives `const` access to a variable of that is a model of `HyperbolicFaceData`.
      \cgalAdvancedFunction
    */
    const HyperbolicFaceData& hyperbolic_data() const;
  /// @}

};
