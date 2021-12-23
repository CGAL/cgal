// Copyright (c) 2016-2017 INRIA Nancy - Grand Est (France).
// All rights reserved.

/*!
\ingroup PkgHyperbolicTriangulation2Concepts
\cgalConcept

\cgalRefines TriangulationFaceBase_2

The concept `HyperbolicTriangulationFaceBase_2` describes the requirements for the base
face class of a hyperbolic triangulation data structure.

This concept provides an interface for the functionality needed in faces to compute
Delaunay triangulations in the hyperbolic plane. The function `tds_data()` is used
internally by the triangulation class during the insertion of points in the triangulation.

\cgalHasModel `CGAL::Hyperbolic_triangulation_face_base_2`

\sa `TriangulationDataStructure_2`

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
      This function gives non-`const` access to a variable of type `CGAL::Object`.
      \cgalAdvancedFunction
    */
    CGAL::Object& tds_data();

    /*!
      This function gives `const` access to a variable of type `CGAL::Object`.
      \cgalAdvancedFunction
    */
    const CGAL::Object& tds_data() const;
  /// @}

};
