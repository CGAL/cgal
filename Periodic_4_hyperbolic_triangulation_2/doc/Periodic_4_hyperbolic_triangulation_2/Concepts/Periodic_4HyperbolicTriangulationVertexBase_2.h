// Copyright (c) 2016-2017 INRIA Nancy - Grand Est (France).
// All rights reserved.

/*!
\ingroup PkgPeriodic4HyperbolicTriangulation2Concepts
\cgalConcept

\cgalRefines{TriangulationVertexBase_2}

A refinement of the concept `TriangulationVertexBase_2` that adds an interface for hyperbolic translations.

For periodic hyperbolic triangulations, the vertex base class needs to temporarily store a hyperbolic
translation during the insertion process.
A boolean flag indicates whether the face stores a translation or not. The value of the flag is automatically
set when storing or removing a translation.

\cgalHasModelsBegin
\cgalHasModels{CGAL::Periodic_4_hyperbolic_triangulation_vertex_base_2}
\cgalHasModelsEnd

\sa `TriangulationDataStructure_2`
\sa `Periodic_4HyperbolicTriangulationFaceBase_2`

*/


class Periodic_4HyperbolicTriangulationVertexBase_2 {
public:

  /// \name Types
  /// @{

  typedef TriangulationVertexBase_2::Face_handle                            Face_handle;
  typedef Periodic_4HyperbolicDelaunayTriangulationTraits_2    Geom_traits;
        typedef Geom_traits::Hyperbolic_point_2                             Point;
  typedef Geom_traits::Hyperbolic_translation                  Hyperbolic_translation;
  /// @}

  /// \name Creation
  /// @{
  /*!
    Default constructor.
  */
        Periodic_4HyperbolicTriangulationVertexBase_2();

  /*!
    Construct a vertex that stores the point `p`.
  */
        Periodic_4HyperbolicTriangulationVertexBase_2(const Point & p);

  /*!
    Constructs a vertex that stores the point `p` and is incident to the face `fh`.
  */
  Periodic_4HyperbolicTriangulationVertexBase_2(const Point & p, Face_handle fh);

  /*!
    Constructs a vertex that is incident to the face `fh`.
  */
        Periodic_4HyperbolicTriangulationVertexBase_2(const Face_handle& fh);
  /// @}


  /// \name Hyperbolic translations API
  /// @{
  /*!
    Stores the translation `tr` in the vertex, and sets the translation flag to `true`.
  */
        void set_translation(const Hyperbolic_translation& tr);

  /*!
    Returns the translation stored in the vertex.
  */
        Hyperbolic_translation translation();

  /*!
    Removes the translation stored in the vertex, and sets the translation flag to `false`.
  */
        void clear_translation();

  /*!
    Returns the value of the translation flag.
  */
        bool get_translation_flag();
  ///@}

};
