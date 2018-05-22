// Copyright (c) 2016-2017 INRIA Nancy - Grand Est (France).
// All rights reserved.

/*!
\ingroup PkgPeriodic4HyperbolicTriangulation2Concepts
\cgalConcept

A refinement of the concept TriangulationVertexBase_2 that adds an interface for hyperbolic translations.

At the base level, a vertex stores an input point and a handle to one of its incident faces.

For periodic hyperbolic triangulations, the vertex base class needs to temporarily store a hyperbolic 
translation during the insertion process.

\cgalHasModel CGAL::Periodic_4_hyperbolic_triangulation_ds_vertex_base_2

\sa `TriangulationDataStructure_2`
\sa `Periodic_4HyperbolicTriangulationDSFaceBase_2`

*/


class Periodic_4HyperbolicTriangulationDSVertexBase_2 : public TriangulationDSVertexBase_2 {
public:

  /// \name Types
  /// @{
	typedef typename Vb::Triangulation_data_structure   Triangulation_data_structure;
	typedef typename Triangulation_data_structure::Face_handle 		             
                                                      Face_handle;
	typedef typename GT::Point_2			                  Point;
	typedef typename GT::Hyperbolic_translation 			  Hyperbolic_translation;

  /*!
    Type used internally to break cyclical dependencies.
    \sa TriangulationDSFaceBase_2::Rebind_TDS
  */
	typedef unspecified_type                            Rebind_TDS;
  /// @}

  /// \name Creation
  /// @{
  /*!
    Default constructor.
  */
	Periodic_4HyperbolicTriangulationDSVertexBase_2();

  /*!
    Construct a vertex that stores the point `p`.
  */
	Periodic_4HyperbolicTriangulationDSVertexBase_2(const Point & p);

  /*!
    Constructs a vertex that stores the point `p` and is incident to the face `fh`.
  */
  Periodic_4HyperbolicTriangulationDSVertexBase_2(const Point & p, Face_handle fh);

  /*!
    Constructs a vertex that is incident to the face `fh`.
  */
	Periodic_4HyperbolicTriangulationDSVertexBase_2(const Face_handle& fh);
  /// @}


  /// \name Hyperbolic translations API
  /// @{
  /*!
    Stores the translation `tr` in the vertex, and sets a flag indicating that the vertex stores a translation.
  */
	void set_translation(Hyperbolic_translation tr);

  /*!
    Returns the translation stored in the vertex.
  */
	Hyperbolic_translation translation();

  /*!
    Removes the translation stored in the vertex, and resets the translation flag.
  */
	void clear_translation();

  /*!
    Returns the value of a flag, indicating whether the vertex stores a translation or not.
  */
	bool get_translation_flag();
  ///@}

};
