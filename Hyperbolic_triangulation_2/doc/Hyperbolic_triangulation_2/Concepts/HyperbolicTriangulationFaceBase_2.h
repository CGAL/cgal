// Copyright (c) 2016-2017 INRIA Nancy - Grand Est (France).
// All rights reserved.

/*!
\ingroup PkgHyperbolicTriangulation2Concepts
\cgalConcept

The concept HyperbolicTriangulationFaceBase_2 describes the requirements for the base 
face class of a hyperbolic triangulation data structure. 

This concept provides access to vertices incident to the face, as well as to its neighboring 
faces. Compare with Section \ref Section_2D_Triangulations_Software_Design in the 2D 
Triangulation User Manual.

\cgalHasModel CGAL::Hyperbolic_triangulation_face_base_2
\cgalHasModel CGAL::Hyperbolic_triangulation_face_base_with_info_2

\sa `TriangulationDataStructure_2`

*/


class HyperbolicTriangulationFaceBase_2 {

public:

  /// \name Types
  /// @{
    typedef unspecified_type    Vertex_handle;
    typedef unspecified_type    Face_handle;
  // @}

  /// \name Creation
  /// @{

    /*!
    Default constructor.
    */
    Hyperbolic_triangulation_face_base_2();

    /*!
      Creates a face to which the vertices `v0, v1, v2` are incident.
    */
    Hyperbolic_triangulation_face_base_2(Vertex_handle v0, 
  			                                 Vertex_handle v1, 
  			                                 Vertex_handle v2);

    /*!
      Creates a face to which the vertices `v0, v1, v2` are incident, and the faces `n0, n1, n2` are neighbors.
    */
    Hyperbolic_triangulation_face_base_2(Vertex_handle v0, 
                    			               Vertex_handle v1, 
                    			               Vertex_handle v2,
                    			               Face_handle n0, 
                    			               Face_handle n1, 
                    			               Face_handle n2);
  /// @}


  /// \name Access functions
  /// @{

    /*!
      Returns `true` if the face is finite, but not hyperbolic, i.e., if its circumscribing circle intersects the circle at infinity
     */
    bool is_finite_non_hyperbolic() const;
    
    /*!
      Sets a flag indicating whether the face is finite but not hyperbolic
    */
    void set_finite_non_hyperbolic(bool is_finite_non_hyperbolic);
    
    /*!
      Returns `true` if the face has a non-hyperbolic edge, i.e., an edge whose circumscribing disk intersects the circle at infinity
    */
    bool has_non_hyperbolic_edge() const;
    
    /*! 
      Returns the index of the vertex opposite to the non-nyperbolic edge of the face
      \pre The face is non-hyperbolic
    */
    unsigned char get_non_hyperbolic_edge() const;
    
    /*!
      Sets the index of the non-hyperbolic edge of the face
      \pre The face is non-hyperbolic
    */
    void set_non_hyperbolic_edge(unsigned char non_hyperbolic_edge);
  /// @}

};
