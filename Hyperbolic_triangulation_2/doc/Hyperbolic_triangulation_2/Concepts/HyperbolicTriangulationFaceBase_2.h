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
    \cgalModifBegin
      Sets a boolean flag stored in the face.
      This is a low-level function used internally by the hyperbolic triangulation.
    \cgalModifEnd
    */
    void set_flag(bool flag);

    /*!
    \cgalModifBegin
      Retrieves the value of a boolean flag stored in the face.
      This is a low-level function used internally by the hyperbolic triangulation.
    \cgalModifEnd
    */
    bool get_flag() const;
    

    /*!
    \cgalModifBegin
      Sets the value of an unsigned char stored in the face.
      This is a low-level function used internally by the hyperbolic triangulation.
    \cgalModifEnd
    */
    void set_char(unsigned char uschar);

    /*! 
    \cgalModifBegin
      Retrieves the value of an unsigned char stored in the face.
      This is a low-level function used internally by the hyperbolic triangulation.
      \cgalModifEnd
    */
    unsigned char get_char() const;
    
  /// @}

};
