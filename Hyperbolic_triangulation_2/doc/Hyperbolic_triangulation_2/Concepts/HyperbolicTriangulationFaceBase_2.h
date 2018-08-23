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


  /// \name Internal access functions
  /// \cgalAdvancedBegin
  /// These functions are used internally by the hyperbolic Delaunay triangulation.
  /// The user is not encouraged to use them directly as they may change in the future.
  /// \cgalAdvancedEnd
  /// @{
    
    /*!
    \cgalModifBegin
      This function gives non-`const` access to a variable of type `CGAL::Object`.
      \cgalAdvancedFunction
    \cgalModifEnd
    */
    CGAL::Object& tds_data();

    /*!
    \cgalModifBegin
      This function gives `const` access to a variable of type `CGAL::Object`.
      \cgalAdvancedFunction
    \cgalModifEnd
    */
    const CGAL::Object& tds_data() const;
  /// @}
  
};
