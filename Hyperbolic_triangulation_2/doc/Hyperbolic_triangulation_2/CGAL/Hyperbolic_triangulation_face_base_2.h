// Copyright (c) 2016-2017 INRIA Nancy - Grand Est (France).
// All rights reserved.

namespace CGAL { 

/*!
\ingroup PkgHyperbolicTriangulation2VertexFaceClasses

The class `Hyperbolic_triangulation_face_base_2` is designed as one of the default models for the
concept `HyperbolicTriangulationFaceBase_2`.

\tparam Gt must be a model of `HyperbolicDelaunayTriangulationTraits_2`.
\tparam Fb must be a model of  `TriangulationDSFaceBase_2`. Defaults to `Triangulation_ds_face_base_2<>`.


\sa Hyperbolic_triangulation_face_base_with_info_2

\cgalModels HyperbolicTriangulationFaceBase_2
*/


template < typename Gt, typename Fb >
class Hyperbolic_triangulation_face_base_2 : public Fb {

public:

  /// \name Creation
  /// @{
    /*!
      Default constructor
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
      Sets the flag indicating whether the face is finite but not hyperbolic
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


} //namespace CGAL 

