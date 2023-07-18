/*!
\ingroup PkgSurfaceMeshParameterizationConcepts
\cgalConcept

`Parameterizer_3` is a concept of parameterization object for a given type
of mesh, `TriangleMesh`, which must be a model of the `FaceGraph` concept.

Border parameterizers are also models of this concept but they only parameterize
the border of a given mesh.

\cgalHeading{Creation}

Construction and destruction are undefined.

\cgalHasModelsBegin
\cgalHasModels{CGAL::Surface_mesh_parameterization::Fixed_border_parameterizer_3<TriangleMesh, BorderParameterizer, SolverTraits>}
\cgalHasModels{CGAL::Surface_mesh_parameterization::ARAP_parameterizer_3<TriangleMesh, BorderParameterizer, SolverTraits>}
\cgalHasModels{CGAL::Surface_mesh_parameterization::Barycentric_mapping_parameterizer_3<TriangleMesh, BorderParameterizer, SolverTraits>}
\cgalHasModels{CGAL::Surface_mesh_parameterization::Discrete_authalic_parameterizer_3<TriangleMesh, BorderParameterizer, SolverTraits>}
\cgalHasModels{CGAL::Surface_mesh_parameterization::Discrete_conformal_map_parameterizer_3<TriangleMesh, BorderParameterizer, SolverTraits>}
\cgalHasModels{CGAL::Surface_mesh_parameterization::LSCM_parameterizer_3<TriangleMesh, BorderParameterizer, SolverTraits>}
\cgalHasModels{CGAL::Surface_mesh_parameterization::Mean_value_coordinates_parameterizer_3<TriangleMesh, BorderParameterizer, SolverTraits>}
\cgalHasModels{CGAL::Surface_mesh_parameterization::Circular_border_parameterizer_3<TriangleMesh>}
\cgalHasModels{CGAL::Surface_mesh_parameterization::Square_border_parameterizer_3<TriangleMesh>}
\cgalHasModels{CGAL::Surface_mesh_parameterization::Two_vertices_parameterizer_3<TriangleMesh>}
\cgalHasModelsEnd

\sa `CGAL::Surface_mesh_parameterization::Orbifold_Tutte_parameterizer_3<SeamMesh, SolverTraits>`
*/

class Parameterizer_3
{
public:

  /// \name Types
  /// @{

  /// A given polygon mesh type, TriangleMesh, which is a model of the `FaceGraph` concept.
  typedef unspecified_type TriangleMesh;
  typedef boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;
  /// @}


  /// \name Operations
  /// @{

  /// Assign a 2D position (i.e.\ a `(u, v)` pair) on the shape to (some of)
  /// the vertices of the mesh. Mark them as <I>parameterized</I>.
  ///
  /// \tparam VertexUVmap must be a model of `ReadWritePropertyMap` with
  ///         `boost::graph_traits<TriangleMesh>::%vertex_descriptor` as key type and
  ///         %Point_2 (type deduced from `TriangleMesh` using `Kernel_traits`)
  ///         as value type.
  /// \tparam VertexIndexMap must be a model of `ReadablePropertyMap` with
  ///         `boost::graph_traits<TriangleMesh>::%vertex_descriptor` as key type and
  ///         a unique integer as value type.
  /// \tparam VertexParameterizedMap must be a model of `ReadWritePropertyMap` with
  ///         `boost::graph_traits<TriangleMesh>::%vertex_descriptor` as key type and
  ///         a Boolean as value type.
  ///
  /// \param mesh a triangulated surface.
  /// \param bhd a halfedge descriptor on the boundary of `mesh`.
  /// \param uvmap an instantiation of the class `VertexUVmap`.
  /// \param vimap an instantiation of the class `VertexIndexMap`.
  /// \param vpmap an instantiation of the class `VertexParameterizedMap`.
  ///
  /// \pre `mesh` must be a triangular mesh.
  /// \pre The vertices must be indexed (`vimap` must be initialized)
  ///
  template<typename VertexUVMap, typename VertexIndexMap, typename VertexParameterizedMap>
  Error_code parameterize(const TriangleMesh& mesh,
                          halfedge_descriptor bhd,
                          VertexUVMap uvmap,
                          VertexIndexMap vimap,
                          VertexParameterizedMap vpmap);

  /// @}

}; /* end Parameterizer_3 */
