namespace CGAL{
//! \ingroup PkgMesh_3IOFunctions
//! 
//! \brief write_polydata_3 writes the content of a triangulated surface mesh in the .vtp
//! XML format.
//! 
//! \tparam TriangleMesh a model of `FaceListGraph` with triangle faces.
//! \tparam NamedParameters a sequence of \ref pmp_namedparameters "Named Parameters"
//! 
//! \param os a `std::ostream`.
//! \param mesh an instance of `TriangleMesh` to be written.
//! \param binary decides if the data should be written in binary(`true`)
//!   or in ASCII(`false`).
//! \param np optional sequence of \ref pmp_namedparameters "Named Parameters" among the 
//! ones listed below
//!
//! \cgalNamedParamsBegin
//!    \cgalParamBegin{vertex_point_map} the property map with the points associated to
//! the vertices of `mesh`. If this parameter is omitted, an internal property map for
//!       `CGAL::vertex_point_t` must be available in `TriangleMesh`.
//!     \cgalParamEnd
//!    \cgalParamBegin{vertex_index_map} the property map with the indices associated to 
//! the vertices of `mesh`. If this parameter is omitted, an internal property map for
//!       `CGAL::vertex_index_t` must be available in `TriangleMesh`.
//!     \cgalParamEnd
//! \cgalNamedParamsEnd
template<class TriangleMesh, 
         class NamedParameters>
void write_polydata(std::ostream& os,
                    const TriangleMesh& mesh,
                    bool binary,
                    const NamedParameters& np);

//! \ingroup PkgMesh_3IOFunctions
//! 
//! \brief write_unstructured_grid_3 writes the content of a `C3t3` in the .vtu
//! XML format.
//! 
//! \tparam C3T3 a model of `MeshComplexWithFeatures_3InTriangulation_3`.
//! 
//! \param os a `std::ostream`.
//! \param c3t3 an instance of `C3T3` to be written.
//!
template <class C3T3>
void write_unstructured_grid_3(std::ostream& os,
                             const C3T3& c3t3);
}
