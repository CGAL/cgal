namespace CGAL {
//! \ingroup PkgMesh_3Functions
//!
//! Gets reconstructed surface out of a `MeshComplexWithFeatures_3InTriangulation_3` object.
//!
//! This variant exports the surface as a `TriangleMesh` and appends it to `graph`, using
//! `orient_polygon_soup()`.
//!
//! @tparam C3T3 model of the `MeshComplexWithFeatures_3InTriangulation_3` concept.
//! @tparam TriangleMesh a model of `MutableFaceGraph` with an internal point property map
//!
//! @param c3t3 an instance of a `C3T3`.
//! @param graph an instance of `TriangleMesh`.
template<class C3T3, class FaceGraph>
void output_c3t3_to_facegraph(const C3T3& c3t3, FaceGraph& graph);
}
