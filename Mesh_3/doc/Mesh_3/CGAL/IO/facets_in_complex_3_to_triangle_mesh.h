namespace CGAL {
//! \ingroup PkgMesh3Functions
//!
//! \brief builds a `TriangleMesh` from the surface facets, with a consistent orientation at the interface of two subdomains.
//!
//! This function exports the surface as a `TriangleMesh` and appends it to `graph`, using
//! `orient_polygon_soup()`.
//!
//! @tparam C3T3 a model of `MeshComplexWithFeatures_3InTriangulation_3`.
//! @tparam TriangleMesh a model of `MutableFaceGraph` with an internal point property map. The point type should be compatible with the one used in `C3T3`.
//!
//! @param c3t3 an instance of `C3T3`.
//! @param graph an instance of `TriangleMesh`.
template<class C3T3, class TriangleMesh>
void facets_in_complex_3_to_triangle_mesh(const C3T3& c3t3, TriangleMesh& graph);
}
