namespace CGAL{
//!\ingroup PkgMesh2
//! \brief writes the faces of a domain and its constrained edges embedded in
//! a 2D constrained Delaunay triangulation using the `PolyData` XML format.
//!
//! The faces output are those for which `DelaunayMeshFaceBase_2::is_in_domain()` returns `true`,
//! the edges are those for which `ConstrainedTriangulationFaceBase_2::is_constained()` returns `true`.
//! \tparam CDT a `Constrained_Delaunay_triangulation_2` with face type model of `DelaunayMeshFaceBase_2`.
//!
//! \param os the stream used for writting.
//! \param tr the triangulated domain to be written.
//! \param binary decides if the data should be written in binary (`true`)
//!   or in ASCII (`false`).
//!
template <class CDT>
void write_VTU(std::ostream& os,
               const CDT& tr,
               bool binary = true);
}
