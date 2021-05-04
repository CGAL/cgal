namespace CGAL{
namespace IO {
//!\ingroup PkgMesh2IO
//! \brief writes the faces of a domain and its constrained edges embedded in
//! a 2D constrained Delaunay triangulation using the `PolyData` XML format.
//!
//! The faces output are those for which `DelaunayMeshFaceBase_2::is_in_domain()` returns `true`,
//! the edges are those for which `ConstrainedTriangulationFaceBase_2::is_constrained()` returns `true`.
//! \tparam CDT a `Constrained_Delaunay_triangulation_2` with face type model of `DelaunayMeshFaceBase_2`.
//!
//! \param os the stream used for writing.
//! \param tr the triangulated domain to be written.
//! \param mode decides if the data should be written in binary (`BINARY`)
//!   or in ASCII (`ASCII`).
//!
template <class CDT>
void write_VTU(std::ostream& os,
               const CDT& tr,
               IO::Mode mode = IO::BINARY);
} }
