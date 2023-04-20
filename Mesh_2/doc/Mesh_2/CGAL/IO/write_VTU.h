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
//!             or in \ascii (`ASCII`).
//!
template <class CDT>
void write_VTU(std::ostream& os,
               const CDT& tr,
               IO::Mode mode = IO::BINARY);


//!\ingroup PkgMesh2IO
//! \brief writes the faces of a domain and its constrained edges embedded in
//! a 2D constrained Delaunay triangulation using the `PolyData` XML
//! format.
  //! The faces output are those for which `get(ipm, f)` returns
  //! `true` where `f` is a  `CDT::Face_handle`,
//! the edges are those for which `ConstrainedTriangulationFaceBase_2::is_constrained()` returns `true`.
//! \tparam CDT a `Constrained_Delaunay_triangulation_2` with face
//! type model of `DelaunayMeshFaceBase_2`.
//! \tparam InDomainPmap a class model of `ReadWritePropertyMap` with
//! `CDT::Face_handle` as key type and `bool` as value type.
//!
//! \param os the stream used for writing.
//! \param tr the triangulated domain to be written.
//! \param ipm the property map storing if a face is in the domain.
//! \param mode decides if the data should be written in binary (`BINARY`)
//!             or in \ascii (`ASCII`).
//!
template <class CDT, typename InDomainPmap>
void write_VTU(std::ostream& os,
               const CDT& tr,
               InDomainPmap ipm,
               IO::Mode mode = IO::BINARY);


} }
