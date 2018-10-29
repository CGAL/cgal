namespace CGAL{
//!\ingroup PkgMesh2
//! \brief writes the content of a `CDT` in the .vtu
//! XML format. 
//! 
//! The triangles inside the domain and the constrained edges will be outputted.
//! 
//! \tparam CDT a `Constrained_Delaunay_triangulation_2`.
//! 
//! \param os a `std::ostream`.
//! \param tr an instance of `CDT` to be written.
//! \param binary decides if the data should be written in binary(`true`)
//!   or in ASCII(`false`).
//!
template <class CDT>
void write_VTU(std::ostream& os,
                               const CDT& tr,
                               bool binary = true);
}
