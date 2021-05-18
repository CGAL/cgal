namespace CGAL{
namespace IO {
//! \ingroup PkgMesh3IOFunctions
//!
//! \brief writes a tetrahedron mesh using the `UnstructuredGrid` XML format.
//!
//! \tparam C3T3 a model of `MeshComplexWithFeatures_3InTriangulation_3`.
//!
//! \param os the stream used for writing.
//! \param c3t3 the instance of `C3T3` to be written.
//! \param mode decides if the data should be written in binary (`IO::BINARY`)
//!   or in ASCII (`IO::ASCII`).
//!
template <class C3T3>
void output_to_vtu(std::ostream& os,
                   const C3T3& c3t3,
                   IO::Mode mode = BINARY);
} }
