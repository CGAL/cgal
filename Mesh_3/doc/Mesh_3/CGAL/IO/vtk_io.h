namespace CGAL{
//! \ingroup PkgMesh_3IOFunctions
//!
//! \brief writes a tetrahedron mesh using the `UnstructuredGrid` XML format.
//!
//! \tparam C3T3 a model of `MeshComplexWithFeatures_3InTriangulation_3`.
//!
//! \param os the stream used for writting.
//! \param c3t3 the instance of `C3T3` to be written.
//! \param binary decides if the data should be written in binary (`true`)
//!   or in ASCII (`false`).
//!
template <class C3T3>
void output_to_vtu(std::ostream& os,
                   const C3T3& c3t3,
                   bool binary = true);
}
