namespace CGAL{
//! \ingroup PkgMesh_3IOFunctions
//! 
//! \brief writes the content of a `C3t3` in the .vtu
//! XML format.
//! 
//! \tparam C3T3 a model of `MeshComplexWithFeatures_3InTriangulation_3`.
//! 
//! \param os a `std::ostream`.
//! \param c3t3 an instance of `C3T3` to be written.
//! \param binary decides if the data should be written in binary(`true`)
//!   or in ASCII(`false`).
//!
template <class C3T3>
void output_to_vtu(std::ostream& os,
                   const C3T3& c3t3,
                   bool binary = true);
}
