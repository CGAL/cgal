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
//!
template <class C3T3>
void write_VTU(std::ostream& os,
                             const C3T3& c3t3);
}
