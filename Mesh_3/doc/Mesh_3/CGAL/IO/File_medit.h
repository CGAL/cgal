namespace CGAL {
namespace IO {
/// \ingroup PkgMesh3IOFunctions
///
/// \brief outputs a mesh complex to the medit (`.mesh`) file format.
///        See \cgalCite{frey:inria-00069921} for a comprehensive description of this file format.
///
/// \param os the output stream
/// \param c3t3 the mesh complex
/// \param rebind if set to `true`, labels of cells are rebinded into `[1..nb_of_labels]`
/// \param show_patches if set to `true`, patches are labeled with different labels than
///                     cells. If set to `false`, each surface facet is written twice,
///                     using the label of each adjacent cell.
///
template <class C3T3>
void output_to_medit(std::ostream& os,
                     const C3T3& c3t3,
                     bool rebind = false,
                     bool show_patches = false);

}} // end namespace CGAL::IO
