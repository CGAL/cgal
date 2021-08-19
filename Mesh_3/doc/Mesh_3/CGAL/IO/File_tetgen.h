namespace CGAL{
namespace IO {
/**
 * \ingroup PkgMesh3IOFunctions
 * @brief outputs a mesh complex to tetgen format
 * @param filename the path to the output file
 * @param c3t3 the mesh
 * @param rebind if true, labels of cells are rebinded into [1..nb_of_labels]
 * @param show_patches if true, patches are labeled with different labels than
 * cells. If false, each surface facet is written twice, using label of
 * each adjacent cell.
 * \see \ref IOStreamTetgen
 */
template <class C3T3>
void
output_to_tetgen(std::string filename,
                 const C3T3& c3t3,
                 bool rebind = false,
                 bool show_patches = false);
} }
