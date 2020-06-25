namespace CGAL {

/// \ingroup PkgPeriodic3Mesh3IOFunctions
///
/// \brief outputs a periodic mesh to the medit (`.mesh`) file format.
///        See \cgalCite{frey:inria-00069921} for a comprehensive description of this file format.
///
/// \param os the output stream
/// \param c3t3 the mesh complex
/// \param occurrence_count the number of copies that are drawn
/// \param distinguish_copies if set to `true`, each copy is assigned a unique color.
///                           Otherwise, all domains are drawn with subdomain index-based colors.
/// \param rebind if set to `true`, labels of cells are rebinded into `[1..nb_of_labels]`
/// \param show_patches if set to `true`, patches are labeled with different labels than
///                     cells. If `false`, each surface facet is written twice,
///                     using the label of each adjacent cell.
///
template <class C3T3>
void output_periodic_mesh_to_medit(std::ostream& os,
                                   const C3T3& c3t3,
                                   int occurrence_count = 8,
                                   bool distinguish_copies = true,
                                   bool rebind = false,
                                   bool show_patches = false);

} // namespace CGAL
