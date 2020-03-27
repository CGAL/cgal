namespace CGAL{
/*!
 *Exports `tr` to `os` in the OFF format.
 * \see \ref IOStreamOFF
 */
template < class GT, class TDS >
std::ostream &
export_triangulation_to_off(std::ostream & os,
                            const Triangulation<GT,TDS> & tr,
                            bool in_3D_export_surface_only = false);
}
