
namespace CGAL {


/*!
\ingroup PkgMesh3MeshClasses
<!-- Meta-comment: this class cannot be deprecated by
Compact_mesh_cell_base_3, because the later has a different API.
-- Laurent Rineau, 2013/10/16
\deprecated This class is deprecated since \cgal 4.3. Use
`CGAL::Compact_mesh_cell_base_3<Gt,MD,Tds>` instead.
-->

The class `Mesh_cell_base_3<Gt, MD, Cb>` is a model of the concept `MeshCellBase_3`.
It is designed to serve as cell base class for the 3D triangulation
used in the 3D mesh generation process.

\tparam Gt is the geometric traits class.
It has to be a model of the concept `MeshTriangulationTraits_3`.

\tparam MD provides the types of indices used to identify
the faces of the input complex. It has to be a model
of the concept `MeshDomain_3`.

\tparam Cb is the cell base class. It has to be a model
of the concept `RegularTriangulationCellBaseWithWeightedCircumcenter_3` and defaults to
`Regular_triangulation_cell_base_with_weighted_circumcenter_3<Gt>`.

\cgalModels `MeshCellBase_3`

\sa `CGAL::Mesh_complex_3_in_triangulation_3<Tr,CornerIndex,CurveIndex>`
\sa `CGAL::Compact_mesh_cell_base_3<Gt, MD, Tds>`

*/
template< typename Gt,  typename MD, typename Cb >
class Mesh_cell_base_3 : public Cb {
public:

}; /* end Mesh_cell_base_3 */
} /* end namespace CGAL */
