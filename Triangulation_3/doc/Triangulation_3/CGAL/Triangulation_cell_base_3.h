
namespace CGAL {

/*!
\ingroup PkgTriangulation3VertexCellClasses

The class `Triangulation_cell_base_3` is a model of the concept 
`TriangulationCellBase_3`, the base cell of a 3D-triangulation. 

This class can be used directly or can serve as a base to derive other classes 
with some additional attributes (a color for example) tuned for a specific 
application. 


\tparam TriangulationTraits_3 is the geometric traits class. It is actually not used by this class. 

\tparam TriangulationDSCellBase_3 is a combinatorial cell base class from which 
`Triangulation_cell_base_3` derives. 
It has the default value `Triangulation_ds_cell_base_3<>`. 

Note that this model does not store the circumcenter, but computes it 
every time the circumcenter function is called. See 
`Triangulation_cell_base_with_circumcenter_3` for a way to cache the 
circumcenter computation. 

\cgalModels ::TriangulationCellBase_3 

\sa `CGAL::Triangulation_ds_cell_base_3` 
\sa `CGAL::Triangulation_cell_base_with_info_3` 
\sa `CGAL::Triangulation_cell_base_with_circumcenter_3` 
\sa `CGAL::Triangulation_vertex_base_3` 

*/
template< typename TriangulationTraits_3, typename TriangulationDSCellBase_3 >
class Triangulation_cell_base_3 : public TriangulationDSCellBase_3 {
public:

/// @}

}; /* end Triangulation_cell_base_3 */
} /* end namespace CGAL */
