
namespace CGAL {

/*!
\ingroup PkgTriangulation3VertexCellClasses

The class `Triangulation_vertex_base_3` is a model of the concept 
`TriangulationVertexBase_3`, the base vertex of a 3D-triangulation. 
This class stores a point. 

This class can be used directly or can serve as a base to derive other classes 
with some additional attributes (a color for example) tuned for a specific 
application. 


\tparam TriangulationTraits_3 is the geometric traits class 
 which provides the point type `Point_3`. 
Users of the geometric triangulations are strongly advised to use the same 
geometric traits class `TriangulationTraits_3` as the one used for 
`Triangulation_3`. This way, the point type defined by the base vertex is 
the same as the point type defined by the geometric traits class. 

\tparam TriangulationDSVertexBase_3_ is a combinatorial vertex base class from which 
`Triangulation_vertex_base_3` derives.
It must be a model of the `TriangulationDSVertexBase_3` concept.
It has the default value `Triangulation_ds_vertex_base_3<TDS>`.

\cgalModels `TriangulationVertexBase_3`

\sa `CGAL::Triangulation_cell_base_3` 
\sa `CGAL::Triangulation_ds_vertex_base_3` 
\sa `CGAL::Triangulation_vertex_base_with_info_3` 

*/
template< typename TriangulationTraits_3, typename TriangulationDSVertexBase_3_ >
class Triangulation_vertex_base_3 : public TriangulationDSVertexBase_3_ {
public:

/// \name Types 
/// @{

/*!

*/ 
typedef TriangulationTraits_3::Point_3 Point; 

/// @}

}; /* end Triangulation_vertex_base_3 */
} /* end namespace CGAL */
