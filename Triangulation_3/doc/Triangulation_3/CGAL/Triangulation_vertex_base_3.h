
namespace CGAL {

/*!
\ingroup PkgTriangulation3VertexCellClasses

The class `Triangulation_vertex_base_3` is a model of the concept 
`TriangulationVertexBase_3`, the base vertex of a 3D-triangulation. 
This class stores a point. 

This class can be used directly or can serve as a base to derive other classes 
with some additional attributes (a color for example) tuned for a specific 
application. 

### Parameters ###

The first template argument is the geometric traits class 
`TriangulationTraits_3` which provides the point type, `Point_3`. 
Users of the geometric triangulations (Section \ref TDS3secdesign and 
Chapter \ref chapterTriangulation3) are strongly advised to use the same 
geometric traits class `TriangulationTraits_3` as the one used for 
`Triangulation_3`. This way, the point type defined by the base vertex is 
the same as the point type defined by the geometric traits class. 

The second template argument is a combinatorial vertex base class from which 
`Triangulation_vertex_base_3` derives. 
It has the default value `Triangulation_ds_vertex_base_3<>`. 

\models ::TriangulationVertexBase_3 

\sa `CGAL::Triangulation_cell_base_3` 
\sa `CGAL::Triangulation_ds_vertex_base_3` 
\sa `CGAL::Triangulation_vertex_base_with_info_3` 

*/
template< typename TriangulationTraits_3, typename TriangulationDSVertexBase_3 >
class Triangulation_vertex_base_3 : public TriangulationDSVertexBase_3 {
public:

/// \name Types 
/// @{

/*! 

*/ 
typedef TriangulationTraits_3::Point_3 Point; 

/// @}

}; /* end Triangulation_vertex_base_3 */
} /* end namespace CGAL */
