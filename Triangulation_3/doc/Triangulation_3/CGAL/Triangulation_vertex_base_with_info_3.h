
namespace CGAL {

/*!
\ingroup PkgTriangulation3

The class `Triangulation_vertex_base_with_info_3` is a model of the concept 
`TriangulationVertexBase_3`, the base vertex of a 3D-triangulation. 
It provides an easy way to add some user defined information in vertices. 
Note that input/output operators discard this additional information. 

Parameters 
-------------- 

The first template argument is the information the user would like to add 
to a vertex. It has to be `DefaultConstructible` and `Assignable`. 

The second template argument is the geometric traits class 
`TriangulationTraits_3` which provides the `Point_3`. 

The third template argument is a vertex base class from which 
`Triangulation_vertex_base_with_info_3` derives. It has the default 
value `Triangulation_vertex_base_3<TriangulationTraits_3>`. 

\models ::TriangulationVertexBase_3 
\models ::TriangulationVertexBaseWithInfo_3 

\sa `CGAL::Triangulation_cell_base_with_info_3` 
\sa `CGAL::Triangulation_vertex_base_3` 

*/
template< typename Info, typename TriangulationTraits_3, typename TriangulationVertexBase_3 >
class Triangulation_vertex_base_with_info_3 : public TriangulationVertexBase_3 {
public:

/// \name Types 
/// @{

/*! 

*/ 
typedef Info Info; 

/// @} 

/// \name Access Functions 
/// @{

/*! 
Returns a const reference to the object of type `Info` stored in the 
vertex. 
*/ 
const Info& info() const; 

/*! 
Returns a reference to the object of type `Info` stored in the vertex. 
*/ 
Info& info(); 

/// @}

}; /* end Triangulation_vertex_base_with_info_3 */
} /* end namespace CGAL */
