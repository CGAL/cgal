
namespace CGAL {

/*!
\ingroup PkgTriangulation2VertexFaceClasses

The class `Triangulation_vertex_base_with_info_2` is designed to be used as a base vertex class
of a triangulation. It provides an easy way to add some user defined information in vertices.

\tparam Info is the information the user would like to add
to a vertex. It has to be `DefaultConstructible` and `Assignable`.

\tparam Traits is the geometric traits class
which provides the `Point_2`. It is strongly
recommended to instantiate this parameter
with the traits class used for the triangulation.
This ensures that the point type defined by `Triangulation_vertex_base_with_info_2` matches the point type defined by
the triangulation.

\tparam Vb must be a vertex base class from which
`Triangulation_vertex_base_with_info_3` derives. By default
this parameter is instantiated by
`Triangulation_vertex_base_2<Traits>`.

\cgalModelsBareBegin
\cgalModelsBare{`TriangulationVertexBaseWithInfo_2`}
\cgalModelsBare{The parameter `Vb` is a model of some vertex base concept.
`Triangulation_vertex_base_with_info_2` derives from `Vb` and will be a model of the
same vertex base concept: `TriangulationVertexBase_2`, or `RegularTriangulationVertexBase_2`.}
\cgalModelsBareEnd

\sa `CGAL::Triangulation_face_base_with_info_2<Info,Traits,Fb>`
\sa `CGAL::Triangulation_vertex_base_2<Traits,Vb>`
\sa `CGAL::Regular_triangulation_vertex_base_2<Traits,Vb>`

*/
template< typename Info, typename Traits, typename Vb >
class Triangulation_vertex_base_with_info_2 {
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

}; /* end Triangulation_vertex_base_with_info_2 */
} /* end namespace CGAL */
