
namespace CGAL {

/*!
\ingroup PkgTDS2Ref

The class `Triangulation_ds_vertex_base_2` can be used as the base vertex
for a `Triangulation_data_structure_2`, it is a model of the concept
`TriangulationDSVertexBase_2`.

This base class can be used directly or can serve as a base to derive
other base classes with some additional attributes (a color for
example) tuned for a specific application.

Note that if the `Triangulation_data_structure_2`
is used as a parameter of a
geometric triangulation, there are additional geometric requirements
to be fulfilled by the vertex base class,
and `Triangulation_ds_vertex_base_2` cannot be plugged in.

\tparam TDS should not be specified (see Section \ref TDS_2TheRebindMechanism and examples)

\cgalModels `TriangulationDSVertexBase_2`

\sa `CGAL::Triangulation_vertex_base_2<Traits,Vb>`
\sa `CGAL::Triangulation_ds_face_base_2<TDS>`

*/
template< typename TDS = void >
class Triangulation_ds_vertex_base_2 {
public:

}; /* end Triangulation_ds_vertex_base_2 */
} /* end namespace CGAL */
