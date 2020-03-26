
namespace CGAL {

/*!
\ingroup PkgTriangulation2VertexFaceClasses

The class `Triangulation_face_base_2` is a model for the concept
`TriangulationFaceBase_2`. It is the default face base class
for basic and Delaunay triangulation.

These default base class can be used directly or can serve as a base to derive
other base classes with some additional attribute (a color for example)
tuned for specific applications.

\tparam Traits must be a geometric traits class.
The geometric traits is actually not used by the class.

\tparam Fb has to be a model of the concept `TriangulationDSFaceBase_2`
and will serve as a base class for `Triangulation_face_base_2` .
\cgal provides a default instantiation for this parameter which is
`Triangulation_ds_face_base_2<>`.

\cgalModels `TriangulationFaceBase_2`

\sa `CGAL::Triangulation_ds_face_base_2<Tds>`
\sa `CGAL::Triangulation_vertex_base_2<Traits,Vb>`
\sa `CGAL::Triangulation_2<Traits,Tds>`

*/
template< typename Traits, typename Fb >
class Triangulation_face_base_2 : public Fb {
public:

}; /* end Triangulation_face_base_2 */
} /* end namespace CGAL */
