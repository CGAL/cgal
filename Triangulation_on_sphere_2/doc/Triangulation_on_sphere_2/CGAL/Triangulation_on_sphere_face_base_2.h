namespace CGAL {

/*!
\ingroup PkgTriangulationOnSphere2VertexFaceClasses

\cgalModels `TriangulationOnSphereFaceBase_2`

The class `Triangulation_on_sphere_face_base_2` is a model for the concept
`TriangulationOnSphereFaceBase_2`. It is the default face base class
for triangulations on the sphere.

`Triangulation_on_sphere_face_base_2` can be simply plugged into the triangulation data
structure of a triangulation, or used as a base class to derive
other base face classes tuned for specific applications.

\tparam Traits must be a geometric traits class.
The geometric traits is actually not used by the class.

\tparam Fb has to be a model of the concept `TriangulationDSFaceBase_2`
and will serve as a base class for `Triangulation_on_sphere_face_base_2`.
By default this parameter is instantiated by `Triangulation_ds_face_base_2<>`.

\sa `CGAL::Triangulation_ds_face_base_2<Tds>`
\sa `CGAL::Triangulation_face_base_with_info_2<Traits,Tds>`
\sa `CGAL::Triangulation_on_sphere_vertex_base_2<Traits,Vb>`
\sa `CGAL::Triangulation_on_sphere_2<Traits,Tds>`
*/
template< typename Traits, typename Fb >
class Triangulation_on_sphere_face_base_2
  : public Fb
{
public:

};

} // namespace CGAL
