namespace CGAL {

/*!
\ingroup PkgTriangulationOnSphere2VertexFaceClasses

\cgalModels `TriangulationOnSphereVertexBase_2`

The class `Triangulation_on_sphere_vertex_base_2` is the default model
for the concept `TriangulationOnSphereVertexBase_2` and
the default vertex base class for triangulations on the sphere.

`Triangulation_on_sphere_vertex_base_2` can be simply plugged into the triangulation data
structure of a triangulation, or used as a base class to derive
other base vertex classes tuned for specific applications.

\tparam Traits must be a geometric traits class which provides the type `Point_on_sphere_2`.
It is strongly recommended to instantiate this traits class with the model used
for the triangulation traits class. This ensures that the point type defined by
`Triangulation_on_sphere_vertex_base_2` is the same as the point type defined by
the triangulation.

\tparam Vb must be a model of the concept `TriangulationDSVertexBase_2`.
By default this parameter is instantiated by `Triangulation_ds_vertex_base_2<>`.

\sa `CGAL::Triangulation_ds_vertex_base_2<Tds>`
\sa `CGAL::Triangulation_on_sphere_face_base_2<Traits,Fb>`
\sa `CGAL::Triangulation_vertex_base_with_info_2<Info,Traits,Vb>`
\sa `CGAL::Triangulation_on_sphere_2<Traits,Tds>`
*/
template< typename Traits, typename Vb >
class Triangulation_on_sphere_vertex_base_2
  : public Vb
{
public:

  /// \name Types
  /// @{

  /*!

  */
  typedef Traits::Point_on_sphere_2 Point;

  /// @}
};

} // namespace CGAL
