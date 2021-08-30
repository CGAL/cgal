
namespace CGAL {

/*!
\ingroup PkgTriangulation2VertexFaceClasses

The class `Constrained_triangulation_face_base_2` is the default model for the concept
`ConstrainedTriangulationFaceBase_2` to be used as base face class
of constrained triangulations.

\cgalModels `ConstrainedTriangulationFaceBase_2`


\tparam Traits must be a geometric traits.

\tparam Fb must be a model
of the concept `TriangulationFaceBase_2`.
Its default is `Triangulation_face_base_2<Traits>`

The class `Constrained_triangulation_face_base_2` derives from `Fb`
and adds three Boolean to deal with information about constrained edges.

The member functions `cw(int i)`, `ccw(int i)`
and `reorient` are overloaded to update
information about constrained edges.

\sa `TriangulationFaceBase_2`
\sa `ConstrainedTriangulationFaceBase_2`
\sa `CGAL::Constrained_triangulation_2<Traits,Tds>`
\sa `CGAL::Triangulation_face_base_2<Traits>`

*/
template< typename Traits, typename Fb >
class Constrained_triangulation_face_base_2 : public Fb {
public:

}; /* end Constrained_triangulation_face_base_2 */
} /* end namespace CGAL */
