
namespace CGAL {

/*!
\ingroup PkgTriangulation2VertexFaceClasses

The class `Regular_triangulation_face_base_2` is a model of the concept 
`RegularTriangulationFaceBase_2`. It is the default face base class 
of regular triangulations. 

\tparam Traits has to be a model 
of `RegularTriangulationTraits_2`. 

\tparam Fb has to be a model 
of `TriangulationFaceBase_2`. By default, this parameter is 
instantiated by 
`Triangulation_face_base_2<Traits>`. 

\cgalModels `RegularTriangulationFaceBase_2`

\sa `RegularTriangulationFaceBase_2` 
\sa `RegularTriangulationTraits_2` 
\sa `CGAL::Regular_triangulation_2<Traits,Tds>` 
\sa `CGAL::Regular_triangulation_vertex_base_2<Traits,Vb>` 

*/
template< typename Traits, typename Fb >
class Regular_triangulation_face_base_2 : public Fb {
public:

/// @}

}; /* end Regular_triangulation_face_base_2 */
} /* end namespace CGAL */
