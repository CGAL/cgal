
namespace CGAL {

/*!
\ingroup PkgTriangulation2VertexFaceClasses

The class `Regular_triangulation_vertex_base_2` is a model of the concept 
`RegularTriangulationVertexBase_2`. It is the default 
vertex base class of regular triangulations. 

\tparam Traits has to be a model 
of `RegularTriangulationTraits_2`. 

\tparam Vb has to be a model 
of the concept `TriangulationVertexBase_2` and is by default 
instantiated by `Triangulation_vertex_base_2<Traits>`. 

\cgalModels `RegularTriangulationVertexBase_2`

\sa `CGA::Triangulation_vertex_base_2<Traits,Vb>` 
\sa `CGAL::Regular_triangulation_2<Traits,Tds>` 
\sa `CGAL::Regular_triangulation_face_base_2<Traits,Fb>` 

*/
template< typename Traits, typename Vb >
class Regular_triangulation_vertex_base_2 : public Vb {
public:

/// @}

}; /* end Regular_triangulation_vertex_base_2 */
} /* end namespace CGAL */
