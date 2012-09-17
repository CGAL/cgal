
namespace CGAL {

/*!
\ingroup PkgTriangulation2VertexFaceClasses

The class `Triangulation_vertex_base_2` is the default model for the concept 
`TriangulationVertexBase_2`. 

`Triangulation_vertex_base_2` can be simply plugged in the triangulation data 
structure 
of a triangulation, or used as a base class to derive 
other base vertex classes tuned for specific applications. 

### Parameters ###

`Triangulation_vertex_base_2` is templated by a geometric traits class which provide the type 
`Point`. It is strongly recommended to instantiate this 
traits class with the model used for the triangulation traits class. 
This ensures that the point type defined by `Triangulation_vertex_base_2` is the same as the point type defined by 
the triangulation. 

The second template parameter of `Triangulation_vertex_base_2` 
has to be a model of the concept `TriangulationDSVertexBase_2` 
By default this parameter is 
instantiated by `CGAL::Triangulation_ds_vertex_base_2<>`. 

\models ::TriangulationVertexBase_2 

\sa `CGAL::Triangulation_ds_vertex_base_2<Tds>` 
\sa `CGAL::Triangulation_face_base_2<Traits,Fb>` 
\sa `CGAL::Regular_triangulation_vertex_base_2<Traits,Vb>` 
\sa `CGAL::Triangulation_vertex_base_with_info_2<Info,Traits,Vb>` 

*/
template< typename Traits, typename Vb >
class Triangulation_vertex_base_2 : public Vb {
public:

/// @}

}; /* end Triangulation_vertex_base_2 */
} /* end namespace CGAL */
