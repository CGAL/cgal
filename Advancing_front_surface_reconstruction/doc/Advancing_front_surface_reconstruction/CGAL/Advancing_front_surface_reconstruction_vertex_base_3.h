
namespace CGAL {

/*!
\ingroup PkgAdvancingFrontSurfaceReconstruction

The class `Advancing_front_surface_reconstruction_vertex_base_3` is the default
vertex type for the class  `Advancing_front_surface_reconstruction`.

\tparam Traits has to be a model of `DelaunayTriangulationTraits_3`. 

\tparam Vb has to be a model of `TriangulationVertexBase_3`. 


*/
  template< typename Traits, typename Vb = Triangulation_vertex_base_3<Traits> >
class Advancing_front_surface_reconstruction_vertex_base_3 : public Vb {
public:

/// @}

}; /* end Advancing_front_surface_reconstruction_vertex_base_3 */
} /* end namespace CGAL */
