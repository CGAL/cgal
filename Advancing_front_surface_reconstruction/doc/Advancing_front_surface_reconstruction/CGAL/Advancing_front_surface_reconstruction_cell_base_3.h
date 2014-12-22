
namespace CGAL {

/*!
\ingroup PkgAdvancingFrontSurfaceReconstruction

The class `Advancing_front_surface_reconstruction_cell_base_3` is the default 
cell type for the class  `Advancing_front_surface_reconstruction`.

\tparam Traits has to be a model of `DelaunayTriangulationTraits_3`. 

\tparam Cb has to be a model of `TriangulationCellBase_3`.  

*/
template< typename Traits, typename Cb =  Triangulation_cell_base_3<Traits> >
class Advancing_front_surface_reconstruction_cell_base_3 : public Cb {
public:

/// @}

}; /* end Advancing_front_surface_reconstruction_cell_base_3 */
} /* end namespace CGAL */
