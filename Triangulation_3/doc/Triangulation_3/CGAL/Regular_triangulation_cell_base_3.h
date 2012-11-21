
namespace CGAL {

/*!
\ingroup PkgTriangulation3VertexCellClasses

The class `Regular_triangulation_cell_base_3` is a model of the concept 
`RegularTriangulationCellBase_3`. It is the default face base class 
of regular triangulations. 


\tparam Traits must be a model of `RegularTriangulationTraits_3`. 

\tparam Cb must be a model of `TriangulationCellBase_3`. 
By default, this parameter is instantiated by `Triangulation_cell_base_3<Traits>`. 

\cgalModels ::RegularTriangulationCellBase_3 

\sa `RegularTriangulationCellBase_3` 
\sa `RegularTriangulationTraits_3` 
\sa `CGAL::Regular_triangulation_3<Traits,Tds>` 

*/
template< typename Traits, typename Cb >
class Regular_triangulation_cell_base_3 : public Cb {
}; /* end Regular_triangulation_cell_base_3 */
} /* end namespace CGAL */
