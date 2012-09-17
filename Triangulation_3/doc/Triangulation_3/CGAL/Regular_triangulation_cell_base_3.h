
namespace CGAL {

/*!
\ingroup PkgTriangulation3VertexFaceClasses

The class `Regular_triangulation_cell_base_3` is a model of the concept 
`RegularTriangulationCellBase_3`. It is the default face base class 
of regular triangulations. 

### Parameters ###

The template parameters `Traits` has to be a model 
of `RegularTriangulationTraits_3`. 

The template parameter `Cb` has to be a model 
of `TriangulationCellBase_3`. By default, this parameter is 
instantiated by 
`CGAL::Triangulation_cell_base_3<Traits>`. 

\models ::RegularTriangulationCellBase_3 

\sa `RegularTriangulationCellBase_3` 
\sa `RegularTriangulationTraits_3` 
\sa `CGAL::Regular_triangulation_3<Traits,Tds>` 

*/
template< typename Traits, typename Cb >
class Regular_triangulation_cell_base_3 : public Cb {
}; /* end Regular_triangulation_cell_base_3 */
} /* end namespace CGAL */
