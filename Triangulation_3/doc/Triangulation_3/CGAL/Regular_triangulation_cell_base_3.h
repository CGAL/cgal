
namespace CGAL {

/*!
\ingroup PkgTriangulation3VertexCellClasses

The class `Regular_triangulation_cell_base_3` is a model of the concept 
`RegularTriangulationCellBase_3`. It is the default cell base class 
of regular triangulations. 


\tparam Traits must be a model of `RegularTriangulationTraits_3`. 

\tparam Cb must be a model of `TriangulationCellBase_3`. 
By default, this parameter is instantiated by `Triangulation_cell_base_3<Traits>`. 

\cgalModels `RegularTriangulationCellBase_3`

\sa `RegularTriangulationCellBase_3` 
\sa `RegularTriangulationTraits_3` 
\sa `CGAL::Regular_triangulation_3<Traits,Tds>` 
\sa `CGAL::Triangulation_cell_base_with_circumcenter_3<RegularTriangulationTraits_3, CellBase_3>`

*/

template< typename Traits, typename Cb >
class Regular_triangulation_cell_base_3 : public Cb {
public:

/*! \name Access function 

As a model of the concept `RegularTriangulationCellBase_3`, 
`Regular_triangulation_cell_base_3` 
provides a `circumcenter()` member fonction. 

In this model, the `circumcenter()` member fonction returns the <b>weighted circumcenter</b>
of the cell, computed by the `ConstructWeightedCircumcenter` constructor of the traits class. 
However, this point has no weight.

A class for the cells of regular triangulations with
cached weighted circumcenters can simply be obtained by plugging 
`Regular_triangulation_cell_base_3` in the second template parameter of 
`CGAL::Triangulation_cell_base_with_circumcenter_3<RegularTriangulationTraits_3, CellBase_3>`
*/

/// @{
/*! 
Returns the weighted circumcenter of the cell.
Be careful that the return type is `Traits::Bare_point`, and the radius of the  weighted 
circumcenter is not supposed to be computed
by the constructor `ConstructWeightedCircumcenter` of the traits
class, so the returned point has no weight.
*/ 
const Traits::Bare_point& circumcenter(const Traits& gt = Traits()) const; 

/// @}
	
}; /* end Regular_triangulation_cell_base_3 */
} /* end namespace CGAL */
