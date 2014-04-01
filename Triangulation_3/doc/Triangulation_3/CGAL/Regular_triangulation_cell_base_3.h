
namespace CGAL {

/*!
\ingroup PkgTriangulation3VertexCellClasses

The class `Regular_triangulation_cell_base_with_weighted_circumcenter_3` derives from 
`Cb`, a cell base class of a 3D triangulation.
It is the default cell base class of regular triangulations.

\tparam RegularTriangulationTraits_3 is the geometric traits class.

\tparam Cb is a cell base class from which `Regular_triangulation_cell_base_3`
derives. It must be a model of `TriangulationCellBase_3`. 
By default, this parameter is instantiated by 
`Triangulation_cell_base_3<RegularTriangulationTraits_3>`. 

\cgalModels `RegularTriangulationCellBase_3`

\sa `RegularTriangulationCellBase_3` 
\sa `CGAL::Regular_triangulation_3` 
\sa `CGAL::Regular_triangulation_cell_base_with_weighted_circumcenter_3`

*/

template< typename RegularTriangulationTraits_3, typename Cb >
class Regular_triangulation_cell_base_3 : public Cb {
public:

/// \name Types 
/// @{
typedef RegularTriangulationTraits_3::Bare_point Bare_point;
/// @}

/*! \name Access function 

As a model of the concept `RegularTriangulationCellBase_3`, 
`Regular_triangulation_cell_base_3` 
provides a `weighted_circumcenter()` member fonction. 

In this model, the `weighted_circumcenter()` member fonction returns the 
<b>weighted circumcenter</b> of the cell, computed 
by the `ConstructWeightedCircumcenter` constructor of the traits class.

Note that this point has no weight.
*/

/// @{
/*! 
Returns the weighted circumcenter of the cell.
Be careful that the return type is `RegularTriangulationTraits_3::Bare_point`,
and the radius of the weighted 
circumcenter is not supposed to be computed
by the constructor `Construct_weighted_circumcenter_3` of the traits
class, so the returned point has no weight.
*/ 
const Bare_point& weighted_circumcenter(
  const RegularTriangulationTraits_3& gt = RegularTriangulationTraits_3()) const; 

/// @}

}; /* end Regular_triangulation_cell_base_3 */
} /* end namespace CGAL */
