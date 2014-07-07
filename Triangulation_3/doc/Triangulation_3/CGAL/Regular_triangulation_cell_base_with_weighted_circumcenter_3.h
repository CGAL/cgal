
namespace CGAL {

/*!
\ingroup PkgTriangulation3VertexCellClasses

The class `Regular_triangulation_cell_base_with_weighted_circumcenter_3` derives from 
`Cb`, a cell base class of a 3D triangulation.
It provides an easy way to cache the computation of the weighted circumcenters of 
tetrahedra. 
Note that input/output operators discard this additional information. 

All functions modifying the vertices of the cell invalidate the cached 
circumcenter. 

\tparam RegularTriangulationTraits_3 is the geometric traits class.

\tparam Cb is a cell base class from which 
`Regular_triangulation_cell_base_with_weighted_circumcenter_3` derives. Cb should
be a model of `RegularTriangulationCellBase_3`.
It has the default value `Regular_triangulation_cell_base_3<RegularTriangulationTraits_3>`.

\cgalModels `RegularTriangulationCellBase_3`

\sa `CGAL::Triangulation_cell_base_3` 
\sa `CGAL::Triangulation_cell_base_with_info_3` 
\sa `CGAL::Regular_triangulation_cell_base_3`

*/
template< typename RegularTriangulationTraits_3, typename Cb >
class Regular_triangulation_cell_base_with_weighted_circumcenter_3 : public Cb {
public:
	
/// \name Types 
/// @{
typedef RegularTriangulationTraits_3::Bare_point Bare_point;
/// @}

/*! \name Access function 

As a model of the concept `RegularTriangulationCellBase_3`, 
`Regular_triangulation_cell_base_3` 
provides a `weighted_circumcenter()` member fonction. 

In this model, the `weighted_circumcenter()` member fonction returns the <b>weighted circumcenter</b>
of the cell.
This `Bare_point` is computed by the `Construct_weighted_circumcenter_3` constructor of the traits class
when this function is first called.
In the next calls, the cached value is returned.

Note that this point has no weight.
*/

/// @{

/*!
Computes the weighted circumcenter of the tetrahedron, or retrieves it if already 
computed.

The returned point has no weight.
*/ 
const Bare_point& weighted_circumcenter( 
	const RegularTriangulationTraits_3&gt = RegularTriangulationTraits_3()) const; 

/*!
Swaps the Regular_triangulation_cell_base_with_weighted_circumcenter_3 and other.
Should be preferred to an assignment or copy constructor when other is deleted after that.
*/
void swap (Regular_triangulation_cell_base_with_weighted_circumcenter_3& other) throw();

/// @}

}; /* end Regular_triangulation_cell_base_with_weighted_circumcenter_3 */
} /* end namespace CGAL */
