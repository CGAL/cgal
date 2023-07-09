
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

\tparam Traits is the geometric traits class and must be a model of `RegularTriangulationTraits_3`.

\tparam Cb is a cell base class from which
`Regular_triangulation_cell_base_with_weighted_circumcenter_3` derives. Cb should
be a model of `RegularTriangulationCellBase_3`.
It has the default value `Regular_triangulation_cell_base_3<RT>`.

\cgalModels{RegularTriangulationCellBaseWithWeightedCircumcenter_3}

\sa `CGAL::Triangulation_cell_base_3`
\sa `CGAL::Triangulation_cell_base_with_info_3`
\sa `CGAL::Regular_triangulation_cell_base_3`

*/
template< typename Traits, typename Cb >
class Regular_triangulation_cell_base_with_weighted_circumcenter_3 : public Cb {
public:

/// \name Types
/// @{
typedef Traits::Point_3 Point_3;

typedef Traits::Weighted_point_3 Point;
/// @}

/*! \name Access function

As a model of the concept `RegularTriangulationCellBase_3`,
`Regular_triangulation_cell_base_with_weighted_circumcenter_3`
provides a `weighted_circumcenter()` member function.

In this model, the `weighted_circumcenter()` member function returns the <b>weighted circumcenter</b>
of the cell.
This `Point_3` is computed using the `Construct_weighted_circumcenter_3` functor
of the traits class when this function is first called and its value is stored.
In the next calls, the cached value is returned.

Note that the returned point has no weight.
*/

/// @{

/*!
Computes the weighted circumcenter of the tetrahedron, or retrieves it if it has
already been computed.

The returned point has no weight.
*/
const Point_3& weighted_circumcenter(const Traits&gt = Traits()) const;

/*!
Swaps the Regular_triangulation_cell_base_with_weighted_circumcenter_3 and `other`.
This function should be preferred to an assignment or the copy constructor
if `other` is deleted thereafter.
*/
void swap (Regular_triangulation_cell_base_with_weighted_circumcenter_3& other) noexcept;

/// @}

}; /* end Regular_triangulation_cell_base_with_weighted_circumcenter_3 */
} /* end namespace CGAL */
