
namespace CGAL {

/*!
\ingroup PkgTriangulation3VertexCellClasses

The class `Regular_triangulation_cell_base_with_weighted_circumcenter_3` derives from 
`Cb`, a cell base class of a 3D triangulation.
It is the default cell base class of regular triangulations.

\tparam Traits is the geometric traits class. It must be a model of `RegularTriangulationTraits_3`.

\tparam Cb is a cell base class from which `Regular_triangulation_cell_base_3`
derives. It must be a model of `TriangulationCellBase_3`. 
By default, this parameter is instantiated by 
`Triangulation_cell_base_3<Traits>`.

\cgalModels `RegularTriangulationCellBase_3`

\sa `RegularTriangulationCellBase_3` 
\sa `CGAL::Regular_triangulation_3` 
\sa `CGAL::Regular_triangulation_cell_base_with_weighted_circumcenter_3`

*/

template< typename Traits, typename Cb >
class Regular_triangulation_cell_base_3
  : public Cb
{
public:

/// \name Types 
/// @{
typedef Traits::Point_3 Point_3;

typedef Traits::Weighted_point_3 Point;

typedef std::list<Point> Point_container;

typedef Point_container::iterator Point_iterator;

/// @}

/// \name Hidden points-related functions
/// Not every weighted point inserted in a regular triangulation necessarily
/// appears in the trinagulation. If the weight of a point is too small compared
/// to other points, it might be <I>hidden</I>. These hidden vertices are stored
/// in a unique corresponding cell (defined through v->cell()). The following
/// functions provide set and get functionalities.
/// @{

/*!
Returns an iterator pointing to the first hidden point.
*/
Point_iterator hidden_points_begin();

/*!
Returns a past-the-end iterator.
*/
Point_iterator hidden_points_end();

/*!
Adds `p` to the set of hidden points of the cell.
*/
void hide_point(const Point & p);

/// @}

/// @{
/*! 
Returns the weighted circumcenter of the cell.
Be careful that the return type is `Point_3`,
and the radius of the weighted 
circumcenter is not supposed to be computed
by the constructor `Construct_weighted_circumcenter_3` of the traits
class, hence the returned point has no weight.
*/ 
const Point_3& weighted_circumcenter(const Traits& gt = Traits()) const;

/// @}

}; /* end Regular_triangulation_cell_base_3 */
} /* end namespace CGAL */
