
namespace CGAL {

/*!
\ingroup PkgMesh2Ref



The class `Delaunay_mesh_criteria_2` is a model for the `MeshingCriteria_2` concept.
The shape criterion on triangles is given by a bound \f$ B\f$ such that for good
triangles \f$ \frac{r}{l} \le B\f$ where \f$ l\f$ is the shortest edge length
and \f$ r\f$ is the circumradius of the triangle. By default, \f$ B=\sqrt{2}\f$,
which is the best bound one can use with the guarantee that the refinement
algorithm will terminate. The upper bound \f$ B\f$ is related to a lower bound
\f$ \alpha_{min}\f$ on the minimum angle in the triangle:
\f[
\sin{ \alpha_{min} } = \frac{1}{2 B}
\f]
so \f$ B=\sqrt{2}\f$ corresponds to \f$ \alpha_{min} \ge 20.7\f$ degrees.

\tparam CDT must be a 2D constrained Delaunay triangulation.

\cgalModels{MeshingCriteria_2}


*/
template< typename CDT >
class Delaunay_mesh_criteria_2 {
public:

/// \name Creation
/// @{

/*!
%Default constructor with bound \f$ B=\sqrt{2}\f$.
*/
Delaunay_mesh_criteria_2();





/*!
Construct
a traits class with bound \f$ B=\sqrt{\frac{1}{4 b}}\f$.
*/
Delaunay_mesh_criteria_2(double b = 0.125);





/// @}

}; /* end Delaunay_mesh_criteria_2 */
} /* end namespace CGAL */
