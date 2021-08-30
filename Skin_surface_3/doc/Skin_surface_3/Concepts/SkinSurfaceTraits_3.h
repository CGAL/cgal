
/*!
\ingroup PkgSkinSurface3Concepts
\cgalConcept

Required types and member functions for the `SkinSurfaceTraits_3` concept.
This geometric traits concept is used for the construction of a
polyhedral mesh approximating a skin surface
`CGAL::Skin_surface_3`.

\cgalRefines `RegularTriangulationTraits_3`

\cgalHasModel `CGAL::Skin_surface_traits_3<K>`

\sa `CGAL::Skin_surface_3<SkinSurfaceTraits_3>`
\sa `CGAL::Skin_surface_traits_3<K>`

*/

class SkinSurfaceTraits_3 {
public:

/// \name Types
/// @{

/*!
The number type.
*/
typedef unspecified_type RT;

/*!
A constructor object which
provides the following function operator:

`Point operator()(const Point &center_del, const
Point &center_vor);`

Constructs the anchor point in between the Delaunay and Voronoi
centers, i.e., the point \f$ (1-s)\cdot\f$`center_del` +
\f$ s\cdot\f$`center_vor`, where \f$ s\f$ is the shrink factor.
*/
typedef unspecified_type Construct_anchor_point_3;

/// @}

/// \name Creation
/// @{

/*!
Constructor that takes the shrink factor as argument.
For meshing the boundary of the union
of a set of balls, the shrink factor is discarded.
*/
SkinSurfaceTraits_3(RT s=.5);

/// @}

/// \name Operations
/// The following functions give access to the constructor objects:
/// @{

/*!
Returns the shrink factor.
*/
Regular_RT shrink_factor() const;

/*!
Returns a `Construct_anchor_point_3` object.
*/
Construct_anchor_point_3 construct_anchor_point_3_object()
const;

/// @}

}; /* end SkinSurfaceTraits_3 */

