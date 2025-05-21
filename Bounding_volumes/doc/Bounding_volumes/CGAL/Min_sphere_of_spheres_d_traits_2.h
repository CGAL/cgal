namespace CGAL {

/*!
\ingroup PkgBoundingVolumesRef

The class
`Min_sphere_of_spheres_d_traits_2<K,FT,UseSqrt,Algorithm>` is a
model for concept `MinSphereOfSpheresTraits`. It uses a pair of \cgal
`Point_2` and `FT` to represent circles.

\cgalModels{MinSphereOfSpheresTraits}


The last two template parameters, `UseSqrt` and `Algorithm`, have
default arguments, namely `CGAL::Tag_false` and
`CGAL::Default_algorithm`, respectively. The template parameters
of class `Min_sphere_of_points_d_traits_2<K,FT,UseSqrt,Algorithm>`
must fulfill the following requirements:

\tparam K must be a model for `Kernel`.

\tparam FT is a number type, which fulfills the requirements of
type `FT` of concept `MinSphereOfSpheresTraits`: It must be
either `double` or `float`, or an exact number type.

\tparam UseSqrt fulfills the
requirements of type `Use_square_roots` of concept
`MinSphereOfSpheresTraits`: It must be either `Tag_true` or `Tag_false`,
and its default is `Tag_false`.

\tparam Algorithm fulfills the requirements of type `Algorithm` of
concept `MinSphereOfSpheresTraits`: It must be either
`Default_algorithm`, `LP_algorithm` or `Farthest_first_heuristic`,
and its default is `Default_algorithm`.
*/
template< typename K, typename FT, typename UseSqrt, typename Algorithm >
class Min_sphere_of_spheres_d_traits_2 {
public:
/// \name Constants
/// @{

/*!
is the constant 2, i.e.\ the dimension of \f$ \mathbb{R}^2\f$.
*/
static const int D;

/// @}

/// \name Types
/// In addition to the types required by the concept
/// `MinSphereOfSpheresTraits`, this model also defines the types
/// `Radius` and `Point`. Here's the complete list of defined types:
/// @{

/*!

*/
typedef unspecified_type FT;

/*!

*/
typedef unspecified_type Use_square_roots;

/*!

*/
typedef unspecified_type Algorithm;

/*!
is a typedef to the template parameter `FT`
*/
typedef unspecified_type Radius;

/*!
is a typedef to `K::Point_2`.
*/
typedef unspecified_type Point;

/*!
is a typedef to
`std::pair<Point,Radius>`.
*/
typedef unspecified_type Sphere;

/*!
is a typedef to
`K::Cartesian_const_iterator_2`.
*/
typedef unspecified_type Cartesian_const_iterator;

/// @}

/// \name Access Functions
/// The class provides the access functions required by the concept
/// `MinSphereOfSpheresTraits`; they simply map to the corresponding
/// routines of class `K::Point_2`:
/// @{

/*!
maps to `s.second`.
*/
FT radius(const
Sphere& s);

/*!
maps to `s.first.cartesian_begin()`.
*/
Cartesian_const_iterator center_cartesian_begin(const
Sphere& s);

/// @}

}; /* end Min_sphere_of_spheres_d_traits_2 */
} /* end namespace CGAL */
