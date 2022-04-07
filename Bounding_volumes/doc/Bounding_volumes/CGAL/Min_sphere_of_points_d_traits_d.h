namespace CGAL {

/*!
\ingroup PkgBoundingVolumesRef

The class
`Min_sphere_of_points_d_traits_d<K,FT,Dim,UseSqrt,Algorithm>` is
a model for concept `MinSphereOfSpheresTraits`. It uses the
\cgal type `Point_d` to represent circles.

\cgalModels `MinSphereOfSpheresTraits`


\tparam K is a model for `Kernel`.

\tparam FT is a number type, which fulfills the requirements of
type `FT` of concept `MinSphereOfSpheresTraits`: It must be
either `double` or `float`, or an exact number type.

\tparam UseSqrt fulfills the
requirements of type `Use_square_roots` of concept
`MinSphereOfSpheresTraits`: It must be either `Tag_true` or `Tag_false`,
and its default is  `Tag_false`.

\tparam Algorithm fulfills the requirements of type `Algorithm` of
concept `MinSphereOfSpheresTraits`: It must be either
`Default_algorithm`, `LP_algorithm` or `Farthest_first_heuristic`,
and its default is `Default_algorithm`.

*/
template< typename K, typename FT, typename Dim, typename UseSqrt, typename Algorithm >
class Min_sphere_of_points_d_traits_d {
public:

/// \name Constants
/// @{

/*!
is the constant `Dim`.
*/
typedef unspecified_type D;

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
is a typedef to `K::Point_d`.
*/
typedef unspecified_type Point;

/*!
is a typedef to
`Point`.
*/
typedef unspecified_type Sphere;

/*!
is a typedef to
`K::Cartesian_const_iterator_d`.
*/
typedef unspecified_type Cartesian_const_iterator;

/// @}

/// \name Access Functions
/// The class provides the access functions required by the concept
/// `MinSphereOfSpheresTraits`; they simply map to the corresponding
/// routines of class `K::Point_d`:
/// @{

/*!
returns `0`.
*/
FT radius(const
Sphere& s);

/*!
maps to `s.cartesian_begin()`.
*/
Cartesian_const_iterator center_cartesian_begin(const
Sphere& s);

/// @}

}; /* end Min_sphere_of_points_d_traits_d */
} /* end namespace CGAL */
