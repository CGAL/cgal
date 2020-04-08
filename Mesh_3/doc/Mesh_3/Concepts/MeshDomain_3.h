/*!
\ingroup PkgMesh3Concepts
\cgalConcept

The concept `MeshDomain_3` describes the knowledge required on the
object to be discretized.
The concept `MeshDomain_3` is the concept to be used when the input
domain does not have \f$ 0\f$ or \f$ 1\f$-dimensional features that need to be
accurately approximated by the mesh.
In such a case, the queries issued by the meshing process concern only the faces of the input domain with dimension \f$ 3\f$
and \f$ 2\f$, that are respectively called <I>subdomains</I> and <I>surface patches</I>.

More specifically the concept `MeshDomain_3` provides
a method to localize a point with respect to the input domain
and its subdomains.
Moreover, as the concept `SurfaceMeshTraits_3`,
it also provides
predicates to test whether a query segment
(or a ray, or a line) intersects the boundary of the domain or
of the subdomains, and constructors
to compute some intersection point if any.
It also includes a method able to provide
a small set of initial points on the boundary.

In the following we consider only proper intersection with the domain and
subdomain boundaries.
A segment, ray or line is said to intersect properly the domain boundary
(resp. a subdomain boundary)
if it includes points which are strictly inside
and strictly outside the domain (resp. the subdomain).

\cgalHasModel `CGAL::Polyhedral_mesh_domain_3<Polyhedron,IGT,TriangleAccessor>`
\cgalHasModel `CGAL::Labeled_mesh_domain_3<BGT>`

\sa `MeshVertexBase_3`
\sa `MeshCellBase_3`
\sa `CGAL::make_mesh_3()`
\sa `CGAL::refine_mesh_3()`

*/

class MeshDomain_3 {
public:

/// \name Types
/// @{

/*!
Geometric traits class. This type is defined to ensure compatibility with
`CGAL::Kernel_traits<T>`.
*/
typedef unspecified_type R;

/*!
Point type.
*/
typedef unspecified_type Point_3;

/*!
Segment type.
*/
typedef unspecified_type Segment_3;

/*!
Ray type.
*/
typedef unspecified_type Ray_3;

/*!
Line type.
*/
typedef unspecified_type Line_3;

/*!
A type to distinguish
`MeshDomain_3` models from `MeshDomainWithFeatures_3` models.
*/
typedef CGAL::Tag_false Has_features;

/*!
Type of indices for subdomains of the
input domain. Must be a model of CopyConstructible,
Assignable, DefaultConstructible and EqualityComparable.
The default constructed value must match the label of the exterior of
the domain (which contains at least the unbounded component).
*/
typedef unspecified_type Subdomain_index;

/*!
Type of indices for surface patches
(boundaries and interfaces) of the
input domain. Must be a model of CopyConstructible,
Assignable, DefaultConstructible and EqualityComparable.
The default constructed value must be the index value assigned
to a non surface facet.
*/
typedef unspecified_type Surface_patch_index;

/*!
Type of indices to be stored at mesh vertices
to characterize the lowest dimensional face of the input complex
on which the vertex lies. Must be a model of CopyConstructible,
Assignable, DefaultConstructible and EqualityComparable.

*/
typedef unspecified_type Index;

/*!
Returns type of `Construct_intersection` queries.
`int` represents the
dimension of the lower dimensional face of the input complex on which the intersection
point lies and `%Index` is the index of this face.
*/
typedef std::tuple<Point_3, Index, int> Intersection;

/*!
A function object to construct
a set of initial points on the surface of the domain. Provides the
following operators:

`template<typename OutputIterator>`
<br>
`OutputIterator operator()(OutputIterator pts)`

`template<typename OutputIterator>`
<br>
`OutputIterator operator()(int n, OutputIterator pts)`

Those two operators output a set of (`n`) surface points to the
output iterator `pts`, as objects of type `std::pair<Point_3,
%Index>`. If `n` is not given, the functor must provide enough
points to initialize the mesh generation process.
*/
typedef unspecified_type Construct_initial_points;

/*!
A function object to query whether a point is in
the input domain or not. In the positive case, it outputs the
subdomain which includes the query point. Provides the operator:

`boost::optional<Subdomain_index> operator()(Point_3 p)`
*/
typedef unspecified_type Is_in_domain;

/*!
A function object which answers
intersection queries between the surface patches of the domain and
objects of type `Segment_3`, `Ray_3` or
`Line_3`. Provides the operators:

`boost::optional<Surface_patch_index> operator()(Segment_3 s)`

`boost::optional<Surface_patch_index> operator()(Ray_3 r)`

`boost::optional<Surface_patch_index> operator()(Line_3 l)`

The return type of the operators tell whether or not the query intersects a
surface patch. In the positive case, it provides (through operator*()) the
`Surface_patch_index` of one of the intersected surface patches.
*/
typedef unspecified_type Do_intersect_surface;

/*!
A function object to construct the
intersection between an object of type `Segment_3`, `Ray_3` or
`Line_3` and an interface. Provides the operators:

`Intersection operator()(Segment_3 s)`

`Intersection operator()(Ray_3 r)`

`Intersection operator()(Line_3 l)`
\pre do_intersect_surface(s/r/l) == true
*/
typedef unspecified_type Construct_intersection;

/// @}

/// \name Bounding box
/// Since CGAL-4.8, a model of `MeshDomain_3` must provide a function
/// providing a bounding box of the domain.
/// @{

/// Returns a bounding box of the domain
Bbox_3 bbox() const;
/// @}

/// \name Operations
/// The following functions give access to the function objects:
/// @{

/*!

*/
Construct_initial_points construct_initial_points_object();

/*!

*/
Is_in_domain is_in_domain_object();

/*!

*/
Do_intersect_surface do_intersect_surface_object();

/*!

*/
Construct_intersection construct_intersection_object();

/// @}

/// \name Index Conversion
/// These methods are designed to convert indices:
/// @{

/*!
Returns
the index to be stored at a vertex lying on the surface patch identified by `surface_patch_index`.
*/
Index index_from_surface_patch_index(Surface_patch_index surface_patch_index);

/*!
Returns
the index to be stored at a vertex lying in the subdomain identified by `subdomain_index`.
*/
Index index_from_subdomain_index(Subdomain_index subdomain_index);

/*!
Returns the `Surface_patch_index` of the surface patch
where lies a vertex with dimension 2 and index `index`.
*/
Surface_patch_index surface_patch_index(Index index);

/*!
Returns the index
of the subdomain containing a vertex with dimension 3 and index `index`.
*/
Subdomain_index subdomain_index(Index index);

/// @}

}; /* end MeshDomain_3 */
