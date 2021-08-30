namespace CGAL {

/*!
\ingroup PkgMesh3Domains

The class `Polyhedral_mesh_domain_with_features_3` implements a domain whose
boundary is a simplicial polyhedral surface.
This surface must be free of intersection.
It can either be closed,
included inside another polyhedral surface which is closed and free of intersection,
or open. In the latter case, the meshing process will only take care of the quality
of the 1D (features and boundaries) and 2D (surfaces) components of the mesh.

It is a model of the concept `MeshDomainWithFeatures_3`. It also
provides a member function to automatically detect sharp features and boundaries from
the input polyhedral surface(s).

\tparam IGT stands for a geometric traits class providing the types
and functors required to implement the intersection tests and intersection computations
for polyhedral boundary surfaces. This parameter has to be
instantiated with a model of the concept `IntersectionGeometricTraits_3`.

\cgalModels `MeshDomainWithFeatures_3`

\sa `CGAL::Mesh_domain_with_polyline_features_3<MD>`
\sa `CGAL::Polyhedral_mesh_domain_3<Polyhedron,IGT,TriangleAccessor>`
\sa `CGAL::Mesh_polyhedron_3<IGT>`
*/
template< typename IGT >
class Polyhedral_mesh_domain_with_features_3
  : public CGAL::Mesh_domain_with_polyline_features_3<
  CGAL::Polyhedral_mesh_domain_3< CGAL::Mesh_polyhedron_3<IGT>::type, IGT> >
 {
public:

/// \name Types
/// @{

/*!
Numerical type.
*/
typedef unspecified_type FT;

/// @}

/// \name Creation
/// @{

/*!
Constructs a `Polyhedral_mesh_domain_with_features_3` from a polyhedral surface of type `Polyhedron`.
The only requirement on type `Polyhedron` is that `CGAL::Mesh_polyhedron_3<IGT>::%type` should
be constructible from `Polyhedron`.
No feature detection is done at this level. Note that a copy of `bounding_polyhedron` will be done.
The polyhedron `bounding_polyhedron` has to be closed and free of intersections.
Its interior of `bounding_polyhedron` will be meshed.
*/
template <typename Polyhedron>
Polyhedral_mesh_domain_with_features_3(const Polyhedron& bounding_polyhedron);


/*!
Constructs a `Polyhedral_mesh_domain_with_features_3` from a polyhedral surface, and a bounding polyhedral surface.
`CGAL::Mesh_polyhedron_3<IGT>::%type` should be constructible from `Polyhedron`.
The first polyhedron should be entirely included inside `bounding_polyhedron`, which has to be closed
and free of intersections.
Using this constructor allows to mesh a polyhedral surface which is not closed, or has holes.
The inside of `bounding_polyhedron` will be meshed.
*/
template <typename Polyhedron>
Polyhedral_mesh_domain_with_features_3(const Polyhedron& polyhedron,
                                       const Polyhedron& bounding_polyhedron);

/*!
\deprecated Constructs a `Polyhedral_mesh_domain_with_features_3` from an off file. No feature
detection is done at this level. Users must read the file into a `Polyhedron_3` and call the
constructor above.
*/
Polyhedral_mesh_domain_with_features_3(const std::string& filename);

/// @}

/// \name Operations
/// @{

/*!
Detects sharp features and boundaries of the internal bounding polyhedron (and the potential internal polyhedron)
and inserts them as features of the domain. `angle_bound` gives the maximum
angle (in degrees) between the two normal vectors of adjacent triangles.
For an edge of the polyhedron, if the angle between the two normal vectors of its
incident facets is bigger than the given bound, then the edge is considered as
a feature edge.
*/
void detect_features(FT angle_bound=60);


/*!
Detects border edges of the bounding polyhedron and inserts them as features of the domain.
This function should be called alone only, and not before or after `detect_features()`.
*/
   void detect_borders();

/// @}

}; /* end Polyhedral_mesh_domain_with_features_3 */
} /* end namespace CGAL */
