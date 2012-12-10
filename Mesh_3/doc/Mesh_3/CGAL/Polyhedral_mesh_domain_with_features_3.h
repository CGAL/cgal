namespace CGAL {

/*!
\ingroup PkgMesh_3Domains

The class `Polyhedral_mesh_domain_with_features_3` implements a domain whose 
boundary is a simplicial polyhedral surface. 
This surface must be closed and free of intersection. 
It is a model of the concept `MeshDomainWithFeatures_3`. It also 
provides a member function to automatically detect sharp features from 
the input polyhedral surface. 


\tparam IGT stands for a geometric traits class providing the types 
and functors required to implement the intersection tests and intersection computations 
for polyhedral boundary surfaces. This parameter has to be 
instantiated with a model of the concept `IntersectionGeometricTraits_3`. 

\cgalModels `MeshDomainWithFeatures_3` 

\sa `CGAL::Mesh_domain_with_polyline_features_3<MeshDomain>` 
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
typedef Hidden_type FT; 

/// @} 

/// \name Creation 
/// @{

/*! 

Constructs a `Polyhedral_mesh_domain_with_features_3` from a `Polyhedron`. 
The only requirement on type `Polyhedron` is that `CGAL::Mesh_polyhedron_3<IGT>::type` should 
be constructible from `Polyhedron`. 
No feature detection is done at this level. Note that a copy of `p` will be done. 
*/ 
template <typename Polyhedron> 
Polyhedral_mesh_domain_with_features_3(const Polyhedron& p); 

/*! 

Constructs a `Polyhedral_mesh_domain_with_features_3` from an off file. No feature 
detection is done at this level. 
*/ 
Polyhedral_mesh_domain_with_features_3(const std::string& filename); 

/// @} 

/// \name Operations 
/// @{

/*! 
Detects sharp features of the internal polyhedron 
and inserts them in as features of the domain. `angle_bound` gives the maximum dihedral angle 
(in degrees) between two triangles of the internal polyhedron which is used to consider that the edge between 
those triangles is a feature edge. 
*/ 
void detect_features(FT angle_bound=120); 

/// @}

}; /* end Polyhedral_mesh_domain_with_features_3 */
} /* end namespace CGAL */
