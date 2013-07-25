namespace CGAL {

/*!
\ingroup PkgMesh_3Domains

The class `Mesh_domain_with_polyline_features_3` is designed to allow the user 
to add some 0- and 1-dimensional 
features into any model of the `MeshDomain_3` concept. 
The 1-dimensional features are described as polylines 
whose endpoints are the added corners. 

\tparam MeshDomain_3 is the type 
of the domain which should be extended. 
It has to be a model of the `MeshDomain_3` concept. 

\cgalModels `MeshDomainWithFeatures_3` 

\sa `MeshDomain_3` 
\sa `MeshPolyline_3` 
\sa `CGAL::Implicit_mesh_domain_3<Function,BGT>` 
\sa `CGAL::Polyhedral_mesh_domain_3<Polyhedron,IGT,TriangleAccessor>` 
\sa `CGAL::Labeled_image_mesh_domain_3<Image,BGT>` 

*/
template< typename MeshDomain_3 >
class Mesh_domain_with_polyline_features_3 : public MeshDomain_3 {
public:

/// \name Types 
/// @{

/*!
`Corner_index` type. 
*/ 
typedef int Corner_index; 

/*!
`Curve_segment_index` type. 
*/ 
typedef int Curve_segment_index; 

/// @} 

/// \name Creation 
/// @{

/*!
Constructor. Forwards the arguments to the constructor 
of the base class. 
*/ 
template <typename ...T> 
Mesh_domain_with_polyline_features_3(T ...t); 

/// @} 

/// \name Operations 
/// @{

/*!

Add 1-dimensional features in the domain. `InputIterator` value type must 
be a model of the concept `MeshPolyline_3`. 
*/ 
template <typename InputIterator> 
void add_features(InputIterator begin, InputIterator beyond); 

/*!

Add 1-dimensional features in the domain with their incidences with 2-dimensional 
features of the domain. The `InputIterator` value type must be 
`std::pair<Polyline, std::pair<InputSurfacePatchIndexIterator, InputSurfacePatchIndexIterator> >` 
where `Polyline` must be a model of the concept `MeshPolyline_3` 
and the internal pair gives a range on surface patches indices which are incident 
to the polyline. 
*/ 
template <typename InputIterator> 
void add_features_and_incidences(InputIterator begin, InputIterator beyond); 

/// @}

}; /* end Mesh_domain_with_polyline_features_3 */
} /* end namespace CGAL */
