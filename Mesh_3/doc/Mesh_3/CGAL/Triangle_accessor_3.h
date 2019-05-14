namespace CGAL {

/*!
\ingroup PkgMesh_3Domains

The class `Triangle_accessor_3` is a model for the concept `TriangleAccessor_3`. It is 
designed to serve as accessor for objects of type `Polyhedron_3<K>`. 

\attention Actually, the class `Triangle_accessor_3` is a partial specialization of the class 
template `template<typename Polyhedron, typename K> 
Triangle_accessor_3<Polyhedron, K>`. One may give another partial 
specialization of this class to handle one's own polyhedron data structure. 


\tparam K is the geometric traits class. 

\cgalModels `TriangleAccessor_3`

\sa `CGAL::Polyhedral_mesh_domain_3<Polyhedron,IGT,TriangleAccessor>` 

*/
template< CGAL::Polyhedron<typename K>, typename K >
class Triangle_accessor_3 {
public:

/// \name Types 
/// @{

/*!
Triangle iterator. 
*/ 
typedef Polyhedron_3<K>::Facet_const_iterator 
Triangle_iterator; 

/*!
Triangle 
handle. 
*/ 
typedef Polyhedron_3<K>::Facet_const_handle Triangle_handle; 

/*!
Triangle type. 
*/ 
typedef K::Triangle_3 Triangle_3; 

/// @}

}; /* end Triangle_accessor_3 */
} /* end namespace CGAL */
