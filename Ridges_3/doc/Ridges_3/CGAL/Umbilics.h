namespace CGAL {

/*!
\ingroup PkgRidges_3Enums
The enum `Umbilic_type` describes the types for the class `Umbilic`.
\sa `Umbilic`
*/
enum Umbilic_type { NON_GENERIC_UMBILIC, ELLIPTIC_UMBILIC, HYPERBOLIC_UMBILIC };


/*!
\ingroup PkgRidges_3

The function `compute_umbilics` is a shortcut to the method `compute` of 
the class `Umbilic_approximation`. 

\sa `Umbilic_approximation` 

*/
template < class TriangulatedSurfaceMesh, 
class Vertex2FTPropertyMap,
class Vertex2VectorPropertyMap,
class OutputIterator>
OutputIterator compute_umbilics(const TriangulatedSurfaceMesh &P,
const Vertex2FTPropertyMap& vertex2k1_pm, 
const Vertex2FTPropertyMap& vertex2k2_pm,
const Vertex2VectorPropertyMap& vertex2d1_pm, 
const Vertex2VectorPropertyMap& vertex2d2_pm,
OutputIterator it, 
double size);

} /* namespace CGAL */


namespace CGAL {

/*!
\ingroup PkgRidges_3

The class `Umbilic_approximation` computes the approximation of 
umbilics on a triangular polyhedral surface. 

\tparam TriangulatedSurfaceMesh is the surface type. 
\tparam Vertex2FTPropertyMap, Vertex2VectorPropertyMap provide 
the differential properties of the surface associated to its vertices. 

Requirements (checked at compile time) : the types 
`TriangulatedSurfaceMesh::Traits::FT` and 
`Vertex2FTPropertyMap::value_type` must coincide; the types 
`TriangulatedSurfaceMesh::Traits::Vector_3` and 
`Vertex2VectorPropertyMap::value_type` must coincide; the types 
`TriangulatedSurfaceMesh::Vertex_handle`, 
`Vertex2FTPropertyMap::key_type` and 
`Vertex2VectorPropertyMap::key_type` must coincide; 

\sa `Umbilic` 
\sa `TriangulatedSurfaceMesh` 
\sa `Vertex2FTPropertyMap` 
\sa `Vertex2VectorPropertyMap` 

*/
template< typename TriangulatedSurfaceMesh, typename Vertex2FTPropertyMap, typename Vertex2VectorPropertyMap >
class Umbilic_approximation {
public:

/// \name Types 
/// @{

/*! 

*/ 
typedef typename TriangulatedSurfaceMesh::Traits::FT FT; 

/// @} 

/// \name Creation 
/// @{

/*! 
default constructor. 
*/ 
Umbilic_approximation(const TriangulatedSurfaceMesh& P, 
const Vertex2FTPropertyMap& vertex2k1_pm, 
const Vertex2FTPropertyMap& vertex2k2_pm, 
const Vertex2VectorPropertyMap& vertex2d1_pm, 
const Vertex2VectorPropertyMap& vertex2d2_pm); 

/// @} 

/// \name Operations 
/// @{

/*! 
Performs the approximation, `size` determines the size of the 
patches around vertices, taken as `size` times the size of the 
1-ring. Umbilics are inserted into the `OutputIterator` `it` with value type `Umbilic*`. 
*/ 
template <class OutputIterator> OutputIterator compute(OutputIterator it, FT size); 

/// @}

}; /* end Umbilic_approximation */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgRidges_3

The class `Umbilic` stores the description of an umbilic. 

### Operations ###

The insert operator `<<` is overloaded for `Umbilic`, it 
gives the location (3d coordinates of the vertex) and the type. 

\sa `Umbilic_approximation` 

*/
template< typename TriangulatedSurfaceMesh >
class Umbilic {
public:

/// \name Types 
/// @{

/*! 

*/ 
typedef typename TriangulatedSurfaceMesh::Vertex_handle Vertex_handle; 

/*! 

*/ 
typedef typename TriangulatedSurfaceMesh::Halfedge_handle Halfedge_handle; 

/// @} 

/// \name Access Functions 
/// @{

/*! 

*/ 
const Vertex_handle vertex(); 

/*! 

*/ 
const Umbilic_type umbilic_type(); 

/*! 

*/ 
const std::list<Halfedge_handle>& contour_list() ; 

/// @}

}; /* end Umbilic */
} /* end namespace CGAL */
