namespace CGAL {

/*!
\ingroup PkgRidges_3
The enum `Ridge_order` describes the order of differential quantities
used to distinguish elliptic and hyperbolic ridges. Third or fourth
order quantities may be used as explained in section  \ref ridgemesh
of the user manual.
\sa `Ridge_approximation`
*/
enum Ridge_order { Ridge_order_3 = 3, Ridge_order_4};

/*!
\ingroup PkgRidges_3
The enum `Ridge_type` describes the types for the class `Ridge_line`.
\sa `Ridge_line`
*/
enum Ridge_type { 
  MAX_ELLIPTIC_RIDGE = 1, MAX_HYPERBOLIC_RIDGE,
  MAX_CREST_RIDGE, MIN_ELLIPTIC_RIDGE,
  MIN_HYPERBOLIC_RIDGE, MIN_CREST_RIDGE 
};


/*!
\ingroup PkgRidges_3

The function `compute_crest_ridges` is a shortcut to the method of the same name 
of the class `Ridge_approximation`. 

The operator `<<` is overloaded for this class and returns the
line type, strength, sharpness and coordinates of the points of the
polyline.

\sa `Ridge_approximation` 

*/
template < class TriangulatedSurfaceMesh, 
class Vertex2FTPropertyMap,
class Vertex2VectorPropertyMap,
class OutputIterator>
OutputIterator compute_crest_ridges(const TriangulatedSurfaceMesh &P,
const Vertex2FTPropertyMap& vertex2k1_pm, const
Vertex2FTPropertyMap& vertex2k2_pm, const
Vertex2FTPropertyMap& vertex2b0_pm, const
Vertex2FTPropertyMap& vertex2b3_pm, const
Vertex2VectorPropertyMap& vertex2d1_pm, const
Vertex2VectorPropertyMap& vertex2d2_pm, const
Vertex2FTPropertyMap& vertex2P1_pm, const
Vertex2FTPropertyMap& vertex2P2_pm,
OutputIterator it, 
CGAL::Ridge_order order = CGAL::Ridge_order_3);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgRidges_3

The function `compute_max_ridges` is a shortcut to the method of the same name 
of the class `Ridge_approximation`. 

\sa `Ridge_approximation` 

*/
template < class TriangulatedSurfaceMesh, 
class Vertex2FTPropertyMap,
class Vertex2VectorPropertyMap,
class OutputIterator>
OutputIterator compute_max_ridges(const TriangulatedSurfaceMesh &P,
const Vertex2FTPropertyMap& vertex2k1_pm, const
Vertex2FTPropertyMap& vertex2k2_pm, const
Vertex2FTPropertyMap& vertex2b0_pm, const
Vertex2FTPropertyMap& vertex2b3_pm, const
Vertex2VectorPropertyMap& vertex2d1_pm, const
Vertex2VectorPropertyMap& vertex2d2_pm, const
Vertex2FTPropertyMap& vertex2P1_pm, const
Vertex2FTPropertyMap& vertex2P2_pm,
OutputIterator it, 
CGAL::Ridge_order order = CGAL::Ridge_order_3);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgRidges_3

The function `compute_min_ridges` is a shortcut to the method of the same name 
of the class `Ridge_approximation`. 

\sa `Ridge_approximation` 

*/
template < class TriangulatedSurfaceMesh, 
class Vertex2FTPropertyMap,
class Vertex2VectorPropertyMap,
class OutputIterator>
OutputIterator compute_min_ridges(const TriangulatedSurfaceMesh &P,
const Vertex2FTPropertyMap& vertex2k1_pm, const
Vertex2FTPropertyMap& vertex2k2_pm, const
Vertex2FTPropertyMap& vertex2b0_pm, const
Vertex2FTPropertyMap& vertex2b3_pm, const
Vertex2VectorPropertyMap& vertex2d1_pm, const
Vertex2VectorPropertyMap& vertex2d2_pm, const
Vertex2FTPropertyMap& vertex2P1_pm, const
Vertex2FTPropertyMap& vertex2P2_pm,
OutputIterator it, 
CGAL::Ridge_order order = CGAL::Ridge_order_3);

} /* namespace CGAL */


namespace CGAL {

/*!
\ingroup PkgRidges_3

The class `Ridge_approximation` computes the approximation of 
ridges of a triangular polyhedral surface. 

\tparam TriangulatedSurfaceMesh is the surface type. 
\tparam Vertex2FTPropertyMap, Vertex2VectorPropertyMap provide the differential properties of 
the surface associated to its vertices. 

Requirements (checked at compile time): the types 
`TriangulatedSurfaceMesh::Traits::FT` and 
`Vertex2FTPropertyMap::value_type` must coincide; the types 
`TriangulatedSurfaceMesh::Traits::Vector_3` and 
`Vertex2VectorPropertyMap::value_type` must coincide; the types 
`TriangulatedSurfaceMesh::Vertex_handle`, 
`Vertex2FTPropertyMap::key_type` and 
`Vertex2VectorPropertyMap::key_type` must coincide; 

\sa `Ridge_line` 

*/
template< typename TriangulatedSurfaceMesh, typename Vertex2FTPropertyMap, typename Vertex2VectorPropertyMap >
class Ridge_approximation {
public:

/// \name Types 
/// @{

/*! 
Order of differential 
quantities used to distinguish elliptic and hyperbolic ridges. Third 
(`Tag_3`) or fourth (`Tag_4`) order quantities may be used as 
explained in section \ref ridgemesh of the user manual. 
*/ 
enum Tag_order {Tag_3, Tag_4}; 

/// @} 

/// \name Creation 
/// @{

/*! 
The two last property maps may 
not be used if computations are performed with 
the parameter `Tag_3`, in which case these 
property maps shall be initialized with their 
default constructors. 
*/ 
Ridge_approximation(const TriangulatedSurfaceMesh &P, 
const Vertex2FTPropertyMap& vertex2k1_pm, const 
Vertex2FTPropertyMap& vertex2k2_pm, const 
Vertex2FTPropertyMap& vertex2b0_pm, const 
Vertex2FTPropertyMap& vertex2b3_pm, const 
Vertex2VectorPropertyMap& vertex2d1_pm, const 
Vertex2VectorPropertyMap& vertex2d2_pm, const 
Vertex2FTPropertyMap& vertex2P1_pm, const 
Vertex2FTPropertyMap& vertex2P2_pm); 

/*! 
Outputs ridges of types `MAX_ELLIPTIC_RIDGE` and `MAX_HYPERBOLIC_RIDGE`. 
Parameter `it` is an output iterator whose 
value type is `Ridge_line*`. 
*/ 
template <class OutputIterator> OutputIterator compute_max_ridges(OutputIterator it, Tag_order ord = Tag_3); 

/*! 
Outputs ridges of types `MIN_ELLIPTIC_RIDGE` and `MIN_HYPERBOLIC_RIDGE`. 
Parameter `it` is an output iterator whose 
value type is `Ridge_line*`. 
*/ 
template <class OutputIterator> OutputIterator compute_min_ridges(OutputIterator it, Tag_order ord = Tag_3); 

/*! 
Outputs ridges of types `MAX_CREST_RIDGE` and `MIN_CREST_RIDGE`. 
Parameter `it` is an output iterator whose 
value type is `Ridge_line*`. 
*/ 
template <class OutputIterator> OutputIterator compute_crest_ridges(OutputIterator it, Tag_order ord = Tag_3); 

/// @}

}; /* end Ridge_approximation */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgRidges_3

The class `Ridge_line` stores the description of a ridge line. The list 
of halfedges defines a connected sequence of edges (but not as 
oriented halfedges). The scalar \f$ b\f$ paired with a halfedge \f$ pq\f$ is the 
barycentric coordinate of the crossing point \f$ r\f$ with the ridge: 
\f$ r = b\times p + (1-b)\times q\f$. 

\sa `Ridge_approximation`

*/
template< typename TriangulatedSurfaceMesh >
class Ridge_line {
public:

/// \name Types 
/// @{

/*! 

*/ 
typedef typename TriangulatedSurfaceMesh::Traits::FT FT; 

/*! 

*/ 
typedef typename TriangulatedSurfaceMesh::Halfedge_handle Halfedge_handle; 

/*! 
A halfedge crossed by a ridge is paired with the barycentric 
coordinate of the crossing point. 
*/ 
typedef std::pair< Halfedge_handle, FT> Ridge_halfhedge; 

/// @} 

/// \name Creation 
/// @{

/*! 
default constructor. 
*/ 
Ridge_line(); 

/// @} 

/// \name Access Functions 
/// @{

/*! 

*/ 
const Ridge_type line_type(); 

/*! 

*/ 
const FT strength(); 

/*! 

*/ 
const FT sharpness(); 

/*! 

*/ 
const std::list<Ridge_halfhedge>* line(); 

/*!
Writes the line type, strength, sharpness and coordinates of the
points of the polyline to `o`.
*/
template< typename TriangulatedSurfaceMesh >
std::ostream& operator<<(std::ostream& o, const Ridge_line<TriangulatedSurfaceMesh>&);

/// @}

}; /* end Ridge_line */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgRidges_3

The class `Vertex2Data_Property_Map_with_std_map` is a model of the concepts 
`Vertex2FTPropertyMap` and `Vertex2VectorPropertyMap` to be used for 
`Ridge_approximation`. The property maps are 
created with the `boost::associative_property_map` adaptor from 
`std::map`. 

\cgalModels `Vertex2FTPropertyMap`
\cgalModels `Vertex2VectorPropertyMap`

\sa `Ridge_approximation` 

*/
template< typename TriangulatedSurfaceMesh >
class Vertex2Data_Property_Map_with_std_map {
public:

/// \name Types 
/// @{

/*! 

*/ 
typedef typename TriangulatedSurfaceMesh::Traits::FT FT; 

/*! 

*/ 
typedef typename TriangulatedSurfaceMesh::Traits::Vector_3 Vector_3; 

/*! 

*/ 
typedef typename TriangulatedSurfaceMesh::Vertex_handle Vertex_handle; 

/*!
  \ingroup PkgRidges_3
*/
struct Vertex_cmp{bool operator();}; 

/*! 

*/ 
typedef std::map<Vertex_handle, FT, Vertex_cmp> Vertex2FT_map; 

/*! 

*/ 
typedef boost::associative_property_map< Vertex2FT_map > Vertex2FT_property_map; 

/*! 

*/ 
typedef std::map<Vertex_handle, Vector_3, Vertex_cmp> Vertex2Vector_map; 

/*! 

*/ 
typedef boost::associative_property_map< Vertex2Vector_map > Vertex2Vector_property_map; 

/// @}

}; /* end Vertex2Data_Property_Map_with_std_map */

} /* end namespace CGAL */
