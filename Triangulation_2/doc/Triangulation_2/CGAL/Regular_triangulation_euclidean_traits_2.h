
namespace CGAL {

/*!
\ingroup PkgTriangulation2TraitsClasses

\deprecated The class is deprecated since \cgal 4.10, as the weighted point and the function
objects for weighted points are part of the concept `Kernel`. The class is kept for  backward 
compatibility, but ignores the template parameter `Weight`. 

`Regular_triangulation_euclidean_traits_2` is a model for the concept `RegularTriangulationTraits_2` 

\tparam K must be a model of the `Kernel` concept.
\tparam Weight  This template parameter is ignored, as `Kernel::Weighted_point_2` uses the type `Kernel::FT`.

 
This class inherits from `K` 
and uses a `K::Weighted_point_2` as point type.


\cgalModels `RegularTriangulationTraits_2`

\sa `CGAL::Regular_triangulation_2` 

*/
template< typename K, typename Weight >
class Regular_triangulation_euclidean_traits_2 : public K {
public:
/// \name Types 
/// @{

/*!
The type for point \f$ p\f$ of a weighted point \f$ {p}^{(w)}=(p,w_p)\f$. 
*/ 
typedef K::Point_2 Bare_point; 

/*!
The type for weighted points. 
*/ 
  typedef K::Weighted_point_2 Weighted_point_2; 

/*!
The type for points. 
*/ 
  typedef K::Weighted_point_2 Point_2; 


/// @}

}; /* end Regular_triangulation_euclidean_traits_2 */
} /* end namespace CGAL */
