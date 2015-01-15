namespace CGAL {

/*!
\ingroup RangeQueryItemClasses

The class `Fuzzy_sphere` implements fuzzy `d`-dimensional spheres. 
A fuzzy sphere with radius \f$ r\f$ and fuzziness value \f$ \epsilon\f$ has 
as outer approximation a sphere with radius \f$ r+\epsilon\f$ and 
as inner approximation a sphere with radius \f$ r-\epsilon\f$. 

\cgalHeading{Parameters}

\tparam Traits must be a model of the concept 
`SearchTraits`, for example `CGAL::Cartesian_d<double>`. 

\cgalModels `FuzzyQueryItem`

\sa `FuzzyQueryItem` 

*/
template< typename Traits >
class Fuzzy_sphere {
public:

/// \name Types 
/// @{

/*!
Dimension Tag.
*/
typedef unspecified_type D;

/*!
Point type. 
*/ 
typedef Traits::Point_d Point_d; 

/*!
Number type. 
*/ 
typedef Traits::FT FT; 

/// @} 

/// \name Creation 
/// 
/// @{

/*!
Constructs a fuzzy sphere 
centered at `center` with radius `radius` and fuzziness value `epsilon`. 
*/ 
Fuzzy_sphere(Point_d center, FT radius, FT epsilon=FT(0),Traits t=Traits()); 

/*!
Constructs a fuzzy sphere centered at `center` with radius `radius` and fuzziness value `epsilon`. 
\attention Only available in case `Traits` is `Search_traits_adapter<Key,PointPropertyMap,BaseTraits>`.
*/ 
Fuzzy_sphere(Traits::Base::Point_d center, FT radius, FT epsilon=FT(0), Traits t=Traits()); 

/// @} 

/// \name Operations 
/// @{

/*!
test whether the fuzzy sphere contains `p`. 
*/ 
bool contains(const Point_d& p) const; 

/*!
test whether the inner sphere intersects the rectangle 
associated with a node of a tree. 
*/ 
bool inner_range_intersects(const Kd_tree_rectangle<FT,D>& rectangle) const; 

/*!
test whether the outer sphere encloses the rectangle associated with a node of a tree. 
*/ 
bool outer_range_contains(const Kd_tree_rectangle<FT,D>& rectangle) const; 

/// @}

}; /* end Fuzzy_sphere */
} /* end namespace CGAL */
