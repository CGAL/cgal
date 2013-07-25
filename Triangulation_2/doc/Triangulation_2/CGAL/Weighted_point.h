
namespace CGAL {

/*!
\ingroup PkgTriangulation2Miscellaneous

The class `Weighted_point` provides a type associating 
a point type `Pt` with a weight type `Wt`. 
It is used in the traits classes `Regular_triangulation_euclidean_traits_2` 
and `Regular_triangulation_euclidean_traits_3`. 

\sa `CGAL::Regular_triangulation_euclidean_traits_2<Rep,Weight>` 
\sa `CGAL::Regular_triangulation_euclidean_traits_3<R,Weight>`

*/
template< typename Pt, typename Wt >
class Weighted_point : public Pt {
public:

/// \name Types 
/// @{

/*!
The point type 
*/ 
Pt Point; 

/*!
The weight type. 
*/ 
Wt Weight; 

/// @} 

/// \name Creation 
/// @{

/*!
copy constructor. 
*/ 
Weighted_point(Weighted_point wq); 

/*!

*/ 
Weighted_point(Point p=Point(), Weight w= Weight(0)); 

/*!
Constructs the point from `x` 
and `y` coordinates, with a weight of 0. Requires that the ambient 
dimension be 2. 
*/ 
Weighted_point(FT x, FT y); 

/*!
Constructs the point from 
`x`, `y` and `z` coordinates, with a weight of 0. Requires that 
the ambient dimension be 3. 
*/ 
Weighted_point(FT x, FT y, FT z); 

/// @} 

/// \name Access Functions 
/// @{

/*!

*/ 
Point point() const; 

/*!

*/ 
Weight weight() const; 

/// @}

}; /* end Weighted_point */
} /* end namespace CGAL */
