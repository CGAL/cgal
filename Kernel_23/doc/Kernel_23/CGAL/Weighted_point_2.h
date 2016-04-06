namespace CGAL {

/*!
\ingroup kernel_classes2

An object of the class `Weighted_point_2` is a tuple of a 2D point and a scalar weight. 


\cgalHeading{Operators}

The following operations can be applied on points: 

\sa `Point_2<Kernel>`
\sa `Kernel::WeightedPoint_2` 

*/
template< typename Kernel >
class Weighted_point_2 {
public:


/// \name Creation 
/// @{
  
/*!
introduces a weighted point with %Cartesian coordinates `(0,0,0)` and  weight `0`. 
*/ 
Weighted_point_2(const Origin &ORIGIN); 

/*!
introduces a weighted point from point `p` and  weight `0`. 
*/ 
  Weighted_point_2(const Point_2<Kernel>& p); 

/*!
introduces a weighted point from point `p` and  weight `w`. 
*/ 
  Weighted_point_2(const Point_2<Kernel>& p, Kernel::FT& w); 


/// @} 

/// \name Operations 
/// @{

/*!
returns the weight of the weighted point.
*/ 
  Point_2<Kernel> point() const;


/*!
returns the weight of the weighted point.
*/ 
  Kernel::FT weight() const;
/// @}

};

} /* end namespace CGAL */
