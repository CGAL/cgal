namespace CGAL {

/*!
\ingroup kernel_classes3

An object of the class `Weighted_point_3` is a tuple of a 3D point and a scalar weight. 


\cgalHeading{Operators}

The following operations can be applied on points: 

\sa `Point_3<Kernel>`
\sa `Kernel::WeightedPoint_3` 

*/
template< typename Kernel >
class Weighted_point_3 {
public:


/// \name Creation 
/// @{
  
/*!
introduces a weighted point with %Cartesian coordinates `(0,0,0)` and  weight `0`. 
*/ 
Weighted_point_3(const Origin &ORIGIN); 

/*!
introduces a weighted point from point `p` and  weight `0`. 
*/ 
  Weighted_point_3(const Point_3<Kernel>& p); 

/*!
introduces a weighted point from point `p` and  weight `w`. 
*/ 
  Weighted_point_3(const Point_3<Kernel>& p, Kernel::FT& w); 

/*!
introduces a weighted point with coordinates `x`, `y`, `z`, and  weight 0. 
*/ 
  Weighted_point_3(const Kernel::FT& x, const Kernel::FT& y, const Kernel::FT& z); 


/// @} 

/// \name Operations 
/// @{

/*!
returns the point of the weighted point.
*/ 
  Point_3<Kernel> point() const;


/*!
returns the weight of the weighted point.
*/ 
  Kernel::FT weight() const;
/// @}

};

} /* end namespace CGAL */
