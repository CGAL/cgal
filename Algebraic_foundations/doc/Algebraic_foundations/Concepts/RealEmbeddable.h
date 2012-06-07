/// \addtogroup PkgAlgebraicFoundations Algebraic Foundations
/// @{
/// \addtogroup concepts Concepts
/// @{

 
///  
///  A model of this concepts represents numbers that are embeddable on the real 
///  axis. The type obeys the algebraic structure and compares two values according 
///  to the total order of the real numbers.
///  Moreover, `CGAL::Real_embeddable_traits< RealEmbeddable >` is a model of 
///  `RealEmbeddableTraits`
///  with: 
///  - `CGAL::Real_embeddable_traits::Is_real_embeddable` set to `Tag_true` 
///  and functors :
///  - `CGAL::Real_embeddable_traits::Is_zero`  
///  - `CGAL::Real_embeddable_traits::Abs`  
///  - `CGAL::Real_embeddable_traits::Sgn`  
///  - `CGAL::Real_embeddable_traits::Is_positive`  
///  - `CGAL::Real_embeddable_traits::Is_negative`  
///  - `CGAL::Real_embeddable_traits::Compare`  
///  - `CGAL::Real_embeddable_traits::To_double`  
///  - `CGAL::Real_embeddable_traits::To_interval`  
///  Remark:
///  If a number type is a model of both `IntegralDomainWithoutDivision` and 
///  `RealEmbeddable`, it follows that the ring represented by such a number type 
///  is a sub-ring of the real numbers and hence has characteristic zero.
///  
///  
///  \refines `Equality Comparable`
///  \refines `LessThanComparable`
///  \sa `RealEmbeddableTraits`
class RealEmbeddable {
public:

}; /* concept RealEmbeddable */

/// @}
/// @} 

///
/// \relates RealEmbeddable
/// 
bool operator==(const RealEmbeddable &a, 
                            const RealEmbeddable &b);

///
/// \relates RealEmbeddable
/// 
bool operator!=(const RealEmbeddable &a, 
                            const RealEmbeddable &b);

///
/// \relates RealEmbeddable
/// 
bool operator< (const RealEmbeddable &a, 
                            const RealEmbeddable &b);

///
/// \relates RealEmbeddable
/// 
bool operator<=(const RealEmbeddable &a, 
                            const RealEmbeddable &b);

///
/// \relates RealEmbeddable
/// 
bool operator> (const RealEmbeddable &a, 
                            const RealEmbeddable &b);

///
/// \relates RealEmbeddable
/// 
bool operator>=(const RealEmbeddable &a, 
                            const RealEmbeddable &b);

                   
  

