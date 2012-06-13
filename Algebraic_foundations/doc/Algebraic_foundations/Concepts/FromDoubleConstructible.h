/// \addtogroup PkgAlgebraicFoundations Algebraic Foundations
/// @{
/// \addtogroup PkgAlgebraicFoundationsConcepts Concepts
/// @{

 
///  
///  
///  
///  A model of the concept `FromDoubleConstructible` is required 
///  to be constructible from the type `double`.
///  In case the type is a model of `RealEmbeddable` too, for any double d
///  the identity: ` d == CGAL::to_double(T(d))`, is guaranteed.
///    
///  
///  
///  
class FromDoubleConstructible {
public:

/// \name Creation
/// @{
/*!
 conversion constructor from double. 
*/
FromDoubleConstructible(const double& d);
/// @}

 
}; /* concept FromDoubleConstructible */
/// @}
/// @} 

                   
  

