/// \addtogroup PkgAlgebraicFoundations Algebraic Foundations
/// @{
/// \addtogroup concepts Concepts
/// @{

 
///  
///  Two types `A` and `B` are a model of the concept 
///  `ImplicitInteroperable`, if there is a superior type, such that 
///  binary arithmetic operations involving `A` and `B` result in 
///  this type. This type is `Coercion_traits::Type`.
///  The type `Coercion_traits::Type` is required to be 
///  implicit constructible from `A` and `B`.
///  
///  
///  In this case `Coercion_traits::Are_implicit_interoperable` 
///  is `Tag_true`.
///  
///  
///  \refines `ExplicitInteroperable`
///  \sa `CGAL::Coercion_traits<A,B>`
///  \sa `ExplicitInteroperable`
///  \sa `AlgebraicStructureTraits`
///  \sa `RealEmbeddableTraits`
class ImplicitInteroperable {
public:

}; /* concept ImplicitInteroperable */
/// @}
/// @} 

                   
  

