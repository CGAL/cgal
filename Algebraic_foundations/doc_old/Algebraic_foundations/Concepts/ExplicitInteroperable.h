/// \addtogroup PkgAlgebraicFoundations Algebraic Foundations
/// @{
/// \addtogroup PkgAlgebraicFoundationsConcepts Concepts
/// @{

 
///  
///  Two types `A` and `B` are a model of the `ExplicitInteroperable` 
///  concept, if it is possible to derive a superior type for `A` and `B`,
///  such that both types are embeddable into this type. 
///  This type is `Coercion_traits::Type`.
///  In this case `Coercion_traits::Are_explicit_interoperable` 
///  is `Tag_true`.
///  `A` and `B` are valid argument types for all binary functors in 
///  `Algebraic_structure_traits<Type>` and `Real_embeddable_traits<Type>`.   
///  This is also the case for the respective global functions.
///  \sa `CGAL::Coercion_traits<A,B>`
///  \sa `ImplicitInteroperable`
///  \sa `AlgebraicStructureTraits`
///  \sa `RealEmbeddableTraits`
class ExplicitInteroperable {
public:

}; /* concept ExplicitInteroperable */
/// @}
/// @} 

                   
  

