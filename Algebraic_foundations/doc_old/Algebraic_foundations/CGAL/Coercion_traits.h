namespace CGAL {
/// \addtogroup PkgAlgebraicFoundations Algebraic Foundations
/// @{

 
///  
///  An instance of `Coercion_traits` reflects the type coercion of
///  the types <span class="textsc">A</span> and <span
///  class="textsc">B</span>, it is symmetric in the two template
///  arguments.
///  
///  <span class="textsc">A</span> and <span class="textsc">B</span>,
///  it is symmetric in the two template arguments.
template< class A, class B >
class Coercion_traits {
public:

/// \name Types
/// @{
/*!
  Tag indicating whether the two types A and B are a model of `ExplicitInteroperable` 
          This is either `CGAL::Tag_true` or `CGAL::Tag_false`. 
*/
typedef Hidden_type Are_explicit_interoperable;
/// @}

/// \name Types
/// @{
/*!
  Tag indicating whether the two types A and B are a model of `ImplicitInteroperable` 
          This is either `CGAL::Tag_true` or `CGAL::Tag_false`. 
*/
typedef Hidden_type Are_implicit_interoperable;
/// @}

/// \name Types
/// @{
/*!
 The coercion type of ` A` and ` B`. 
        In case A and B are not `ExplicitInteroperable` this is undefined.  
*/
typedef Hidden_type Type;
/// @}

/// \name Types
/// @{
/*!
 A model of the `AdaptableFunctor` concept, providing the conversion of ` A` or ` B` to ` Type`. 
        In case A and B are not `ExplicitInteroperable` this is undefined.  
*/
typedef Hidden_type Cast;
/// @}

}; /* class Coercion_traits */
/// @}
} // namespace CGAL

                   
  

