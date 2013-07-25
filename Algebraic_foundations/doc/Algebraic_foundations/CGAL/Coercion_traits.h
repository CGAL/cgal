
namespace CGAL {

/*!
\ingroup PkgAlgebraicFoundations

An instance of `Coercion_traits` reflects the type coercion of the types 
<span class="textsc">A</span> and <span class="textsc">B</span>, it is symmetric in the two template arguments. 

\sa `ExplicitInteroperable` 
\sa `ImplicitInteroperable` 

*/
template< typename A, typename B >
class Coercion_traits {
public:

/// \name Types 
/// @{

/*!
Tag indicating whether the two types A and B are a model of `ExplicitInteroperable` 

This is either `CGAL::Tag_true` or `CGAL::Tag_false`. 
*/ 
typedef unspecified_type Are_explicit_interoperable; 

/*!
Tag indicating whether the two types A and B are a model of `ImplicitInteroperable` 

This is either `CGAL::Tag_true` or `CGAL::Tag_false`. 
*/ 
typedef unspecified_type Are_implicit_interoperable; 

/*!
The coercion type of `A` and `B`. 

In case A and B are not `ExplicitInteroperable` this is undefined. 
*/ 
typedef unspecified_type Type; 

/*!
A model of the `AdaptableFunctor` concept, providing the conversion of `A` or `B` to `Type`. 

In case A and B are not `ExplicitInteroperable` this is undefined. 
*/ 
typedef unspecified_type Cast; 

/// @}

}; /* end Coercion_traits */
} /* end namespace CGAL */
