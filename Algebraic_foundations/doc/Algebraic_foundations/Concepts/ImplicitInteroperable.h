
/*!
\ingroup PkgAlgebraicFoundationsInteroperabilityConcepts
\cgalconcept

Two types `A` and `B` are a model of the concept 
`ImplicitInteroperable`, if there is a superior type, such that 
binary arithmetic operations involving `A` and `B` result in 
this type. This type is \link CGAL::Coercion_traits::Type `CGAL::Coercion_traits<A,B>::Type`\endlink.

The type \link CGAL::Coercion_traits::Type `CGAL::Coercion_traits<A,B>::Type`\endlink is required to be 
implicit constructible from `A` and `B`. 

In this case 
\link CGAL::Coercion_traits::Are_implicit_interoperable  `CGAL::Coercion_traits<A,B>::Are_implicit_interoperable`\endlink
is `CGAL::Tag_true`. 

\refines ::ExplicitInteroperable 

\sa `CGAL::Coercion_traits<A,B>` 
\sa `ExplicitInteroperable` 
\sa `AlgebraicStructureTraits` 
\sa `RealEmbeddableTraits` 

*/

class ImplicitInteroperable {
public:

/// @}

}; /* end ImplicitInteroperable */

