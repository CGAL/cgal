
/*!
\ingroup PkgAlgebraicFoundationsAlgebraicStructuresConcepts
\cgalConcept

Two types `A` and `B` are a model of the `ExplicitInteroperable`
concept, if it is possible to derive a superior type for `A` and `B`,
such that both types are embeddable into this type.
This type is \link CGAL::Coercion_traits::Type `CGAL::Coercion_traits<A,B>::Type`\endlink.

In this case
\link CGAL::Coercion_traits::Are_explicit_interoperable  `CGAL::Coercion_traits<A,B>::Are_explicit_interoperable`\endlink
is `Tag_true`.

`A` and `B` are valid argument types for all binary functors in
`CGAL::Algebraic_structure_traits<Type>` and `CGAL::Real_embeddable_traits<Type>`.
This is also the case for the respective global functions.

\sa `CGAL::Coercion_traits<A,B>`
\sa `ImplicitInteroperable`
\sa `AlgebraicStructureTraits`
\sa `RealEmbeddableTraits`

*/

class ExplicitInteroperable {
public:

}; /* end ExplicitInteroperable */
