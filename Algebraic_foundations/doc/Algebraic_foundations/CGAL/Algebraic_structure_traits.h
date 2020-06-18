namespace CGAL {

/*!
\ingroup PkgAlgebraicFoundationsRef

An instance of `Algebraic_structure_traits` is a model of `AlgebraicStructureTraits`, where <span class="textsc">T</span> is the associated type.

\cgalModels `AlgebraicStructureTraits`

*/
template< typename T >
class Algebraic_structure_traits {

}; /* end Algebraic_structure_traits */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgAlgebraicFoundationsRef

Tag indicating that a type is a model of the
`EuclideanRing` concept.

\cgalModels `DefaultConstructible`

\sa `EuclideanRing`
\sa `AlgebraicStructureTraits`

*/

struct Euclidean_ring_tag : public Unique_factorization_domain_tag {

}; /* end Euclidean_ring_tag */

/*!
\ingroup PkgAlgebraicFoundationsRef

Tag indicating that a type is a model of the `Field` concept.

\cgalModels `DefaultConstructible`

\sa `Field`
\sa `AlgebraicStructureTraits`

*/

struct Field_tag : public Integral_domain_tag {

}; /* end Field_tag */

/*!
\ingroup PkgAlgebraicFoundationsRef

Tag indicating that a type is a model of the `FieldWithKthRoot` concept.

\cgalModels `DefaultConstructible`

\sa `FieldWithKthRoot`
\sa `AlgebraicStructureTraits`

*/

struct Field_with_kth_root_tag : public Field_with_sqrt_tag {

}; /* end Field_with_kth_root_tag */

/*!
\ingroup PkgAlgebraicFoundationsRef

Tag indicating that a type is a model of the `FieldWithRootOf` concept.

\cgalModels `DefaultConstructible`

\sa `FieldWithRootOf`
\sa `AlgebraicStructureTraits`

*/

struct Field_with_root_of_tag : public Field_with_kth_root_tag {

}; /* end Field_with_root_of_tag */

/*!
\ingroup PkgAlgebraicFoundationsRef

Tag indicating that a type is a model of the `FieldWithSqrt` concept.

\cgalModels `DefaultConstructible`

\sa `FieldWithSqrt`
\sa `AlgebraicStructureTraits`

*/

struct Field_with_sqrt_tag : public Field_tag {

}; /* end Field_with_sqrt_tag */

/*!
\ingroup PkgAlgebraicFoundationsRef

Tag indicating that a type is a model of the `IntegralDomain` concept.

\cgalModels `DefaultConstructible`

\sa `IntegralDomain`
\sa `AlgebraicStructureTraits`

*/

struct Integral_domain_tag : public Integral_domain_without_division_tag {

}; /* end Integral_domain_tag */

/*!
\ingroup PkgAlgebraicFoundationsRef

Tag indicating that a type is a model of the `IntegralDomainWithoutDivision` concept.

\cgalModels `DefaultConstructible`

\sa `IntegralDomainWithoutDivision`

*/

struct Integral_domain_without_division_tag {

}; /* end Integral_domain_without_division_tag */

/*!
\ingroup PkgAlgebraicFoundationsRef

Tag indicating that a type is a model of the `UniqueFactorizationDomain` concept.

\cgalModels `DefaultConstructible`

\sa `UniqueFactorizationDomain`
\sa `AlgebraicStructureTraits`

*/

struct Unique_factorization_domain_tag : public Integral_domain_tag {

}; /* end Unique_factorization_domain_tag */
} /* end namespace CGAL */
