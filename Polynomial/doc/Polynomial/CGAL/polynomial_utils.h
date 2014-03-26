namespace CGAL {

/*!
\ingroup PkgPolynomialFunctions

For a given `Polynomial_d`, adapts the 
according functor in `Polynomial_traits_d<Polynomial_d>`. 

Adapts `Polynomial_traits_d::Canonicalize`. 

For more details see the concept `PolynomialTraits_d::Canonicalize`. 



\sa `Polynomial_d` 
\sa `PolynomialTraits_d` 
\sa `PolynomialTraits_d::Canonicalize` 
*/

template <class Polynomial_d>
Polynomial_traits_d<Polynomial_d>::Canonicalize::result_type 
canonicalize(const Polynomial_d& p);


/*!
\ingroup PkgPolynomialFunctions

For a given `Polynomial_d`, adapts the 
according functor in `Polynomial_traits_d<Polynomial_d>`. 

Adapts `Polynomial_traits_d::Compare`. 

For more details see the concept `PolynomialTraits_d::Compare`. 

\sa `Polynomial_d` 
\sa `PolynomialTraits_d` 
\sa `PolynomialTraits_d::Compare` 
*/
template <class Polynomial_d>
Polynomial_traits_d<Polynomial_d>::Compare::result_type 
compare(const Polynomial_d& p, const Polynomial_d& q);


/*!
\ingroup PkgPolynomialFunctions

For a given `Polynomial_d`, adapts the 
according functor in `Polynomial_traits_d<Polynomial_d>`. 

Adapts `Polynomial_traits_d::Degree`. 

For more details see the concept `PolynomialTraits_d::Degree`. 

\sa `Polynomial_d` 
\sa `PolynomialTraits_d` 
\sa `PolynomialTraits_d::Degree` 


*/
template <class Polynomial_d>
Polynomial_traits_d<Polynomial_d>::Degree::result_type 
degree(const Polynomial_d& p, 
int i, 
index = Polynomial_traits_d<Polynomial_d>::d-1);


/*!
\ingroup PkgPolynomialFunctions

For a given `Polynomial_d`, adapts the 
according functor in `Polynomial_traits_d<Polynomial_d>`. 

Adapts `Polynomial_traits_d::DegreeVector`. 

For more details see the concept `PolynomialTraits_d::DegreeVector`. 

\sa `Polynomial_d` 
\sa `PolynomialTraits_d` 
\sa `PolynomialTraits_d::Degree_vector` 


*/

template <class Polynomial_d>
Polynomial_traits_d<Polynomial_d>::Degree_vector::result_type 
degree_vector(const Polynomial_d& p);


/*!
\ingroup PkgPolynomialFunctions

For a given `Polynomial_d`, adapts the 
according functor in `Polynomial_traits_d<Polynomial_d>`. 

Adapts `Polynomial_traits_d::Differentiate`. 

For more details see the concept `PolynomialTraits_d::Differentiate`. 

\sa `Polynomial_d` 
\sa `PolynomialTraits_d` 
\sa `PolynomialTraits_d::Differentiate` 


*/

template <class Polynomial_d>
Polynomial_traits_d<Polynomial_d>::Differentiate::result_type 
differentiate(const Polynomial_d& p,
index = Polynomial_traits_d<Polynomial_d>::d-1 );


/*!
\ingroup PkgPolynomialFunctions

For a given `Polynomial_d`, adapts the 
according functor in `Polynomial_traits_d<Polynomial_d>`. 

Adapts `Polynomial_traits_d::Evaluate_homogeneous`. 

For more details see the concept `PolynomialTraits_d::EvaluateHomogeneous`. 

\sa `Polynomial_d` 
\sa `PolynomialTraits_d` 
\sa `PolynomialTraits_d::EvaluateHomogeneous` 


*/

template <class Polynomial_d>
Polynomial_traits_d<Polynomial_d>::Evaluate_homogeneous::result_type 
evaluate_homogeneous(
const Polynomial_d& p,
Polynomial_traits_d<Polynomial_d>::Coefficient_type u,
Polynomial_traits_d<Polynomial_d>::Coefficient_type v);


/*!
\ingroup PkgPolynomialFunctions

For a given `Polynomial_d`, adapts the 
according functor in `Polynomial_traits_d<Polynomial_d>`. 

Adapts `Polynomial_traits_d::Evaluate`. 

For more details see the concept `PolynomialTraits_d::Evaluate`. 

\sa `Polynomial_d` 
\sa `PolynomialTraits_d` 
\sa `PolynomialTraits_d::Evaluate` 


*/

template <class Polynomial_d>
Polynomial_traits_d<Polynomial_d>::Evaluate::result_type 
evaluate(const Polynomial_d& p,
Polynomial_traits_d<Polynomial_d>::Coefficient_type x);


/*!
\ingroup PkgPolynomialFunctions

For a given `Polynomial_d`, adapts the 
according functor in `Polynomial_traits_d<Polynomial_d>`. 

Adapts `Polynomial_traits_d::Gcd_up_to_constant_factor`. 

For more details see the concept `PolynomialTraits_d::GcdUpToConstantFactor`. 

\sa `Polynomial_d` 
\sa `PolynomialTraits_d` 
\sa `PolynomialTraits_d::GcdUpToConstantFactor` 


*/

template <class Polynomial_d>
Polynomial_traits_d<Polynomial_d>::Gcd_up_to_constant_factor::result_type 
gcd_up_to_constant_factor(const Polynomial_d& p, const Polynomial_d& q);


/*!
\ingroup PkgPolynomialFunctions

For a given `Polynomial_d`, adapts the 
according functor in `Polynomial_traits_d<Polynomial_d>`. 

Adapts `Polynomial_traits_d::GetCoefficient`. 

For more details see the concept `PolynomialTraits_d::GetCoefficient`. 

\sa `Polynomial_d` 
\sa `PolynomialTraits_d` 
\sa `PolynomialTraits_d::GetCoefficient` 
\sa `PolynomialTraits_d::GetInnermostCoefficient` 


*/

template <class Polynomial_d>
Polynomial_traits_d<Polynomial_d>::get_coefficient::result_type 
get_coefficient(const Polynomial_d& p, int i);


/*!
\ingroup PkgPolynomialFunctions

For a given `Polynomial_d`, adapts the 
according functor in `Polynomial_traits_d<Polynomial_d>`. 

Adapts `Polynomial_traits_d::GetInnermostCoefficient`. 

For more details see the concept `PolynomialTraits_d::GetInnermostCoefficient`. 

\sa `Polynomial_d` 
\sa `PolynomialTraits_d` 
\sa `PolynomialTraits_d::GetCoefficient` 
\sa `PolynomialTraits_d::GetInnermostCoefficient` 


*/

template <class Polynomial_d>
Polynomial_traits_d<Polynomial_d>::get_innermost_coefficient::result_type 
get_innermost_coefficient(const Polynomial_d& p, Exponent_vector ev);


/*!
\ingroup PkgPolynomialFunctions

For a given `Polynomial_d`, adapts the 
according functor in `Polynomial_traits_d<Polynomial_d>`. 

Adapts `Polynomial_traits_d::InnermostLeadingCoefficient`. 

For more details see the concept 
`PolynomialTraits_d::InnermostLeadingCoefficient`. 

\sa `Polynomial_d` 
\sa `PolynomialTraits_d` 
\sa `PolynomialTraits_d::InnermostLeadingCoefficient` 


*/

template <class Polynomial_d>
Polynomial_traits_d<Polynomial_d>::Innermost_leading_coefficient::result_type 
innermost_leading_coefficient(const Polynomial_d& p);


/*!
\ingroup PkgPolynomialFunctions

For a given `Polynomial_d`, adapts the 
according functor in `Polynomial_traits_d<Polynomial_d>`. 

Adapts `Polynomial_traits_d::Integral_division_up_to_constant_factor`. 

For more details see the concept `PolynomialTraits_d::IntegralDivisionUpToConstantFactor`. 

\sa `Polynomial_d` 
\sa `PolynomialTraits_d` 
\sa `PolynomialTraits_d::IntegralDivisionUpToConstantFactor` 


*/

template <class Polynomial_d>
Polynomial_traits_d<Polynomial_d>::Integral_division_up_to_constant_factor::result_type 
integral_division_up_to_constant_factor(const Polynomial_d& p, const Polynomial_d& q);


/*!
\ingroup PkgPolynomialFunctions

For a given `Polynomial_d`, adapts the 
according functor in `Polynomial_traits_d<Polynomial_d>`. 

Adapts `Polynomial_traits_d::Invert`. 

For more details see the concept `PolynomialTraits_d::Invert`. 

\sa `Polynomial_d` 
\sa `PolynomialTraits_d` 
\sa `PolynomialTraits_d::Invert` 


*/

template <class Polynomial_d>
Polynomial_traits_d<Polynomial_d>::Invert::result_type 
invert(const Polynomial_d& p, int index = Polynomial_traits_d<Polynomial_d>::d-1);


/*!
\ingroup PkgPolynomialFunctions

For a given `Polynomial_d`, adapts the 
according functor in `Polynomial_traits_d<Polynomial_d>`. 

Adapts `Polynomial_traits_d::Is_square_free`. 

For more details see the concept `PolynomialTraits_d::IsSquareFree`. 

\sa `Polynomial_d` 
\sa `PolynomialTraits_d` 
\sa `PolynomialTraits_d::IsSquareFree` 


*/

template <class Polynomial_d>
Polynomial_traits_d<Polynomial_d>::Is_square_free::result_type 
is_square_free(const Polynomial_d& p);


/*!
\ingroup PkgPolynomialFunctions

For a given `Polynomial_d`, adapts the 
according functor in `Polynomial_traits_d<Polynomial_d>`. 

Adapts `Polynomial_traits_d::Is_zero_at_homogeneous`. 

For more details see the concept `PolynomialTraits_d::IsZeroAtHomogeneous`. 

\sa `Polynomial_d` 
\sa `PolynomialTraits_d` 
\sa `PolynomialTraits_d::IsZeroAtHomogeneous` 


*/

template < class Polynomial_d, class InputIterator >
Polynomial_traits_d<Polynomial_d>::Is_zero_at_homogeneous::result_type
is_zero_at_homogeneous(
const Polynomial_d& p, InputIterator begin, InputIterator end);


/*!
\ingroup PkgPolynomialFunctions

For a given `Polynomial_d`, adapts the 
according functor in `Polynomial_traits_d<Polynomial_d>`. 

Adapts `Polynomial_traits_d::Is_zero_at`. 

For more details see the concept `PolynomialTraits_d::IsZeroAt`. 

\sa `Polynomial_d` 
\sa `PolynomialTraits_d` 
\sa `PolynomialTraits_d::IsZeroAt` 


*/

template < class Polynomial_d, class InputIterator >
Polynomial_traits_d<Polynomial_d>::Is_zero_at::result_type
is_zero_at(
const Polynomial_d& p, InputIterator begin, InputIterator end);


/*!
\ingroup PkgPolynomialFunctions

For a given `Polynomial_d`, adapts the 
according functor in `Polynomial_traits_d<Polynomial_d>`. 

Adapts `Polynomial_traits_d::Leading_coefficient`. 

For more details see the concept `PolynomialTraits_d::LeadingCoefficient`. 

\sa `Polynomial_d` 
\sa `PolynomialTraits_d` 
\sa `PolynomialTraits_d::LeadingCoefficient` 


*/

template <class Polynomial_d>
Polynomial_traits_d<Polynomial_d>::Leading_coefficient::result_type 
leading_coefficient(const Polynomial_d& p);


/*!
\ingroup PkgPolynomialFunctions

For a given `Polynomial_d`, adapts the 
according functor in `Polynomial_traits_d<Polynomial_d>`. 

Adapts `Polynomial_traits_d::Make_square_free`. 

For more details see the concept `PolynomialTraits_d::MakeSquareFree`. 

\sa `Polynomial_d` 
\sa `PolynomialTraits_d` 
\sa `PolynomialTraits_d::MakeSquareFree` 


*/

template <class Polynomial_d>
Polynomial_traits_d<Polynomial_d>::Make_square_free::result_type 
make_square_free(const Polynomial_d& p);


/*!
\ingroup PkgPolynomialFunctions

For a given `Polynomial_d`, adapts the 
according functor in `Polynomial_traits_d<Polynomial_d>`. 

Adapts `Polynomial_traits_d::Move`. 

For more details see the concept `PolynomialTraits_d::Move`. 

\sa `Polynomial_d` 
\sa `PolynomialTraits_d` 
\sa `PolynomialTraits_d::Move` 


*/

template <class Polynomial_d>
Polynomial_traits_d<Polynomial_d>::Move::result_type 
move(const Polynomial_d& p, int i, int j);


/*!
\ingroup PkgPolynomialFunctions

For a given `Polynomial_d`, adapts the 
according functor in `Polynomial_traits_d<Polynomial_d>`. 

Adapts `Polynomial_traits_d::Multivariate_content`. 

For more details see the concept `PolynomialTraits_d::MultivariateContent`. 

\sa `Polynomial_d` 
\sa `PolynomialTraits_d` 
\sa `PolynomialTraits_d::MultivariateContent` 


*/

template <class Polynomial_d>
Polynomial_traits_d<Polynomial_d>::Multivariate_content::result_type 
multivariate_content(const Polynomial_d& p);


/*!
\ingroup PkgPolynomialFunctions

For a given `Polynomial_d`, adapts the 
according functor in `Polynomial_traits_d<Polynomial_d>`. 

Adapts `Polynomial_traits_d::Negate`. 

For more details see the concept `PolynomialTraits_d::Negate`. 

\sa `Polynomial_d` 
\sa `PolynomialTraits_d` 
\sa `PolynomialTraits_d::Negate` 


*/

template <class Polynomial_d>
Polynomial_traits_d<Polynomial_d>::Negate::result_type 
negate(const Polynomial_d& p, int index = Polynomial_traits_d<Polynomial_d>::d-1);


/*!
\ingroup PkgPolynomialFunctions

\brief computes the number of distinct real roots of \f$f\f$.

Given a polynomial \f$ f\f$, or a range of values that is interpreted 
as the principal Sturm-Habicht coefficients of \f$ f\f$, the function computes 
\f[ m:=\# \{\alpha\in\mathbb{R}\mid f(\alpha)=0\} \f] 
that is, the number of distinct real roots of \f$ f\f$. 

The coefficient type of the polynomial, 
or the value type of the iterator range, respectively 
must be a model of `RealEmbeddable`. 
In the second version, 
it is not required to pass the exact princiapl Sturm-Habicht coefficients 
to the functions; it is only required that the sign of each element 
corresponds to the sign of the actual principal Sturm-Habicht coefficient. 

\cgalAdvancedBegin
We explain the internals of this function. 
For a sequence \f$ I:=(a_0,\ldots,a_n)\f$ of real numbers with \f$ a_0\neq 0\f$, define 
\f[ C(I)=\ccSum{i=1}{s}\epsilon_i \f] 
where \f$ s\f$ is the number of subsequences of \f$ I\f$ of the form 

\image html underbrace.png
\image latex underbrace.png

with \f$ a\neq 0,b\neq 0, k\geq 0\f$. 

For the \f$ i\f$-th subsequence of \f$ I\f$, define 

\f[
\epsilon_i:=\begin{array}{cc}
0 & \mbox{if $k$ is odd},\\
(-1)^{k/2}\mathrm{sign}(ab) & \mbox{if $k$ is even}.
\end{array}
\f]

For \f$ f\in\mathbb{R}[x]\f$ with \f$ \deg f=n\f$, we have: 
\f[ C(\mathrm{stha}_n(f),\ldots,\mathrm{stha}_0(f)) = \#\{\alpha\in\R\mid f(\alpha)=0\} \f] 
In other words, the signs of the principal Sturm-Habicht coefficients 
determine the number of distinct real roots of \f$ f\f$. 
\cgalAdvancedEnd

\sa `Polynomial_d` 
\sa `PolynomialTraits_d` 
\sa `PolynomialTraits_d::PrincipalSturmHabichtSequence` 
*/
template<typename Polynomial_d>
int number_of_real_roots(Polynomial_d f);

/*!
\ingroup PkgPolynomialFunctions

\brief computes the number of distinct real roots of \f$ f\f$ whose principal
Sturm-Habicht coefficients are passed by the iterator range.

Given a polynomial \f$ f\f$, or a range of values that is interpreted 
as the principal Sturm-Habicht coefficients of \f$ f\f$, the function computes 
\f[ m:=\# \{\alpha\in\mathbb{R}\mid f(\alpha)=0\} \f] 
that is, the number of distinct real roots of \f$ f\f$. 

The coefficient type of the polynomial, 
or the value type of the iterator range, respectively 
must be a model of `RealEmbeddable`. 
In the second version, 
it is not required to pass the exact princiapl Sturm-Habicht coefficients 
to the functions; it is only required that the sign of each element 
corresponds to the sign of the actual principal Sturm-Habicht coefficient. 

\cgalAdvancedBegin
We explain the internals of this function. 
For a sequence \f$ I:=(a_0,\ldots,a_n)\f$ of real numbers with \f$ a_0\neq 0\f$, define 
\f[ C(I)=\ccSum{i=1}{s}\epsilon_i \f] 
where \f$ s\f$ is the number of subsequences of \f$ I\f$ of the form 

\image html underbrace.png
\image latex underbrace.png

with \f$ a\neq 0,b\neq 0, k\geq 0\f$. 

For the \f$ i\f$-th subsequence of \f$ I\f$, define 

\f[
\epsilon_i:=\begin{array}{cc}
0 & \mbox{if $k$ is odd},\\
(-1)^{k/2}\mathrm{sign}(ab) & \mbox{if $k$ is even}.
\end{array}
\f]

For \f$ f\in\mathbb{R}[x]\f$ with \f$ \deg f=n\f$, we have: 
\f[ C(\mathrm{stha}_n(f),\ldots,\mathrm{stha}_0(f)) = \#\{\alpha\in\R\mid f(\alpha)=0\} \f] 
In other words, the signs of the principal Sturm-Habicht coefficients 
determine the number of distinct real roots of \f$ f\f$. 
\cgalAdvancedEnd

\sa `Polynomial_d` 
\sa `PolynomialTraits_d` 
\sa `PolynomialTraits_d::PrincipalSturmHabichtSequence` 
*/
template<typename InputIterator>
int number_of_real_roots(InputIterator start,InputIterator end);


/*!
\ingroup PkgPolynomialFunctions

For a given `Polynomial_d`, adapts the 
according functor in `Polynomial_traits_d<Polynomial_d>`. 

Adapts `Polynomial_traits_d::Permute`.

For more details see the concept `PolynomialTraits_d::Permute`. 

\sa `Polynomial_d` 
\sa `PolynomialTraits_d` 
\sa `PolynomialTraits_d::Permute` 
*/

template <class Polynomial_d, class InputIterator >
Polynomial_traits_d<Polynomial_d>::Permute::result_type
permute(const Polynomial_d& p, InputIterator begin, InputIterator end );


/*!
\ingroup PkgPolynomialFunctions

computes the polynomial subresultants of \f$ p\f$ and \f$ q\f$, 
with respect to the outermost variable. Each element is of type
`Polynomial_d`.

For a given `Polynomial_d`, adapts the 
according functor in `Polynomial_traits_d<Polynomial_d>`. 

For more details see the concept 
`PolynomialTraits_d::PolynomialSubresultants`. 

\sa `Polynomial_d` 
\sa `PolynomialTraits_d` 
\sa `PolynomialTraits_d::PolynomialSubresultants` 

*/
template<typename Polynomial_d,typename OutputIterator>
OutputIterator polynomial_subresultants
(Polynomial_d p,
Polynomial_d q,
OutputIterator out);


/*!
\ingroup PkgPolynomialFunctions

computes the polynomial subresultants of \f$ p\f$ and \f$ q\f$, 
`sres_out`, with respect to the outermost variable, and
the cofactors for \f$ P\f$, `coP_out` and \f$ Q\f$, `coQ_out`. 
The elements of each output range are of type
`Polynomial_d`.

For a given `Polynomial_d`, adapts the 
according functor in `Polynomial_traits_d<Polynomial_d>`. 

For more details see the concept 
`PolynomialTraits_d::PolynomialSubresultantsWithCofactors`. 

\sa `Polynomial_d` 
\sa `PolynomialTraits_d` 
\sa `PolynomialTraits_d::PolynomialSubresultantsWithCofactors` 

*/
template<typename Polynomial_d,
typename OutputIterator1, 
typename OutputIterator2,
typename OutputIterator3>
OutputIterator1 polynomial_subresultants_with_cofactors
(Polynomial_d p,
Polynomial_d q,
OutputIterator1 sres_out,
OutputIterator2 coP_out,
OutputIterator3 coQ_out);


/*!
\ingroup PkgPolynomialFunctions

computes the principal Sturm-Habicht coefficients of \f$ f\f$
with respect to the outermost variable. Each element is of type
`Polynomial_traits_d::Coefficient_type`b.

For a given `Polynomial_d`, adapts the 
according functor in `Polynomial_traits_d<Polynomial_d>`. 

For more details see the concept 
`PolynomialTraits_d::PrincipalSturmHabichtSequence`. 

\sa `Polynomial_d` 
\sa `PolynomialTraits_d` 
\sa `PolynomialTraits_d::PrincipalSturmHabichtSequence` 

*/
template <typename Polynomial_d,typename OutputIterator> inline
OutputIterator
principal_sturm_habicht_sequence
(typename Polynomial_d f, 
OutputIterator out);


/*!
\ingroup PkgPolynomialFunctions

computes the principal subresultants of \f$ p\f$ and \f$ q\f$, 
with respect to the outermost variable. Each element is of type
`Polynomial_traits_d::Coefficient_type`.

For a given `Polynomial_d`, adapts the 
according functor in `Polynomial_traits_d<Polynomial_d>`. 

For more details see the concept 
`PolynomialTraits_d::PrincipalSubresultants`. 

\sa `Polynomial_d` 
\sa `PolynomialTraits_d` 
\sa `PolynomialTraits_d::PrincipalSubresultants` 

*/
template<typename Polynomial_d,typename OutputIterator>
OutputIterator principal_subresultants
(Polynomial_d p,
Polynomial_d q,
OutputIterator out);


/*!
\ingroup PkgPolynomialFunctions

For a given `Polynomial_d`, adapts the 
according functor in `Polynomial_traits_d<Polynomial_d>`. 

Adapts `Polynomial_traits_d::Pseudo_division`. 

For more details see the concept `PolynomialTraits_d::PseudoDivision`. 

\sa `Polynomial_d` 
\sa `PolynomialTraits_d` 
\sa `PolynomialTraits_d::PseudoDivision` 
*/

template <class Polynomial_d>
void
pseudo_division(
const Polynomial_d& f, const Polynomial_d& g, 
Polynomial_d& q, Polynomial_d& r, Polynomial_traits_d<Polynomial_d>::Coefficient_type& D );


/*!
\ingroup PkgPolynomialFunctions

For a given `Polynomial_d`, adapts the 
according functor in `Polynomial_traits_d<Polynomial_d>`. 

Adapts `Polynomial_traits_d::Pseudo_division_quotient`.

For more details see the concept `PolynomialTraits_d::PseudoDivisionQuotient`. 

\sa `Polynomial_d` 
\sa `PolynomialTraits_d` 
\sa `PolynomialTraits_d::PseudoDivisionQuotient` 



*/

template <class Polynomial_d>
Polynomial_traits_d<Polynomial_d>::Pseudo_division_quotient::result_type 
pseudo_division_quotient(const Polynomial_d& p, const Polynomial_d& q);


/*!
\ingroup PkgPolynomialFunctions

For a given `Polynomial_d`, adapts the 
according functor in `Polynomial_traits_d<Polynomial_d>`. 

Adapts `Polynomial_traits_d::Pseudo_division_remainder`. 

For more details see the concept `PolynomialTraits_d::PseudoDivisionRemainder`. 

\sa `Polynomial_d` 
\sa `PolynomialTraits_d` 
\sa `PolynomialTraits_d::PseudoDivisionRemainder` 

*/

template <class Polynomial_d>
Polynomial_traits_d<Polynomial_d>::Pseudo_division_remainder::result_type 
pseudo_division_remainder(const Polynomial_d& p, const Polynomial_d& q);


/*!
\ingroup PkgPolynomialFunctions

For a given `Polynomial_d`, adapts the 
according functor in `Polynomial_traits_d<Polynomial_d>`. 

Adapts `Polynomial_traits_d::Resultant`. 

For more details see the concept `PolynomialTraits_d::Resultant`. 

\sa `Polynomial_d` 
\sa `PolynomialTraits_d` 
\sa `PolynomialTraits_d::Resultant` 

*/

template <class Polynomial_d>
Polynomial_traits_d<Polynomial_d>::Resultant::result_type 
resultant(const Polynomial_d& p, const Polynomial_d& q);


/*!
\ingroup PkgPolynomialFunctions

For a given `Polynomial_d`, adapts the 
according functor in `Polynomial_traits_d<Polynomial_d>`. 

Adapts `Polynomial_traits_d::Scale_homogeneous`. 

For more details see the concept `PolynomialTraits_d::ScaleHomogeneous`. 

\sa `Polynomial_d` 
\sa `PolynomialTraits_d` 
\sa `PolynomialTraits_d::ScaleHomogeneous` 

*/

template <class Polynomial_d>
Polynomial_traits_d<Polynomial_d>::Scale_homogeneous::result_type 
scale_homogeneous(
const Polynomial_d& p, 
const Polynomial_traits_d<Polynomial_d>::Innermost_coefficient_type& u, 
const Polynomial_traits_d<Polynomial_d>::Innermost_coefficient_type& v, 
int index = Polynomial_traits_d<Polynomial_d>::d-1);


/*!
\ingroup PkgPolynomialFunctions

For a given `Polynomial_d`, adapts the 
according functor in `Polynomial_traits_d<Polynomial_d>`. 

Adapts `Polynomial_traits_d::Scale`. 

For more details see the concept `PolynomialTraits_d::Scale`. 

\sa `Polynomial_d` 
\sa `PolynomialTraits_d` 
\sa `PolynomialTraits_d::Scale` 

*/

template <class Polynomial_d>
Polynomial_traits_d<Polynomial_d>::Scale::result_type 
scale(
const Polynomial_d& p, 
const Polynomial_traits_d<Polynomial_d>::Innermost_coefficient_type& a, 
int index = Polynomial_traits_d<Polynomial_d>::d-1);


/*!
\ingroup PkgPolynomialFunctions

For a given `Polynomial_d`, adapts the 
according functor in `Polynomial_traits_d<Polynomial_d>`. 

Adapts `Polynomial_traits_d::Shift`. 

For more details see the concept `PolynomialTraits_d::Shift`. 

\sa `Polynomial_d` 
\sa `PolynomialTraits_d` 
\sa `PolynomialTraits_d::Shift` 
*/

template <class Polynomial_d>
Polynomial_traits_d<Polynomial_d>::Shift::result_type 
shift(const Polynomial_d& p, int i, int index = Polynomial_traits_d<Polynomial_d>::d-1);


/*!
\ingroup PkgPolynomialFunctions

For a given `Polynomial_d`, adapts the 
according functor in `Polynomial_traits_d<Polynomial_d>`. 

Adapts `Polynomial_traits_d::Sign_at_homogeneous`. 

For more details see the concept `PolynomialTraits_d::SignAtHomogeneous`. 

\sa `Polynomial_d` 
\sa `PolynomialTraits_d` 
\sa `PolynomialTraits_d::SignAtHomogeneous` 


*/

template < class Polynomial_d, class InputIterator >
Polynomial_traits_d<Polynomial_d>::Sign_at_homogeneous::result_type
sign_at_homogeneous(
const Polynomial_d& p, InputIterator begin, InputIterator end);


/*!
\ingroup PkgPolynomialFunctions

For a given `Polynomial_d`, adapts the 
according functor in `Polynomial_traits_d<Polynomial_d>`. 

Adapts `Polynomial_traits_d::Sign_at`. 

For more details see the concept `PolynomialTraits_d::SignAt`. 

\sa `Polynomial_d` 
\sa `PolynomialTraits_d` 
\sa `PolynomialTraits_d::SignAt` 


*/

template < class Polynomial_d, class InputIterator >
Polynomial_traits_d<Polynomial_d>::Sign_at::result_type
sign_at(
const Polynomial_d& p, InputIterator begin, InputIterator end);


/*!
\ingroup PkgPolynomialFunctions

For a given `Polynomial_d`, adapts the 
according functor in `Polynomial_traits_d<Polynomial_d>`. 

Adapts `Polynomial_traits_d::Square_free_factorize`. 

For more details see the concept `PolynomialTraits_d::SquareFreeFactorize`. 

\sa `Polynomial_d` 
\sa `PolynomialTraits_d` 
\sa `PolynomialTraits_d::SquareFreeFactorize` 


*/

template <class Polynomial_d, class OutputIterator > 
OutputIterator
square_free_factorize(
const Polynomial_d& p, 
OutputIterator it,
Polynomial_traits_d<Polynomial>::Innermost_coefficient& a);

/*!
\ingroup PkgPolynomialFunctions

For a given `Polynomial_d`, adapts the 
according functor in `Polynomial_traits_d<Polynomial_d>`. 

Adapts `Polynomial_traits_d::Square_free_factorize`. 

For more details see the concept `PolynomialTraits_d::SquareFreeFactorize`. 

\sa `Polynomial_d` 
\sa `PolynomialTraits_d` 
\sa `PolynomialTraits_d::SquareFreeFactorize` 


*/

template <class Polynomial_d, class OutputIterator > 
OutputIterator
square_free_factorize(const Polynomial_d& p, OutputIterator it);


/*!
\ingroup PkgPolynomialFunctions

For a given `Polynomial_d`, adapts the 
according functor in `Polynomial_traits_d<Polynomial_d>`. 

Adapts `Polynomial_traits_d::Square_free_factorize_up_to_constant_factor`. 

For more details see the concept `PolynomialTraits_d::SquareFreeFactorizeUpToConstantFactor`. 

\sa `Polynomial_d` 
\sa `PolynomialTraits_d` 
\sa `PolynomialTraits_d::SquareFreeFactorizeUpToConstantFactor` 


*/

template <class Polynomial_d, class OutputIterator > 
OutputIterator
square_free_factorize_up_to_constant_factor(const Polynomial_d& p, OutputIterator it);


/*!
\ingroup PkgPolynomialFunctions

computes the Sturm-Habicht-sequence of \f$ f\f$
with respect to the outermost variable. Each element is of type
`Polynomial_d`.

For a given `Polynomial_d`, adapts the 
according functor in `Polynomial_traits_d<Polynomial_d>`. 

For more details see the concept 
`PolynomialTraits_d::SturmHabichtSequence`. 

\sa `Polynomial_d` 
\sa `PolynomialTraits_d` 
\sa `PolynomialTraits_d::SturmHabichtSequence` 

*/
template<typename Polynomial_d,typename OutputIterator> OutputIterator
sturm_habicht_sequence(Polynomial_d f, 
OutputIterator out);


/*!
\ingroup PkgPolynomialFunctions

computes the Sturm-Habicht sequence of \f$ f\f$
`stha_out`, with respect to the outermost variable, and
the cofactors for \f$ f\f$, `cof_out` and \f$ f'\f$, `cofx_out`. 
The elements of each output range are of type
`Polynomial_d`.

For a given `Polynomial_d`, adapts the 
according functor in `Polynomial_traits_d<Polynomial_d>`. 

For more details see the concept 
`PolynomialTraits_d::SturmHabichtSequenceWithCofactors`. 

\sa `Polynomial_d` 
\sa `PolynomialTraits_d` 
\sa `PolynomialTraits_d::SturmHabichtSequenceWithCofactors` 
*/
template<typename Polynomial_d,
typename OutputIterator1,
typename OutputIterator2,
typename OutputIterator3> 
OutputIterator1
sturm_habicht_sequence_with_cofactors
(Polynomial_d f,
OutputIterator1 stha_out,
OutputIterator2 cof_out,
OutputIterator3 cofx_out);


/*!
\ingroup PkgPolynomialFunctions

For a given `Polynomial_d`, adapts the 
according functor in `Polynomial_traits_d<Polynomial_d>`. 

Adapts `Polynomial_traits_d::Substitute_homogeneous`. 

For more details see the concept `PolynomialTraits_d::SubstituteHomogeneous`. 

\sa `Polynomial_d` 
\sa `PolynomialTraits_d` 
\sa `PolynomialTraits_d::SubstituteHomogeneous` 

*/

template < class Polynomial_d, class InputIterator >
CGAL::Coercion_traits<
Polynomial_traits_d<Polynomial_d>::Innermost_coefficient,
std::iterator_traits<Input_iterator>::value_type
>::Type
substitute_homogeneous(
const Polynomial_d& p, InputIterator begin, InputIterator end);


/*!
\ingroup PkgPolynomialFunctions

For a given `Polynomial_d`, adapts the 
according functor in `Polynomial_traits_d<Polynomial_d>`. 

Adapts `Polynomial_traits_d::Substitute`. 

For more details see the concept `PolynomialTraits_d::Substitute`. 

\sa `Polynomial_d` 
\sa `PolynomialTraits_d` 
\sa `PolynomialTraits_d::Substitute` 


*/

template < class Polynomial_d, class InputIterator >
CGAL::Coercion_traits<
Polynomial_traits_d<Polynomial_d>::Innermost_coefficient,
std::iterator_traits<Input_iterator>::value_type
>::Type
substitute(
const Polynomial_d& p, InputIterator begin, InputIterator end);


/*!
\ingroup PkgPolynomialFunctions

For a given `Polynomial_d`, adapts the 
according functor in `Polynomial_traits_d<Polynomial_d>`. 

Adapts `Polynomial_traits_d::Swap`. 

For more details see the concept `PolynomialTraits_d::Swap`. 

\sa `Polynomial_d` 
\sa `PolynomialTraits_d` 
\sa `PolynomialTraits_d::Swap` 


*/

template <class Polynomial_d>
Polynomial_traits_d<Polynomial_d>::Swap::result_type 
swap(const Polynomial_d& p, int i, int j);


/*!
\ingroup PkgPolynomialFunctions

For a given `Polynomial_d`, adapts the 
according functor in `Polynomial_traits_d<Polynomial_d>`. 

Adapts `Polynomial_traits_d::Total_degree`. 

For more details see the concept `PolynomialTraits_d::TotalDegree`. 

\sa `Polynomial_d` 
\sa `PolynomialTraits_d` 
\sa `PolynomialTraits_d::TotalDegree` 


*/

template <class Polynomial_d>
Polynomial_traits_d<Polynomial_d>::Total_degree::result_type 
total_degree(const Polynomial_d& p);


/*!
\ingroup PkgPolynomialFunctions

For a given `Polynomial_d`, adapts the 
according functor in `Polynomial_traits_d<Polynomial_d>`. 

Adapts `Polynomial_traits_d::Translate_homogeneous`. 

For more details see the concept `PolynomialTraits_d::TranslateHomogeneous`. 

\sa `Polynomial_d` 
\sa `PolynomialTraits_d` 
\sa `PolynomialTraits_d::TranslateHomogeneous` 


*/

template <class Polynomial_d>
Polynomial_traits_d<Polynomial_d>::Translate_homogeneous::result_type 
translate_homogeneous(
const Polynomial_d& p, 
const Polynomial_traits_d<Polynomial_d>::Innermost_coefficient_type& u, 
const Polynomial_traits_d<Polynomial_d>::Innermost_coefficient_type& v, 
int index = Polynomial_traits_d<Polynomial_d>::d-1);


/*!
\ingroup PkgPolynomialFunctions

For a given `Polynomial_d`, adapts the 
according functor in `Polynomial_traits_d<Polynomial_d>`. 

Adapts `Polynomial_traits_d::Translate`. 

For more details see the concept `PolynomialTraits_d::Translate`. 

\sa `Polynomial_d` 
\sa `PolynomialTraits_d` 
\sa `PolynomialTraits_d::Translate` 


*/

template <class Polynomial_d>
Polynomial_traits_d<Polynomial_d>::Translate::result_type 
translate(
const Polynomial_d& p, 
const Polynomial_traits_d<Polynomial_d>::Innermost_coefficient_type& a, 
int index = Polynomial_traits_d<Polynomial_d>::d-1);


/*!
\ingroup PkgPolynomialFunctions

For a given `Polynomial_d`, adapts the 
according functor in `Polynomial_traits_d<Polynomial_d>`. 

Adapts `Polynomial_traits_d::UnivariateContent`. 

For more details see the concept `PolynomialTraits_d::UnivariateContent`. 

\sa `Polynomial_d` 
\sa `PolynomialTraits_d` 
\sa `PolynomialTraits_d::Univariate_Content` 


*/

template <class Polynomial_d>
Polynomial_traits_d<Polynomial_d>::Univariate_content::result_type 
univariate_content(const Polynomial_d& p);


/*!
\ingroup PkgPolynomialFunctions

For a given `Polynomial_d`, adapts the 
according functor in `Polynomial_traits_d<Polynomial_d>`. 

Adapts `Polynomial_traits_d::Univariate_content_up_to_constant_factor`. 

For more details see the concept `PolynomialTraits_d::UnivariateContentUpToConstantFactor`. 

\sa `Polynomial_d` 
\sa `PolynomialTraits_d` 
\sa `PolynomialTraits_d::UnivariateContentUpToConstantFactor` 
*/

template <class Polynomial_d>
Polynomial_traits_d<Polynomial_d>::Univariate_content_up_to_constant_factor::result_type 
univariate_content_up_to_constant_factor(const Polynomial_d& p);

} /* namespace CGAL */

