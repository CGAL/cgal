
namespace CGAL {

/*!
\ingroup PkgApolloniusGraph2Ref

The class `Apollonius_graph_filtered_traits_2` provides a model for the
`ApolloniusGraphTraits_2` concept.

The class `Apollonius_graph_filtered_traits_2` uses the filtering technique \cgalCite{cgal:bbp-iayed-01}
to achieve traits for the `Apollonius_graph_2<Gt,Agds>` class
with efficient and exact predicates given an exact
kernel `EK` and a filtering kernel `FK`. The geometric
constructions associated provided by this class are equivalent
to those provided by the traits class
`Apollonius_graph_traits_2<CK,CM>`, which means that they may
be inexact.

This class has six template parameters. The first, third and fifth
template parameters must be a models of the `Kernel` concept. The
second, fourth and sixth template parameters correspond to how
predicates are evaluated. There are two predefined possible values for
`Method_tag`, namely `CGAL::Field_with_sqrt_tag` and
`CGAL::Integral_domain_without_division_tag`. The first one must be used when the number type
used in the representation supports the exact evaluation of signs of
expressions involving all four basic operations and square roots,
whereas the second one requires the exact evaluation of signs of
ring-type expressions, i.e., expressions involving only additions,
subtractions and multiplications.
The way the predicates are evaluated is discussed in
\cgalCite{cgal:ke-ppawv-02}, \cgalCite{cgal:ke-rctac-03}.

The default values for the template parameters are as follows:
`CM = CGAL::Integral_domain_without_division_tag`,
`EK = CGAL::Simple_cartesian<CGAL::MP_Float>`,
`EM = CM`,
`FK = CGAL::Simple_cartesian<CGAL::Interval_nt<false> >`,
`FM = CM`.

\cgalModels{ApolloniusGraphTraits_2}

\sa `Kernel`
\sa `ApolloniusGraphTraits_2`
\sa `CGAL::Integral_domain_without_division_tag`
\sa `CGAL::Field_with_sqrt_tag`
\sa `CGAL::Apollonius_graph_2<Gt,Agds>`
\sa `CGAL::Apollonius_graph_traits_2<K,Method_tag>`

*/
template< typename CK, typename CM, typename EK, typename EM, typename FK, typename FM >
class Apollonius_graph_filtered_traits_2 {
public:

/// \name Creation
/// @{

/*!
%Default constructor.
*/
Apollonius_graph_filtered_traits_2<CK,CM,EK,EM,FK,FM>();

/*!
Copy
constructor.
*/
Apollonius_graph_filtered_traits_2
(Apollonius_graph_filtered_traits_2<CK,CM,EK,EM,FK,FM> other);

/*!
Assignment
operator.
*/
Apollonius_graph_filtered_traits_2<CK,CM,EK,EM,FK,FM> operator=(other);

/// @}

}; /* end Apollonius_graph_filtered_traits_2 */
} /* end namespace CGAL */
