
namespace CGAL {

/*!
\ingroup PkgApolloniusGraph2Ref

The class `Apollonius_graph_traits_2` provides a model for the
`ApolloniusGraphTraits_2` concept.
This class has two template parameters. The first template parameter
must be a model of the `Kernel` concept. The second template
parameter corresponds to how predicates are evaluated. There are two
predefined possible values for `Method_tag`, namely
`CGAL::Field_with_sqrt_tag` and `CGAL::Integral_domain_without_division_tag`. The first one
must be used when the number type used in the representation supports
the exact evaluation of signs of expressions involving all four basic
operations and square roots, whereas the second one requires the exact
evaluation of signs of ring-type expressions, i.e., expressions
involving only additions, subtractions and multiplications. The
default value for `Method_tag` is `CGAL::Integral_domain_without_division_tag`.
The way the predicates are evaluated is discussed in
\cgalCite{cgal:ke-ppawv-02}, \cgalCite{cgal:ke-rctac-03}.

\cgalModels `ApolloniusGraphTraits_2`

\sa `CGAL::Apollonius_graph_2<Gt,Agds>`
\sa `CGAL::Apollonius_graph_filtered_traits_2<CK,CM,EK,EM,FK,FM>`
*/
template< typename K, typename Method_tag >
class Apollonius_graph_traits_2 {
public:

/// \name Creation
/// @{

/*!
%Default constructor.
*/
Apollonius_graph_traits_2<K,Method_tag>();

/*!
Copy constructor.
*/
Apollonius_graph_traits_2<K,Method_tag>(const Apollonius_graph_traits_2<K,Method_tag>& other);

/*!
Assignment operator.
*/
Apollonius_graph_traits_2<K,Method_tag>
operator=(const Apollonius_graph_traits_2<K,Method_tag>& other);

/// @}

}; /* end Apollonius_graph_traits_2 */
} /* end namespace CGAL */
