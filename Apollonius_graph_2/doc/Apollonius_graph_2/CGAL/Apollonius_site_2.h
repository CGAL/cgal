
namespace CGAL {

/*!
\ingroup PkgApolloniusGraph2Ref

The class `Apollonius_site_2` is a model for the concept
`ApolloniusSite_2`. It is parametrized by a template parameter
`K` which must be a model of the `Kernel` concept.

\cgalModels{ApolloniusSite_2}

\cgalHeading{Types}

The class `Apollonius_site_2` does not introduce any types in addition to the
concept `ApolloniusSite_2`.

\cgalHeading{I/O}

The I/O operators are defined for `std::iostream`.

The information output in the `std::iostream` is: the point of the
Apollonius site and its weight.

\sa `CGAL::Apollonius_graph_traits_2<K,Method_tag>`
\sa `CGAL::Apollonius_graph_filtered_traits_2<CK,CM,EK,EM,FK,FM>`
*/
template< typename K >
class Apollonius_site_2 {
public:

/// \name Creation
/// @{

/*!

*/
Apollonius_site_2(Point_2 p=Point_2(), Weight w= Weight(0));

/*!
Copy constructor.
*/
Apollonius_site_2(const Apollonius_site_2<K>& other);

/// @}

}; /* end Apollonius_site_2 */

/*!
Inserts the
Apollonius site `s` into the stream `os`.
\pre The insert operator must be defined for `Point_2` and `Weight`.
\relates Apollonius_site_2
*/
std::ostream& operator<<(std::ostream& os, const Apollonius_site_2<K>& s) const;

/*!
Reads an Apollonius site from the stream `is` and assigns it
to `s`.
\pre The extract operator must be defined for `Point_2` and `Weight`.
\relates Apollonius_site_2
*/
std::istream& operator>>(std::istream& is, const Apollonius_site_2<K>& s);

} /* end namespace CGAL */
