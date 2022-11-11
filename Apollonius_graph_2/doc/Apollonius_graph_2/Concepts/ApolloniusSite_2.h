
/*!
\ingroup PkgApolloniusGraph2Concepts
\cgalConcept

The concept `ApolloniusSite_2` provides the requirements for an
Apollonius site class.

\sa `ApolloniusGraphTraits_2`
\sa `CGAL::Apollonius_site_2<K>`
\sa `CGAL::Apollonius_graph_traits_2<K,Method_tag>`
\sa `CGAL::Apollonius_graph_filtered_traits_2<CK,CM,EK,EM,FK,FM>`

*/

class ApolloniusSite_2 {
public:

/// \name Types
/// @{

/*!
The point type.
*/
typedef unspecified_type Point_2;

/*!
The field number type.
*/
typedef unspecified_type FT;

/*!
The ring number type.
*/
typedef unspecified_type RT;

/*!
The weight type.
\pre It must be the same as `FT`.
*/
typedef unspecified_type Weight;

/// @}

/// \name Creation
/// @{

/*!

*/
ApolloniusSite2(Point_2 p=Point_2(), Weight w= Weight(0));

/*!
Copy constructor.
*/
ApolloniusSite_2(ApolloniusSite_2 other);

/// @}

/// \name Access Functions
/// @{

/*!
Returns the center of the Apollonius
site.
*/
Point_2 point() const;

/*!
Returns the weight of the
Apollonius site.
*/
Weight weight() const;

/// @}

}; /* end ApolloniusSite_2 */

