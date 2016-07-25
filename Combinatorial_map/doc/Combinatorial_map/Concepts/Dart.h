
/*!
\ingroup PkgCombinatorialMapsConcepts
\cgalConcept

The concept `Dart` defines a <I>d</I>-dimensional dart. A dart mainly
stores handles to the darts linked with itself by \f$ \beta_i\f$, \f$ \forall\f$<I>i</I>,
0\f$ \leq\f$<I>i</I>\f$ \leq\f$<I>d</I>. Moreover, it stores also handles to each
non void attribute associated with itself.

\cgalHeading{Creation}

A dart `d0` is never constructed directly, but always created
within a combinatorial map `cm` by using the method
\ref CombinatorialMap::create_dart "cm.create_dart()". A new dart is initialized to be <I>i</I>-free,
\f$ \forall\f$<I>i</I>: 0\f$ \leq\f$<I>i</I>\f$ \leq\f$<I>dimension</I>, and having all
its attribute handles initialized to `NULL`, for each non `void` attribute.

\cgalHasModel \ref CGAL::Dart "CGAL::Dart<d,CMap>"

*/

class Dart {
public:

/// \name Constants
/// @{

/*!
The dimension of the dart.
*/
static unsigned int dimension;

/// @}

/// \name Types
/// @{

/*!
%Dart handle type. Must be a model of `Handle` concept.
*/
typedef unspecified_type Dart_handle;

/*!
%Dart const handle type. Must be a model of `ConstHandle` concept.
*/
typedef unspecified_type Dart_const_handle;

/*!
`Attribute_handle<i>::%type` is a handle to `i`-attributes, with 0 \f$ \leq \f$ `i` \f$ \leq \f$ `dimension`. Must be a model of `Handle` concept.
\note It can be implemented using a nested template class.
*/
template <unsigned int i>
using Attribute_handle = unspecified_type;

/*!
`Attribute_const_handle<i>::%type` is a const handle to `i`-attributes, with 0 \f$ \leq \f$ `i` \f$ \leq \f$ `dimension`. Must be a model of `ConstHandle` concept.
\note It can be implemented using a nested template class.
*/
template <unsigned int i>
using Attribute_const_handle = unspecified_type;

/// @}

/// \name Access Member Functions
/// @{

/*!
To simplify a future implementation, it is recommended to not use this function and to use \ref CombinatorialMap::beta "cmap.beta(dh,i)" instead.
Returns \f$ \beta_i\f$(`*this`).
\pre 0\f$ \leq\f$<I>i</I>\f$ \leq\f$<I>dimension</I>.
*/
Dart_handle beta(unsigned int i);

/*!
To simplify a future implementation, it is recommended to not use this function and to use \ref CombinatorialMap::beta "cmap.beta(dh,i)" instead.
Returns \f$ \beta_i\f$(`*this`) when the dart is const.
\pre 0\f$ \leq\f$<I>i</I>\f$ \leq\f$<I>dimension</I>.
*/
Dart_const_handle beta(unsigned int i) const;

/*!
To simplify a future implementation, it is recommended to not use this function and to use \ref CombinatorialMap::beta "cmap.beta(dh,CGAL_BETAINV(i))" instead.
Returns \f$ \beta_i^{-1}\f$(`*this`).
\pre 0\f$ \leq\f$<I>i</I>\f$ \leq\f$<I>dimension</I>.
*/
Dart_handle beta_inv(unsigned int i);

/*!
To simplify a future implementation, it is recommended to not use this function and to use \ref CombinatorialMap::beta "cmap.beta(dh,CGAL_BETAINV(i))" instead.
Returns \f$ \beta_i^{-1}\f$(`*this`) when the dart is const.
\pre 0\f$ \leq\f$<I>i</I>\f$ \leq\f$<I>dimension</I>.
*/
Dart_const_handle beta_inv(unsigned int i) const;

/*!
To simplify a future implementation, it is recommended to not use this function and to use \ref CombinatorialMap::attribute "cmap.attribute(dh)" instead.
Returns a handle to the <I>i</I>-attribute associated to the dart.
\pre 0\f$ \leq\f$<I>i</I>\f$ \leq\f$<I>dimension</I>, and <I>i</I>-attributes are non `void`.
*/
template <unsigned int i> Attribute_handle<i>::type attribute();

/*!
To simplify a future implementation, it is recommended to not use this function and to use \ref CombinatorialMap::attribute "cmap.attribute(dh)" instead.
Returns a const handle to the <I>i</I>-attribute associated to the dart,
when the dart is const.
\pre 0\f$ \leq\f$<I>i</I>\f$ \leq\f$<I>dimension</I>, and <I>i</I>-attributes are non `void`.
*/
template <unsigned int i>
Attribute_const_handle<i>::type attribute() const;

/// @}

}; /* end Dart */
