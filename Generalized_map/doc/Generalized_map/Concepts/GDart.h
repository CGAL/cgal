
/*!
\ingroup PkgGeneralizedMapsConcepts
\cgalConcept

The concept `GDart` defines a <I>d</I>-dimensional dart for generalized maps. A dart mainly
stores handles to the darts linked with itself by \f$ \alpha_i\f$, \f$ \forall\f$<I>i</I>,
0\f$ \leq\f$<I>i</I>\f$ \leq\f$<I>d</I>. Moreover, it stores also handles to each
non void attribute associated with itself.

\cgalHeading{Creation}

A dart `d0` is never constructed directly, but always created within a generalized map `cm` by using the method \ref GeneralizedMap::create_dart "cm.create_dart()". A new dart is initialized to be <I>i</I>-free, \f$ \forall\f$<I>i</I>: 0\f$ \leq\f$<I>i</I>\f$ \leq\f$<I>dimension</I>, and having all its attribute handles initialized to `NULL`, for each non `void` attribute.

\cgalHasModel \ref CGAL::GDart "CGAL::GDart<d,GMap>"
*/

class GDart {
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

}; /* end GDart */
