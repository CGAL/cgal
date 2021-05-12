
namespace CGAL {

/*!
\ingroup PkgCombinatorialMapsClasses

The class `Dart` represents a <I>d</I>D dart in a combinatorial map.

\deprecated This class is deprecated since CGAL 4.9. Dart is now a type defined internally; users can now only define the information associated with darts. All functions defined in this class are now defined as methods of a combinatorial map taking a `Dart_handle` as first parameter. `CGAL_CMAP_DART_DEPRECATED` can be defined to keep the old behavior.

*/
template< typename d, typename CMap >
struct Dart {
/*!

*/
typedef CMap::Dart_handle Dart_handle;

/*!

*/
typedef CMap::Dart_const_handle Dart_const_handle;

/*!

*/
template <unsigned int i>
using Attribute_handle = CMap::Attribute_handle<i>;

/*!

*/
template <unsigned int i>
using Attribute_const_handle = CMap::Attribute_const_handle<i>;

/*!
Returns \f$ \beta_i\f$(`*this`).
\pre 0 \f$ \leq \f$ <I>i</I> \f$ \leq \f$ <I>dimension</I>.
*/
Dart_handle beta(unsigned int i);

/*!
Returns \f$ \beta_i\f$(`*this`) when the dart is const.
\pre 0 \f$ \leq \f$ <I>i</I> \f$ \leq \f$ <I>dimension</I>.
*/
Dart_const_handle beta(unsigned int i) const;

/*!
Returns \f$ \beta_i^{-1}\f$(`*this`).
\pre 0 \f$ \leq \f$ <I>i</I> \f$ \leq \f$ <I>dimension</I>.
*/
Dart_handle beta_inv(unsigned int i);

/*!
Returns \f$ \beta_i^{-1}\f$(`*this`) when the dart is const.
\pre 0 \f$ \leq \f$ <I>i</I> \f$ \leq \f$ <I>dimension</I>.
*/
Dart_const_handle beta_inv(unsigned int i) const;

}; /* end Dart */
} /* end namespace CGAL */
