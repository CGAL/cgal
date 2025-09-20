
/*!
\ingroup PkgAlgebraicFoundationsAlgebraicStructuresConcepts
\cgalConcept

A model of the concept `FromDoubleConstructible` is required
to be constructible from the type `double`.

In case the type is a model of `RealEmbeddable` too, for any double d
the identity: `d == ` \link CGAL::to_double ` CGAL::to_double(T(d))`\endlink, is guaranteed.

*/

class FromDoubleConstructible {
public:

/// \name Creation
/// @{

/*!

conversion constructor from double.

*/
FromDoubleConstructible(const double& d);

/// @}

}; /* end FromDoubleConstructible */
