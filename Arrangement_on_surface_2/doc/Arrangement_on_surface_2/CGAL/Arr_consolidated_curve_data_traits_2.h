
namespace CGAL {

/*!
\ingroup PkgArrangementOnSurface2TraitsClasses

The class `Arr_consolidated_curve_data_traits_2` is a model of the concept `ArrangementTraits_2`,
and serves as a decorator class that enables the extension of the curve
type defined by the `Traits` parameter. The traits class inherits its
point type from `Traits::Point_2`, and defines the types
`Curve_2` and `X_monotone_curve_2` extended with extraneous data
fields of type `Data`.

Each `Curve_2` object is associated with a single data field of type
`Data`, and each `X_monotone_curve_2` object is associated with
a set of unique data objects. When a curve is subdivided into \f$ x\f$-monotone
subcurves, all resulting subcurves are associated with a list containing
a single data object, copied from the inducing curve. When an \f$ x\f$-monotone
curve is split, its data set is duplicated, and inserted into the sets of
both resulting subcurves. In case two (or more) \f$ x\f$-monotone curves
overlap, their data sets are consolidated, and are inserted into the set
of the \f$ x\f$-monotone curve that represents the overlap.

\cgalModels `ArrangementTraits_2`

*/
template< typename Traits, typename Data >
class Arr_consolidated_curve_data_traits_2
  : public Arr_curve_data_traits_2<Traits, _Unique_list<Data>,
                                   _Consolidate_unique_lists<Data>,
                                   Data>
{
public:

/// \name Types
/// @{

/*!
the base traits-class.
*/
typedef Traits Base_traits_2;

/*!
the base curve.
*/
typedef typename Base_traits_2::Curve_2 Base_curve_2;

/*!
the base \f$ x\f$-monotone curve curve.
*/
typedef typename Base_traits_2::X_monotone_curve_2 Base_x_monotone_curve_2;

/*!
a set of data objects that is associated with an \f$ x\f$-monotone curve.
*/
typedef unspecified_type typedef Data_container;

/*!
a non-mutable iterator for the data objects in the data container.
*/
typedef unspecified_type typedef Data_iterator;

/// @}


/*!


The `Data_container` class nested within the consolidated
curve-data traits and associated with the `Traits::X_monotone_curve_2`
type is maintained as a list with unique data objects. This representation is
simple and efficient in terms of memory consumption. It also requires that
the `Data` class supports only the equality operator. Note however that
most set operations require linear time.

*/
class Data_container {
public:

/// \name Creation
/// @{

/*!
default constructor.
*/
Data_container ();

/*!
constructs set containing a single `data` object.
*/
Data_container (const Data& data);

/// @}

/// \name Access Functions
/// @{

/*!
returns the number of data objects in the set.
*/
std::size_t size () const;

/*!
returns an iterator pointing to the first data object.
*/
Data_iterator begin () const;

/*!
returns a past-the-end iterator for the data objects.
*/
Data_iterator end () const;

/*!
returns the first data object inserted into the set.
\pre The number of data objects is not \f$ 0\f$.
*/
const Data& front () const;

/*!
returns the last data object inserted into the set.
\pre The number of data objects is not \f$ 0\f$.
*/
const Data& back () const;

/// @}

/// \name Predicates
/// @{

/*!
check if the two sets contain the same data objects (regardless of order).
*/
bool operator== (const Data_container& other) const;

/*!
find the given `data` object in the set and returns an iterator
for this object, or `end()` if it is not found.
*/
Data_iterator find (const Data& data);

/// @}

/// \name Modifiers
/// @{

/*!
inserts the given `data` object into the set. Returns `true` on
success, or `false` if the set already contains the object.
*/
bool insert (const Data& data);

/*!
erases the given `data` object from the set. Returns `true` on
success, or `false` if the set does not contain the object.
*/
bool erase (const Data& data);

/*!
clears the set.
*/
void clear ();

/// @}

}; /* end Arr_consolidated_curve_data_traits_2::Data_container */



}; /* end Arr_consolidated_curve_data_traits_2 */
} /* end namespace CGAL */
