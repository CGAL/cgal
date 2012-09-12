
namespace CGAL {

/*!
\ingroup PkgArrangement2

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

\models ::ArrangementTraits_2 

CONVERROR 1 missing nested class 
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
typedef Hidden_type typedef Data_container; 

/*! 
a non-mutable iterator for the data objects in the data container. 
*/ 
typedef Hidden_type typedef Data_iterator; 

/// @}

}; /* end Arr_consolidated_curve_data_traits_2 */
} /* end namespace CGAL */
