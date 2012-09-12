
namespace CGAL {

/*!
\ingroup PkgArrangement2

The class `Arr_curve_data_traits_2` is a model of the `ArrangementTraits_2` concept 
and serves as a decorator class that allows the extension of the curves 
defined by the base traits-class (the `Tr` parameter), which serves as a 
geometric traits-class (a model of the `ArrangementTraits_2` concept), with 
extraneous (non-geometric) data fields. 

The traits class inherits its point type from `Traits::Point_2`, 
and defines an extended `Curve_2` and `X_monotone_curve_2` types, 
as detailed below. 

Each `Curve_2` object is associated with a single data field of type 
`CData`, and each `X_monotone_curve_2` object is associated with 
a single data field of type `XData`. When a curve is 
subdivided into \f$ x\f$-monotone subcurves, its data field is converted using 
the conversion functor, which is specified by the `Cnv` template-parameter, 
and the resulting objects is copied to all `X_monotone_curve_2` objects 
induced by this curve. The conversion functor should provide an operator with 
the following prototype: 

`XData operator() (const CData& d) const;` 

By default, the two data types are the same, so the conversion operator 
is trivial: 

<TABLE><TR><TD ALIGN=LEFT VALIGN=TOP NOWRAP> 

`CData` = 
<TD ALIGN=LEFT VALIGN=TOP NOWRAP> 
`XData` 
<TR><TD ALIGN=LEFT VALIGN=TOP NOWRAP> 
`Cnv` = 
<TD ALIGN=LEFT VALIGN=TOP NOWRAP> 
`_Default_convert_functor<CData,XData>` 

</TABLE> 

In case two (or more) \f$ x\f$-monotone curves overlap, their data fields are 
merged to a single field, using the merge functor functor, which is 
specified by the `Mrg` template-parameter. This functor should provide 
an operator with the following prototype: 

`XData operator() (const XData& d1, const XData& d2) const;` 

which returns a single data object that represents the merged data field 
of `d1` and `d2`. The \f$ x\f$-monotone curve that represents the overlap 
is associated with the output of this functor. 

\models ::ArrangementTraits_2 

Inherits From 
-------------- 

CONVERROR Inherits From must be handled manually, e.g. adjust the class decl`Base_traits_2` 

CONVERROR 2 nested classes missing 

*/
template< typename Tr, typename XData, typename Mrg, typename CData, typename Cnv >
class Arr_curve_data_traits_2 {
public:

/// \name Types 
/// @{

/*! 
the base traits-class. 
*/ 
typedef Tr Base_traits_2; 

/*! 
the base curve. 
*/ 
typedef typename Base_traits_2::Curve_2 Base_curve_2; 

/*! 
the base \f$ x\f$-monotone curve curve. 
*/ 
typedef typename Base_traits_2::X_monotone_curve_2 Base_x_monotone_curve_2; 

/*! 
the merge functor. 
*/ 
typedef Mrg Merge; 

/*! 
the conversion functor. 
*/ 
typedef Cnv Convert; 

/*! 
the type of data associated with curves. 
*/ 
typedef CData Curve_data; 

/*! 
the type of data associated with \f$ x\f$-monotone curves. 
*/ 
typedef XData X_monotone_curve_data; 

/// @}

}; /* end Arr_curve_data_traits_2 */
} /* end namespace CGAL */
