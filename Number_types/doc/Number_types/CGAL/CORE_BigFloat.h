
namespace CORE {

/*!
\ingroup nt_core

The class `CORE::BigFloat` is a variable precision floating-point type. 
Rounding mode and precision (i.e.\ mantissa length) of 
`CORE::BigFloat` can be set. 
Since it also carries the error of a computed value. 

This number type is provided by the <span class="textsc">Core</span> library \cgalCite{klpy-clp-99}. 

\cgal defines the necessary functions so that this class complies to the 
requirements on number types. 

\cgalModels `FieldWithKthRoot` 
\cgalModels `RealEmbeddable` 
\cgalModels `FromDoubleConstructible` 

*/

class BigFloat {
}; /* end CORE::BigFloat */
} /* end namespace CORE */
