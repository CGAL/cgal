
/*!
\ingroup nt_leda

The class `leda_real` is a wrapper class that provides the functions needed 
to use the number type `real`, representing exact real numbers 
numbers provided by \leda. 
The class `leda_real` provides exact computation over the subset of real 
numbers that contains integers, and which is closed by the operations 
\f$ +,-,\times,/,\sqrt{}\f$ and \f$\sqrt[k]{}\f$. For \leda version 5.0 or later 
`leda_real` is also able to represent real roots of polynomials. 
Operations and comparisons between objects of this type are guaranteed 
to be exact. 

\cgalModels `FieldWithRootOf` 
\cgalModels `RealEmbeddable` 
\cgalModels `FromDoubleConstructible` 

For more details on the number types of \leda we refer to the \leda manual \cgalCite{cgal:mnsu-lum}. 

*/

class leda_real {
}; /* end leda_real */

