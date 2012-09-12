
namespace CGAL {

/*!
\ingroup PkgArrangement2

\anchor arr_refarr_dcel_base 

The `Arr_dcel_base` class is an important ingredient in the 
definition of <span class="textsc">Dcel</span> data structures. It serves as a basis class for 
any instance of the `Dcel` template parameter of the 
`Arrangement_2` template. In particular it is the basis class of 
the default `Dcel` template parameter, and the basis class of any 
extended <span class="textsc">Dcel</span>. The template parameters `V`, `H`, and `F` 
must be instantiated with models of the concepts 
`ArrangementDcelVertex`, `ArrangementDcelHalfedge`, 
and `ArrangementDcelFace` respectively. 

\models ::ArrangementDcel 
CONVERRORIsModel: CONVERROR 3 nested classes missing 

*/
template< typename V, typename H, typename F >
class Arr_dcel_base {
public:

/// @}

}; /* end Arr_dcel_base */
} /* end namespace CGAL */
