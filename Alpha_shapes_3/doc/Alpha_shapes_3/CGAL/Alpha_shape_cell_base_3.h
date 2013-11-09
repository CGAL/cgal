
namespace CGAL {

/*!
\ingroup PkgAlphaShapes3

The class `Alpha_shape_cell_base_3` is the default model for the concept 
`AlphaShapeCell_3`. 

\tparam Traits is the geometric traits class that is provided 
to the `Alpha_shape_3` class. 
\tparam Cb must be a cell base class instantiated by default 
with `Triangulation_cell_base_3<Traits>`. 
\tparam ExactAlphaComparisonTag is a tag that, when set to 
`Tag_true`, triggers exact comparisons between alpha values. See the description 
provided in the documentation of `Alpha_shape_3` for more details. The default value is `Tag_false`. 
\tparam WeightedTag is used only if `ExactAlphaComparisonTag` is `Tag_true`. It 
must be `Tag_true` if the underlying triangulation of the alpha shape to be used is a Regular triangulation 
and `Tag_false` otherwise. The default is `Tag_false`. 

\cgalModels `AlphaShapeCell_3`

*/
template< typename Traits, typename Cb, typename ExactAlphaComparisonTag, typename WeightedTag >
class Alpha_shape_cell_base_3 : public Fb {
public:

/// @}

}; /* end Alpha_shape_cell_base_3 */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgAlphaShapes3

The class `Alpha_status` is a small data structure to store 
the critical alpha values of faces of an alpha shape. 
Each face has three critical alpha values, called 
`alpha_min`, `alpha_mid` and `alpha_max` in increasing order. 
The face will be exterior for any \f$ \alpha < \f$ `alpha_min`, 
singular for `alpha_min` \f$ \leq \alpha < \f$ `alpha_mid`, 
regular for `alpha_mid` \f$ \leq \alpha < \f$ `alpha_max` 
and interior for `alpha_max` \f$ \leq \alpha\f$. 
The value `alpha_min` is undefined for faces which are not Gabriel 
faces and therefore do not appear in the alpha complex 
without any of their 
including face. The value `alpha_max` is undefined 
for convex hull faces which can never be interior.
The data structure also includes two Boolean to mark 
if the face is a Gabriel face or a convex hull face. 

The class `Alpha_status` is parameterized by a number type `NT`. 

\sa `::AlphaShapeCell_3`
\sa `::AlphaShapeVertex_3`

*/
template< typename NT >
class Alpha_status {
public:

/// \name Creation 
/// @{

/*!
default constructor. 
*/ 
Alpha_status(); 

/// @} 

/// \name Modifiers 
/// @{

/*!
sets Gabriel marker. 
*/ 
void set_is_Gabriel(bool yesorno); 

/*!
sets convex hull marker.
*/ 
void set_is_on_chull(bool yesorno); 

/*!
sets `alpha_min`. 
*/ 
void set_alpha_min(NT alpha); 

/*!
sets `alpha_mid`. 
*/ 
void set_alpha_mid(NT alpha); 

/*!
sets `alpha_max`. 
*/ 
void set_alpha_max(NT alpha); 

/// @} 

/// \name Access Functions 
/// @{

/*!
Returns true for Gabriel faces. 
*/ 
bool is_Gabriel() const ; 

/*!
Returns true for convex hull faces. 
*/ 
bool is_on_chull() const; 

/*!
Returns the `alpha_min`. 
\pre `is_Gabriel()` returns false; 
*/ 
NT alpha_min() const; 

/*!
Returns the `alpha_mid`.
*/ 
NT alpha_mid() const; 

/*!
Returns `alpha_max`. 
\pre `is_on_chull()` returns `false`. 
*/ 
NT alpha_max() const; 

/// @}

}; /* end Alpha_status */
} /* end namespace CGAL */
