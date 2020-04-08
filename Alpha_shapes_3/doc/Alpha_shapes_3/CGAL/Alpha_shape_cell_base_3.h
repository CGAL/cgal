
namespace CGAL {

/*!
\ingroup PkgAlphaShapes3Ref

The class `Alpha_shape_cell_base_3` is the default model for the concept
`AlphaShapeCell_3`.

\tparam Traits is the geometric traits class that is provided
to the `Alpha_shape_3` class.
\tparam Cb must be a cell base class adapted to the type of triangulation that is being used.
        By default, it is instantiated with `Delaunay_triangulation_cell_base_3<Traits>`,
        which is appropriate for basic alpha shapes.
\tparam ExactAlphaComparisonTag is a tag that, when set to
\link Tag_true `Tag_true`\endlink, triggers exact comparisons between alpha values. See the description
provided in the documentation of `Alpha_shape_3` for more details. The default value is \link Tag_false `Tag_false`\endlink.
\tparam WeightedTag is used only if `ExactAlphaComparisonTag` is \link Tag_true `Tag_true`\endlink. It
must be \link Tag_true `Tag_true`\endlink if the underlying triangulation of the alpha shape to be used is a regular triangulation
and \link Tag_false `Tag_false`\endlink otherwise. The default is \link Tag_false `Tag_false`\endlink.

\cgalModels `AlphaShapeCell_3`

\sa `Delaunay_triangulation_cell_base_3`
\sa `Regular_triangulation_cell_base_3`
\sa `Periodic_3_triangulation_ds_cell_base_3`
*/
template< typename Traits, typename Cb, typename ExactAlphaComparisonTag, typename WeightedTag >
class Alpha_shape_cell_base_3 : public Cb {
public:

}; /* end Alpha_shape_cell_base_3 */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgAlphaShapes3Ref

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
returns `true`  for Gabriel faces.
*/
bool is_Gabriel() const ;

/*!
returns `true` for convex hull faces.
*/
bool is_on_chull() const;

/*!
returns the `alpha_min`.
\pre `is_Gabriel()` returns `false`.
*/
NT alpha_min() const;

/*!
returns `alpha_mid`.
*/
NT alpha_mid() const;

/*!
returns `alpha_max`.
\pre `is_on_chull()` returns `false`.
*/
NT alpha_max() const;

/// @}

}; /* end Alpha_status */
} /* end namespace CGAL */
