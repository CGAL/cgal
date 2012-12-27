
namespace CGAL {

/*!
\ingroup PkgAlphaShapes3

The class `Fixed_alpha_shape_cell_base_3` is the default model for the concept 
`FixedAlphaShapeCell_3`. 

\tparam Traits provides the number type for alpha values. 
\tparam Cb must be a cell base class and it is instantiated by default 
with `Triangulation_cell_base_3<Traits>`. 

\cgalModels `FixedAlphaShapeCell_3`

*/
template< typename Traits, typename Cb >
class Fixed_alpha_shape_cell_base_3 : public Cb {
public:

/// @}

}; /* end Fixed_alpha_shape_cell_base_3 */
} /* end namespace CGAL */
