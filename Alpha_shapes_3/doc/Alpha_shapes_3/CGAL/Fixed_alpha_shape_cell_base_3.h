
namespace CGAL {

/*!
\ingroup PkgAlphaShapes3

The class `Fixed_alpha_shape_cell_base_3` is the default model for the concept 
`FixedAlphaShapeCell_3`. 

The class has two parameters. The traits class `Traits` 
provides the number type for alpha values. 
The second parameter `Fb` is a base class instantiated by default 
with `CGAL::Triangulation_cell_base_3<Traits>`. 

\models ::FixedAlphaShapeCell_3 

*/
template< typename Traits, typename Fb >
class Fixed_alpha_shape_cell_base_3 : public Fb {
public:

/// @}

}; /* end Fixed_alpha_shape_cell_base_3 */
} /* end namespace CGAL */
