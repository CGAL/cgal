namespace CGAL {

/*!
\ingroup PkgSkinSurface3

constructs a mesh of the skin surface defined
by the weighted points and the shrink factor. 

The function `make_skin_surface_mesh_3()` constructs a mesh isotopic to the skin 
surface based on the algorithm in \cgalCite{cgal:kv-mssct-05}. It takes 
as input a range of weighted points and a shrink factor and outputs 
the mesh in a `Polyhedron_3` object. A number of subdivision 
steps might be applied to refine the mesh. 

\tparam WP_iterator must be an input iterator with weighted points as value type.
\tparam Polyhedron must be an instance of `Polyhedron_3`.

\pre `Polyhedron::HDS` can be used as the template argument of 
the `Polyhedron_incremental_builder_3<HDS>`.
*/

template <class WP_iterator,
	  class Polyhedron>
void make_skin_surface_mesh_3
(Polyhedron &p, WP_iterator begin, WP_iterator end, double
shrink_factor=.5, int nSubdivisions=0, bool
grow_balls=true);

} /* namespace CGAL */

