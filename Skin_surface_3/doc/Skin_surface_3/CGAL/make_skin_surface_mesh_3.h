namespace CGAL {

/*!
\ingroup PkgSkinSurface3

Constructs a mesh of the skin surface defined
by the weighted points and the shrink factor. \pre `Polyhedron_3::HDS` can be used as the template argument of the `Polyhedron_incremental_builder_3<HDS>`.

The function `make_skin_surface_mesh_3` constructs a mesh isotopic to the skin 
surface based on the algorithm in \cite cgal:kv-mssct-05. It takes 
as input a range of weighted points and a shrink factor and outputs 
the mesh in a `Polyhedron_3` object. A number of subdivision 
steps might be applied to refine the mesh. 

*/
void mesh_skin_surface_3<Polyhedron_3>
(Polyhedron_3 &p, WP_iterator begin, WP_iterator end, double
shrink_factor=.5, int nSubdivisions=0, bool
grow_balls=true);

} /* namespace CGAL */

