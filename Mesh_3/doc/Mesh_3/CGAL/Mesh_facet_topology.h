namespace CGAL {

/*!
\ingroup PkgMesh_3Enum

The enum `Mesh_facet_topology` is designed to tell which constraints have to
be checked on each surface facet during the mesh refinement process.

\sa `CGAL::Mesh_criteria_3<Tr>`, 
\sa `CGAL::Mesh_facet_criteria_3<Tr>`.
*/
enum Mesh_facet_topology {
  FACET_VERTICES_ON_SURFACE = 1, //!< Each vertex of the facet have 
                                 //!< to be on the surface, on a curve segment, or on a corner.
  FACET_VERTICES_ON_SAME_SURFACE_PATCH, //!< The three vertices of a facet belonging 
                                        //!< to a surface patch `s` have to be on 
                                        //!< the same surface patch `s`, on a curve segment or on a corner.
  /*!
    The three vertices of a facet belonging to a surface patch `s`
    have to be on the same surface patch `s`, or on a curve segment
    incident to the surface patch `s` or on a corner incident to the
    surface patch `s`.
  */
  FACET_VERTICES_ON_SAME_SURFACE_PATCH_WITH_ADJACENCY_CHECK  

};

}
