// List of named parameters used in the Polygon Mesh Processing package
CGAL_add_named_parameter(geom_traits_t, geom_traits, geom_traits)

// approximation parameters
CGAL_add_named_parameter(verbose_level_t, verbose_level, verbose_level)
CGAL_add_named_parameter(seeding_method_t, seeding_method, seeding_method)
CGAL_add_named_parameter(max_nb_proxies_t, max_nb_proxies, max_nb_proxies)
CGAL_add_named_parameter(min_error_drop_t, min_error_drop, min_error_drop)
CGAL_add_named_parameter(nb_of_iterations_t, nb_of_iterations, nb_of_iterations)
CGAL_add_named_parameter(nb_of_relaxations_t, nb_of_relaxations, nb_of_relaxations)

// meshing parameters
CGAL_add_named_parameter(mesh_chord_error_t, mesh_chord_error, mesh_chord_error)
CGAL_add_named_parameter(is_relative_to_chord_t, is_relative_to_chord, is_relative_to_chord)
CGAL_add_named_parameter(with_dihedral_angle_t, with_dihedral_angle, with_dihedral_angle)
CGAL_add_named_parameter(optimize_anchor_location_t, optimize_anchor_location, optimize_anchor_location)
CGAL_add_named_parameter(pca_plane_t, pca_plane, pca_plane)

// output parameters
CGAL_add_named_parameter(facet_proxy_map_t, facet_proxy_map, facet_proxy_map)
CGAL_add_named_parameter(proxies_t, proxies, proxies)
CGAL_add_named_parameter(anchors_t, anchors, anchors)
CGAL_add_named_parameter(triangles_t, triangles, triangles)
