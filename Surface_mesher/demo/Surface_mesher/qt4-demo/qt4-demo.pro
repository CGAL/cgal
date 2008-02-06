TEMPLATE = app
DEPENDPATH += ../../../include/CGAL \
              ../../../../Data_structure_for_queries_3/include/CGAL \
              ../../../include/CGAL/IO \
              ../../../include/CGAL/Surface_mesher \
              ../../../../Data_structure_for_queries_3/include/CGAL/Octree \
               ../../../../Marching_cube/include \
               ../../../../Marching_cube/src/mc

INCLUDEPATH += . \
               ../../../include \
               ../../../../Data_structure_for_queries_3/include \
               ../../../../Marching_cube/include

unix:LIBS += -lQGLViewer 
CONFIG += qt cgal
QT += xml opengl script


# Input
HEADERS += get_polyhedral_surface.h \
           polyhedral_surface.h \
           volume.h \
           surface.h \
           viewer.h \
           mainwindow.h \
           isovalues_list.h \
           isovalues_list.cpp \
           colorlisteditor.h \
           ../../../include/CGAL/Complex_2_in_triangulation_3.h \
           ../../../include/CGAL/Complex_2_in_triangulation_cell_base_3.h \
           ../../../include/CGAL/Complex_2_in_triangulation_vertex_base_3.h \
           ../../../include/CGAL/enriched_polyhedron.h \
           ../../../include/CGAL/Gray_level_image_3.h \
           ../../../include/CGAL/Implicit_surface_3.h \
           ../../../include/CGAL/make_piecewise_smooth_surface_mesh.h \
           ../../../include/CGAL/make_surface_mesh.h \
           ../../../include/CGAL/Multi_surface_3.h \
           ../../../include/CGAL/Point_traits.h \
           ../../../include/CGAL/Point_with_surface_index.h \
           ../../../include/CGAL/Point_with_surface_index_geom_traits.h \
           ../../../include/CGAL/Polyhedral_surface_3.h \
           ../../../include/CGAL/pws_loop_subdivision.h \
           ../../../include/CGAL/Robust_circumcenter_traits_3.h \
           ../../../include/CGAL/Surface_mesh_cell_base_3.h \
           ../../../include/CGAL/Surface_mesh_complex_2_in_triangulation_3.h \
           ../../../include/CGAL/Surface_mesh_default_criteria_3.h \
           ../../../include/CGAL/Surface_mesh_default_edges_criteria_3.h \
           ../../../include/CGAL/Surface_mesh_default_triangulation_3.h \
           ../../../include/CGAL/Surface_mesh_traits_generator_3.h \
           ../../../include/CGAL/Surface_mesh_triangulation_generator_3.h \
           ../../../include/CGAL/Surface_mesh_vertex_base_3.h \
           ../../../include/CGAL/Weighted_point_with_surface_index.h \
           ../../../include/CGAL/Weighted_point_with_surface_index_geom_traits.h \
           /home/lrineau/CGAL/Packages/Data_structure_for_queries_3/include/CGAL/Constrained_Edge.h \
           /home/lrineau/CGAL/Packages/Data_structure_for_queries_3/include/CGAL/Constrained_Element.h \
           /home/lrineau/CGAL/Packages/Data_structure_for_queries_3/include/CGAL/Constrained_Facet.h \
           /home/lrineau/CGAL/Packages/Data_structure_for_queries_3/include/CGAL/Constrained_Vertex.h \
           /home/lrineau/CGAL/Packages/Data_structure_for_queries_3/include/CGAL/Data_structure_using_Delaunay_triangulation_3.h \
           /home/lrineau/CGAL/Packages/Data_structure_for_queries_3/include/CGAL/Data_structure_using_octree_3.h \
           /home/lrineau/CGAL/Packages/Data_structure_for_queries_3/include/CGAL/DS_Container.h \
           /home/lrineau/CGAL/Packages/Data_structure_for_queries_3/include/CGAL/Triangulation_2_traits_3.h \
           ../../../include/CGAL/IO/Complex_2_in_triangulation_3_file_writer.h \
           ../../../include/CGAL/Surface_mesher/Combining_oracle.h \
           ../../../include/CGAL/Surface_mesher/Has_edges.h \
           ../../../include/CGAL/Surface_mesher/Implicit_surface_oracle_3.h \
           ../../../include/CGAL/Surface_mesher/Intersection_data_structure_3.h \
           ../../../include/CGAL/Surface_mesher/Null_oracle_visitor.h \
           ../../../include/CGAL/Surface_mesher/Point_surface_indices_oracle_visitor.h \
           ../../../include/CGAL/Surface_mesher/Polyhedral_oracle.h \
           ../../../include/CGAL/Surface_mesher/Simple_map_container.h \
           ../../../include/CGAL/Surface_mesher/Sphere_oracle_3.h \
           ../../../include/CGAL/Surface_mesher/Standard_criteria.h \
           ../../../include/CGAL/Surface_mesher/Surface_mesher.h \
           ../../../include/CGAL/Surface_mesher/Surface_mesher_edges_level.h \
           ../../../include/CGAL/Surface_mesher/Surface_mesher_edges_level_visitor.h \
           ../../../include/CGAL/Surface_mesher/Surface_mesher_manifold.h \
           ../../../include/CGAL/Surface_mesher/Surface_mesher_regular_edges.h \
           ../../../include/CGAL/Surface_mesher/Surface_mesher_visitor.h \
           ../../../include/CGAL/Surface_mesher/Types_generators.h \
           ../../../include/CGAL/Surface_mesher/Verbose_flag.h \
           ../../../include/CGAL/Surface_mesher/Vertices_on_the_same_surface_criterion.h \
           /home/lrineau/CGAL/Packages/Data_structure_for_queries_3/include/CGAL/Octree/Octree_helping_traits_3.h
INTERFACES += ui/mainwindow.ui ui/optionsdialog.ui ui/meshing_bar.ui ui/isovalues_list.ui
SOURCES += surface_mesher.cpp \
           viewer.cpp \
           mainwindow.cpp \
           polyhedral_surface.cpp \
           volume.cpp \
           isovalues_list.cpp \
           colorlisteditor.cpp \
           ../../../../Marching_cube/src/mc/MarchingCubes.cpp \
           ../../../../Marching_cube/src/mc/ply.c
RESOURCES += surface_mesher.qrc isovalues_list.qrc
