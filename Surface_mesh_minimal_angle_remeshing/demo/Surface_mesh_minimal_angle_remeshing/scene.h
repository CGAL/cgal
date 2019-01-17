#ifndef SCENE_H_
#define SCENE_H_

#include <QtOpenGL/qgl.h>
#include <QtCore/qglobal.h>
#include <QMap>
#include <QOpenGLVertexArrayObject>
#include <QOpenGLBuffer>
#include <QOpenGLShaderProgram>
#include <CGAL/Qt/manipulatedFrame.h>
#include <CGAL/Qt/qglviewer.h>

// local
#include "remesh.h"
#include "viewer.h"

class Scene : public QObject {
  Q_OBJECT

 public:

  // life cycle
  Scene();
  virtual ~Scene();

  // disable copy/move construction
  Scene(const Scene &) = delete;
  Scene(const Scene &&) = delete;
  Scene &operator = (const Scene &) = delete;
  Scene &operator = (const Scene &&) = delete;

 public:
  typedef CGAL::qglviewer::ManipulatedFrame ManipulatedFrame;

  // 1) for rendering
  QGLContext *context;      // TODO: do we really need this data?
  Bbox bbox() { return m_bbox; }
  ManipulatedFrame* manipulatedFrame() const { return m_frame; }
  void update_bbox();
  void initGL();
  void draw(CGAL::QGLViewer *viewer);

  // 2) file process
  bool open(QString file_name);
  bool open_input(QString file_name);
  bool open_remesh(QString file_name);
  void save_remesh_as(QString file_name);

  // 3) parameter settings access
  double get_max_error_threshold() const 
      { return remesh.get_max_error_threshold(); }
  void set_max_error_threshold(double value) 
      { remesh.set_max_error_threshold(value); }
  double get_min_angle_threshold() const 
      { return remesh.get_min_angle_threshold(); }
  void set_min_angle_threshold(double value) 
      { remesh.set_min_angle_threshold(value); }
  int get_max_mesh_complexity() const 
      { return remesh.get_max_mesh_complexity(); }
  void set_max_mesh_complexity(int value) 
      { remesh.set_max_mesh_complexity(value); }
  double get_smooth_angle_delta() const
      { return remesh.get_smooth_angle_delta(); }
  void set_smooth_angle_delta(double value) 
      { remesh.set_smooth_angle_delta(value); }
  bool get_apply_edge_flip() const 
      { return remesh.get_apply_edge_flip(); }
  void set_apply_edge_flip(bool value)
      { remesh.set_apply_edge_flip(value); }
  EdgeFlipStrategy get_edge_flip_strategy() const 
      { return remesh.get_edge_flip_strategy(); }
  void set_edge_flip_strategy(EdgeFlipStrategy value) 
      { remesh.set_edge_flip_strategy(value); }
  bool get_flip_after_split_and_collapse() const 
      { return remesh.get_flip_after_split_and_collapse(); }
  void set_flip_after_split_and_collapse(bool value) 
      { remesh.set_flip_after_split_and_collapse(value); }
  bool get_relocate_after_local_operations() const 
      { return remesh.get_relocate_after_local_operations(); }
  void set_relocate_after_local_operations(bool value)
      { remesh.set_relocate_after_local_operations(value); }
  RelocateStrategy get_relocate_strategy() const
      { return remesh.get_relocate_strategy(); }
  void set_relocate_strategy(RelocateStrategy value) 
      { remesh.set_relocate_strategy(value); }
  bool get_keep_vertex_in_one_ring() const
      { return remesh.get_keep_vertex_in_one_ring(); }
  void set_keep_vertex_in_one_ring(bool value) 
      { remesh.set_keep_vertex_in_one_ring(value); }
  bool get_use_local_aabb_tree() const
      { return remesh.get_use_local_aabb_tree(); }
  void set_use_local_aabb_tree(bool value)
      { remesh.set_use_local_aabb_tree(value); }
  int get_collapsed_list_size() const 
      { return remesh.get_collapsed_list_size(); }
  void set_collapsed_list_size(int value)
      { remesh.set_collapsed_list_size(value); }
  bool get_decrease_max_errors() const 
      { return remesh.get_decrease_max_errors(); }
  void set_decrease_max_errors(bool value) 
      { remesh.set_decrease_max_errors(value); }
  bool get_track_information() const
      { return remesh.get_track_information(); }
  void set_track_information(bool value) 
      { remesh.set_track_information(value); }
  bool get_apply_initial_mesh_simplification() const
      { return remesh.get_apply_initial_mesh_simplification(); }
  void set_apply_initial_mesh_simplification(bool value) 
      { remesh.set_apply_initial_mesh_simplification(value); }
  bool get_apply_final_vertex_relocation() const 
      { return remesh.get_apply_final_vertex_relocation(); }
  void set_apply_final_vertex_relocation(bool value) 
      { remesh.set_apply_final_vertex_relocation(value); }
  
  int get_samples_per_facet_in() const 
      { return remesh.get_samples_per_facet_in(); }
  void set_samples_per_facet_in(int value)
      { remesh.set_samples_per_facet_in(value); }
  int get_samples_per_facet_out() const 
      { return remesh.get_samples_per_facet_out(); }
  void set_samples_per_facet_out(int value) 
      { remesh.set_samples_per_facet_out(value); }
  int get_max_samples_per_area() const
      { return remesh.get_max_samples_per_area(); }
  void set_max_samples_per_area(int value)
      { remesh.set_max_samples_per_area(value); }
  int get_min_samples_per_triangle() const
      { return remesh.get_min_samples_per_triangle(); }
  void set_min_samples_per_triangle(int value)
      { remesh.set_min_samples_per_triangle(value); }
  int get_bvd_iteration_count() const
      { return remesh.get_bvd_iteration_count(); }
  void set_bvd_iteration_count(int value)
      { remesh.set_bvd_iteration_count(value); }
  SampleNumberStrategy get_sample_number_strategy() const
      { return remesh.get_sample_number_strategy(); }
  void set_sample_number_strategy(SampleNumberStrategy value) 
      { remesh.set_sample_number_strategy(value); }
  SampleStrategy get_sample_strategy() const
      { return remesh.get_sample_strategy(); }
  void set_sample_strategy(SampleStrategy value)
      { remesh.set_sample_strategy(value); }
  bool get_use_stratified_sampling() const 
      { return remesh.get_use_stratified_sampling(); }
  void set_use_stratified_sampling(bool value)
      { remesh.set_use_stratified_sampling(value); }

  double get_sum_theta() const { return remesh.get_sum_theta(); }
  void set_sum_theta(double value) { remesh.set_sum_theta(value); }
  double get_sum_delta() const { return remesh.get_sum_delta(); }
  void set_sum_delta(double value) { remesh.set_sum_delta(value); }
  double get_dihedral_theta() const { return remesh.get_dihedral_theta(); }
  void set_dihedral_theta(double value) { remesh.set_dihedral_theta(value); }
  double get_dihedral_delta() const { return remesh.get_dihedral_delta(); }
  void set_dihedral_delta(double value) { remesh.set_dihedral_delta(value); }
  double get_feature_difference_delta() const 
      { return remesh.get_feature_difference_delta(); }
  void set_feature_difference_delta(double value)
      { remesh.set_feature_difference_delta(value); }
  double get_feature_control_delta() const 
      { return remesh.get_feature_control_delta(); }
  void set_feature_control_delta(double value)
      { remesh.set_feature_control_delta(value); }
  bool get_inherit_element_types() const
      { return remesh.get_inherit_element_types(); }
  void set_inherit_element_types(bool value)
      { remesh.set_inherit_element_types(value); }
  bool get_use_feature_intensity_weights() const 
      { return remesh.get_use_feature_intensity_weights(); }
  void set_use_feature_intensity_weights(bool value)
      { remesh.set_use_feature_intensity_weights(value); }

  int get_vertex_optimize_count() const
      { return remesh.get_vertex_optimize_count(); }
  void set_vertex_optimize_count(int value) 
      { remesh.set_vertex_optimize_count(value); }
  double get_vertex_optimize_ratio() const 
      { return remesh.get_vertex_optimize_ratio(); }
  void set_vertex_optimize_ratio(double value) 
      { remesh.set_vertex_optimize_ratio(value); }
  int get_stencil_ring_size() const 
      { return remesh.get_stencil_ring_size(); }
  void set_stencil_ring_size(int value)
      { remesh.set_stencil_ring_size(value); }
  OptimizeStrategy get_optimize_strategy() const 
      { return remesh.get_optimize_strategy(); }
  void set_optimize_strategy(OptimizeStrategy value)
      { remesh.set_optimize_strategy(value); }
  OptimizeType get_facet_optimize_type() const
      { return remesh.get_facet_optimize_type(); }
  void set_facet_optimize_type(OptimizeType value) 
      { remesh.set_facet_optimize_type(value); }
  OptimizeType get_edge_optimize_type() const 
      { return remesh.get_edge_optimize_type(); }
  void set_edge_optimize_type(OptimizeType value) 
      { remesh.set_edge_optimize_type(value); }
  OptimizeType get_vertex_optimize_type() const
      { return remesh.get_vertex_optimize_type(); }
  void set_vertex_optimize_type(OptimizeType value)
      { remesh.set_vertex_optimize_type(value); }
  bool get_optimize_after_local_operations() const
      { return remesh.get_optimize_after_local_operations(); }
  void set_optimize_after_local_operations(bool value)
      { remesh.set_optimize_after_local_operations(value); }

  bool get_link_initialized() const { return m_links_initialized; }
  int get_optimize_type_index(OptimizeType ot) const;
  OptimizeType get_optimize_type(int index) const;

  // 4) toggle view option
  void toggle_view_input();
  void toggle_view_remesh();
  void toggle_view_input_remesh();
  void toggle_view_minimal_angle();              // the minimal angle
  void toggle_view_polyhedron_edges();
  void toggle_view_polyhedron_facets();
  void toggle_view_facet_errors();
  void toggle_view_interpolated_feature_intensities();
  void toggle_view_element_classifications();
  void toggle_view_gaussian_curvatures();
  void toggle_view_maximal_normal_dihedrals();
  void toggle_view_normal_dihedrals();
  void toggle_view_facet_in_start_points();           // facet in links
  void toggle_view_facet_in_end_points();
  void toggle_view_facet_in_links();
  void toggle_view_facet_out_start_points();          // facet out links
  void toggle_view_facet_out_end_points();
  void toggle_view_facet_out_links();
  void toggle_view_edge_in_start_points();            // edge in links
  void toggle_view_edge_in_end_points();
  void toggle_view_edge_in_links();
  void toggle_view_edge_out_start_points();           // edge out links
  void toggle_view_edge_out_end_points();
  void toggle_view_edge_out_links();
  void toggle_view_vertex_in_start_points();          // vertex in links
  void toggle_view_vertex_in_end_points();
  void toggle_view_vertex_in_links();
  void toggle_view_vertex_out_start_points();         // vertex out links
  void toggle_view_vertex_out_end_points();
  void toggle_view_vertex_out_links();
  void toggle_view_all_sample_feature_intensities();  // all samples
  void toggle_view_all_sample_capacities();
  void toggle_view_all_sample_weights();
  void toggle_view_vertex_feature_intensities();      // vertex samples
  void toggle_view_vertex_capacities();
  void toggle_view_vertex_weights();
  void toggle_view_edge_feature_intensities();        // edge samples
  void toggle_view_edge_capacities();
  void toggle_view_edge_weights();
  void toggle_view_facet_feature_intensities();       // facet samples
  void toggle_view_facet_capacities();
  void toggle_view_facet_weights();

  // 5) menu functions (may need update)
  void eliminate_degenerations();                       // input menu
  void split_input_long_edges();
  void input_properties();
  void reset_from_input();                              // remesh menu
  void generate_links_and_types();
  void remesh_properties();                        
  void isotropic_remeshing();                           // isotropic menu
  void initial_mesh_simplification();
  void split_local_longest_edge();  
  void increase_minimal_angle();
  void maximize_minimal_angle();
  void final_vertex_relocation();

  // 6) operations (may need update)
  void update_feature_intensities_and_clear_links();

 private:
   // 1) for renderring
  void compile_shaders();
  void attrib_buffers(CGAL::QGLViewer*);
  void initialize_buffers();                      // initialize buffers
  void changed();                                 // compute elements
  void compute_elements();
  void set_draw_render_types(DrawType draw_type, RenderType render_type);
  void reset_draw_render_types();
  // 2) opeartions
  void reset();         // reset from input
  void generate();      // generate links

 private:
  // 1) general data
  // 1.1£© rendeing data
  ManipulatedFrame *m_frame;
  QOpenGLFunctions_2_1 *gl;
  bool gl_init;
  Bbox m_bbox;
  bool are_buffers_initialized;
  // 1.2) member data
  Polyhedron *m_pInput, *m_pRemesh;
  Remesh<Kernel, Polyhedron> remesh;
  Facet_tree m_input_facet_tree, m_remesh_facet_tree;

  // 1) Shaders elements
  enum VAOs {
    ka_Input_faces = 0,           // input facets or cells
    ka_Input_boundaries,          //       cell boundaries
    ka_Input_samples,             //       samples
    ka_Input_normal_edges,        //       edges
    ka_Input_special_edges,       //       minimal radian edges or crease edges
    ka_Remesh_faces,              // remesh
    ka_Remesh_boundaries,
    ka_Remesh_samples,
    ka_Remesh_normal_edges,
    ka_Remesh_special_edges,
    ka_Facet_in_start,             // facet in links
    ka_Facet_in_end,
    ka_Facet_in_links,
    ka_Facet_out_start,            // facet out links
    ka_Facet_out_end,
    ka_Facet_out_links,
    ka_Edge_in_start,              // edge in links
    ka_Edge_in_end,
    ka_Edge_in_links,
    ka_Edge_out_start,             // edge out links
    ka_Edge_out_end,
    ka_Edge_out_links,
    ka_Vertex_in_start,            // vertex in links
    ka_Vertex_in_end,
    ka_Vertex_in_links,
    ka_Vertex_out_start,           // vertex out links
    ka_Vertex_out_end,
    ka_Vertex_out_links,
    ka_NbOfVaos
  };
  enum VBOs {
    kb_Input_face_pos = 0,         // input
    kb_Input_face_normals,
    kb_Input_face_colors,
    kb_Input_boundary_pos,
    kb_Input_sample_pos,
    kb_Input_normal_edge_pos,
    kb_Input_special_edge_pos,
    kb_Remesh_face_pos,            // remesh
    kb_Remesh_face_normals,
    kb_Remesh_face_colors,
    kb_Remesh_boundary_pos,
    kb_Remesh_sample_pos,
    kb_Remesh_normal_edge_pos,
    kb_Remesh_special_edge_pos,
    kb_Facet_in_start_points,      // facet in links
    kb_Facet_in_end_points,
    kb_Facet_in_link_lines,
    kb_Facet_out_start_points,     // facet out links
    kb_Facet_out_end_points,
    kb_Facet_out_link_lines,
    kb_Edge_in_start_points,       // edge in links
    kb_Edge_in_end_points,
    kb_Edge_in_link_lines,
    kb_Edge_out_start_points,      // edge out links
    kb_Edge_out_end_points,
    kb_Edge_out_link_lines,
    kb_Vertex_in_start_points,     // vertex in links
    kb_Vertex_in_end_points,
    kb_Vertex_in_link_lines,
    kb_Vertex_out_start_points,    // vertex out links
    kb_Vertex_out_end_points,
    kb_Vertex_out_link_lines,
    kb_NbOfVbos
  };
  QOpenGLShaderProgram rendering_program, rendering_program_with_light;
  QOpenGLVertexArrayObject vao[VAOs::ka_NbOfVaos];
  QOpenGLBuffer vbo[VBOs::kb_NbOfVbos];
  // 2.1) rendering program
  //int poly_vertexLocation;
  int points_vertexLocation;
  int lines_vertexLocation;
  int mvpLocation;
  int fLocation;
  int colorLocation;
  // 2.2) rendering program with light
  int poly_vertexLocation_with_light;
  int normalLocation_with_light;
  int mvpLocation_with_light;
  int fLocation_with_light;
  int colorLocation_with_light;
  
  // 3) rendering variables
  bool m_view_input;                              // input
  bool m_view_remesh;                             // remesh
  bool m_view_facet_in_start_points;              // facet in links
  bool m_view_facet_in_end_points;
  bool m_view_facet_in_links;
  bool m_view_facet_out_start_points;             // facet out links
  bool m_view_facet_out_end_points;
  bool m_view_facet_out_links;
  bool m_view_edge_in_start_points;               // edge in links
  bool m_view_edge_in_end_points;
  bool m_view_edge_in_links;
  bool m_view_edge_out_start_points;              // edge out links
  bool m_view_edge_out_end_points;
  bool m_view_edge_out_links;
  bool m_view_vertex_in_start_points;             // vertex in links
  bool m_view_vertex_in_end_points;
  bool m_view_vertex_in_links;
  bool m_view_vertex_out_start_points;            // vertex out links
  bool m_view_vertex_out_end_points;
  bool m_view_vertex_out_links;
  bool m_view_minimal_angle;
  bool m_view_polyhedron_edges;
  DrawType m_draw_type;                           // draw options
  RenderType m_render_type;
  std::vector<float> pos_input_faces;             // input
  std::vector<float> pos_input_face_normals;
  std::vector<float> pos_input_face_colors;
  std::vector<float> pos_input_boundaries;
  std::vector<float> pos_input_samples;
  std::vector<float> pos_input_normal_edges;
  std::vector<float> pos_input_special_edges;
  std::vector<float> pos_remesh_faces;            // remesh
  std::vector<float> pos_remesh_face_normals;
  std::vector<float> pos_remesh_face_colors;
  std::vector<float> pos_remesh_boundaries;
  std::vector<float> pos_remesh_samples;
  std::vector<float> pos_remesh_normal_edges;
  std::vector<float> pos_remesh_special_edges;
  std::vector<float> pos_facet_in_start_points;   // facet in links
  std::vector<float> pos_facet_in_end_points;
  std::vector<float> pos_facet_in_links;
  std::vector<float> pos_facet_out_start_points;  // facet out links
  std::vector<float> pos_facet_out_end_points;
  std::vector<float> pos_facet_out_links;
  std::vector<float> pos_edge_in_start_points;    // edge in links
  std::vector<float> pos_edge_in_end_points;
  std::vector<float> pos_edge_in_links;
  std::vector<float> pos_edge_out_start_points;   // edge out links
  std::vector<float> pos_edge_out_end_points;
  std::vector<float> pos_edge_out_links;
  std::vector<float> pos_vertex_in_start_points;  // vertx in links
  std::vector<float> pos_vertex_in_end_points;
  std::vector<float> pos_vertex_in_links;
  std::vector<float> pos_vertex_out_start_points; // vertex out links
  std::vector<float> pos_vertex_out_end_points;
  std::vector<float> pos_vertex_out_links;

  // 4) status data
  bool m_input_aabb_tree_constructed;
  bool m_links_initialized;

}; // end class Scene

#endif // SCENE_H_
