#ifndef CGAL_SCENE_H
#define CGAL_SCENE_H

//#include <QtOpenGL/qgl.h>
//#include <QtCore/qglobal.h>
//#include <QMap>
//#include <QOpenGLVertexArrayObject>
//#include <QOpenGLBuffer>
//#include <QOpenGLShaderProgram>
//#include <CGAL/Qt/manipulatedFrame.h>
//#include <CGAL/Qt/qglviewer.h>

// local
#include "Viewer.h"
#include "types.h"

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
  //QGLContext *context;
  Bbox bbox() { return m_bbox; }
  ManipulatedFrame* manipulatedFrame() const { return m_frame; }
  void initGL();
  void draw(CGAL::QGLViewer *viewer);

  // 2) file process
  bool open(QString file_name);
  bool open_input(QString file_name);
  bool open_remesh(QString file_name);
  void save_remesh_as(QString file_name);

  // 3) parameter settings access
  // isotropic remeshing parameters
  double get_target_edge_length() const { return m_target_edge_length; }

  void set_target_edge_length(double value) { m_target_edge_length = value; }

  int get_smooth_iteration_count() const { return m_smooth_iteration_count; }

  void set_smooth_iteration_count(int value) { m_smooth_iteration_count = value; }

  double get_max_error_threshold() const {
    return m_minangle_remesh.get_remesher()->get_max_error_threshold();
  }

  // min angle remeshing parameters
  void set_max_error_threshold(double value) {
    m_minangle_remesh.get_remesher()->set_max_error_threshold(value);
  }

  double get_min_angle_threshold() const {
    return m_minangle_remesh.get_remesher()->get_min_angle_threshold();
  }

  void set_min_angle_threshold(double value) {
    m_minangle_remesh.get_remesher()->set_min_angle_threshold(value);
  }

  int get_max_mesh_complexity() const {
    return m_minangle_remesh.get_remesher()->get_max_mesh_complexity();
  }

  void set_max_mesh_complexity(int value) {
    m_minangle_remesh.get_remesher()->set_max_mesh_complexity(value);
  }

  double get_smooth_angle_delta() const {
    return m_minangle_remesh.get_remesher()->get_smooth_angle_delta();
  }

  void set_smooth_angle_delta(double value) {
    m_minangle_remesh.get_remesher()->set_smooth_angle_delta(value);
  }

  bool get_apply_edge_flip() const {
    return m_minangle_remesh.get_remesher()->get_apply_edge_flip();
  }

  void set_apply_edge_flip(bool value) {
    m_minangle_remesh.get_remesher()->set_apply_edge_flip(value);
  }

  EdgeFlipStrategy get_edge_flip_strategy() const {
    return m_minangle_remesh.get_remesher()->get_edge_flip_strategy();
  }

  void set_edge_flip_strategy(EdgeFlipStrategy value) {
    m_minangle_remesh.get_remesher()->set_edge_flip_strategy(value);
  }

  bool get_flip_after_split_and_collapse() const {
    return m_minangle_remesh.get_remesher()->get_flip_after_split_and_collapse();
  }

  void set_flip_after_split_and_collapse(bool value) {
    m_minangle_remesh.get_remesher()->set_flip_after_split_and_collapse(value);
  }

  bool get_relocate_after_local_operations() const {
    return m_minangle_remesh.get_remesher()->get_relocate_after_local_operations();
  }

  void set_relocate_after_local_operations(bool value) {
    m_minangle_remesh.get_remesher()->set_relocate_after_local_operations(value);
  }

  RelocateStrategy get_relocate_strategy() const {
    return m_minangle_remesh.get_remesher()->get_relocate_strategy();
  }

  void set_relocate_strategy(RelocateStrategy value) {
    m_minangle_remesh.get_remesher()->set_relocate_strategy(value);
  }

  bool get_keep_vertex_in_one_ring() const {
    return m_minangle_remesh.get_remesher()->get_keep_vertex_in_one_ring();
  }

  void set_keep_vertex_in_one_ring(bool value) {
    m_minangle_remesh.get_remesher()->set_keep_vertex_in_one_ring(value);
  }

  bool get_use_local_aabb_tree() const {
    return m_minangle_remesh.get_remesher()->get_use_local_aabb_tree();
  }

  void set_use_local_aabb_tree(bool value) {
    m_minangle_remesh.get_remesher()->set_use_local_aabb_tree(value);
  }

  int get_collapsed_list_size() const {
    return m_minangle_remesh.get_remesher()->get_collapsed_list_size();
  }

  void set_collapsed_list_size(int value) {
    m_minangle_remesh.get_remesher()->set_collapsed_list_size(value);
  }

  bool get_decrease_max_errors() const {
    return m_minangle_remesh.get_remesher()->get_decrease_max_errors();
  }

  void set_decrease_max_errors(bool value) {
    m_minangle_remesh.get_remesher()->set_decrease_max_errors(value);
  }

  bool get_verbose_progress() const {
    return m_minangle_remesh.get_remesher()->get_verbose_progress();
  }

  void set_verbose_progress(bool value) {
    m_minangle_remesh.get_remesher()->set_verbose_progress(value);
  }

  bool get_apply_initial_mesh_simplification() const {
    return m_minangle_remesh.get_remesher()->get_apply_initial_mesh_simplification();
  }

  void set_apply_initial_mesh_simplification(bool value) {
    m_minangle_remesh.get_remesher()->set_apply_initial_mesh_simplification(value);
  }

  bool get_apply_final_vertex_relocation() const {
    return m_minangle_remesh.get_remesher()->get_apply_final_vertex_relocation();
  }

  void set_apply_final_vertex_relocation(bool value) {
    m_minangle_remesh.get_remesher()->set_apply_final_vertex_relocation(value);
  }
  
  int get_samples_per_face_in() const {
    return m_minangle_remesh.get_remesher()->get_samples_per_face_in();
  }

  void set_samples_per_face_in(int value) {
    m_minangle_remesh.get_remesher()->set_samples_per_face_in(value);
  }

  int get_samples_per_face_out() const {
    return m_minangle_remesh.get_remesher()->get_samples_per_face_out();
  }

  void set_samples_per_face_out(int value) {
    m_minangle_remesh.get_remesher()->set_samples_per_face_out(value);
  }

  int get_max_samples_per_area() const {
    return m_minangle_remesh.get_remesher()->get_max_samples_per_area();
  }

  void set_max_samples_per_area(int value) {
    m_minangle_remesh.get_remesher()->set_max_samples_per_area(value);
  }

  int get_min_samples_per_triangle() const {
    return m_minangle_remesh.get_remesher()->get_min_samples_per_triangle();
  }

  void set_min_samples_per_triangle(int value) {
    m_minangle_remesh.get_remesher()->set_min_samples_per_triangle(value);
  }

  int get_bvd_iteration_count() const {
    return m_minangle_remesh.get_remesher()->get_bvd_iteration_count();
  }

  void set_bvd_iteration_count(int value) {
    m_minangle_remesh.get_remesher()->set_bvd_iteration_count(value);
  }

  SampleNumberStrategy get_sample_number_strategy() const {
    return m_minangle_remesh.get_remesher()->get_sample_number_strategy();
  }

  void set_sample_number_strategy(SampleNumberStrategy value) {
    m_minangle_remesh.get_remesher()->set_sample_number_strategy(value);
  }

  SampleStrategy get_sample_strategy() const {
    return m_minangle_remesh.get_remesher()->get_sample_strategy();
  }

  void set_sample_strategy(SampleStrategy value) {
    m_minangle_remesh.get_remesher()->set_sample_strategy(value);
  }

  bool get_use_stratified_sampling() const {
    return m_minangle_remesh.get_remesher()->get_use_stratified_sampling();
  }

  void set_use_stratified_sampling(bool value) {
    m_minangle_remesh.get_remesher()->set_use_stratified_sampling(value);
  }

  double get_sum_theta() const { 
    return m_minangle_remesh.get_remesher()->get_sum_theta();
  }

  void set_sum_theta(double value) { 
    m_minangle_remesh.get_remesher()->set_sum_theta(value);
  }

  double get_sum_delta() const {
    return m_minangle_remesh.get_remesher()->get_sum_delta();
  }

  void set_sum_delta(double value) {
    m_minangle_remesh.get_remesher()->set_sum_delta(value);
  }

  double get_dihedral_theta() const {
    return m_minangle_remesh.get_remesher()->get_dihedral_theta();
  }

  void set_dihedral_theta(double value) {
    m_minangle_remesh.get_remesher()->set_dihedral_theta(value);
  }

  double get_dihedral_delta() const {
    return m_minangle_remesh.get_remesher()->get_dihedral_delta();
  }

  void set_dihedral_delta(double value) {
    m_minangle_remesh.get_remesher()->set_dihedral_delta(value);
  }

  double get_feature_difference_delta() const {
    return m_minangle_remesh.get_remesher()->get_feature_difference_delta();
  }

  void set_feature_difference_delta(double value) {
    m_minangle_remesh.get_remesher()->set_feature_difference_delta(value);
  }

  double get_feature_control_delta() const {
    return m_minangle_remesh.get_remesher()->get_feature_control_delta();
  }

  void set_feature_control_delta(double value) {
    m_minangle_remesh.get_remesher()->set_feature_control_delta(value);
  }

  bool get_inherit_element_types() const {
    return m_minangle_remesh.get_remesher()->get_inherit_element_types();
  }

  void set_inherit_element_types(bool value) {
    m_minangle_remesh.get_remesher()->set_inherit_element_types(value);
  }

  bool get_use_feature_intensity_weights() const {
    return m_minangle_remesh.get_remesher()->get_use_feature_intensity_weights();
  }

  void set_use_feature_intensity_weights(bool value) {
    m_minangle_remesh.get_remesher()->set_use_feature_intensity_weights(value);
  }

  int get_vertex_optimize_count() const {
    return m_minangle_remesh.get_remesher()->get_vertex_optimize_count();
  }

  void set_vertex_optimize_count(int value) {
    m_minangle_remesh.get_remesher()->set_vertex_optimize_count(value);
  }

  double get_vertex_optimize_ratio() const {
    return m_minangle_remesh.get_remesher()->get_vertex_optimize_ratio();
  }

  void set_vertex_optimize_ratio(double value) {
    m_minangle_remesh.get_remesher()->set_vertex_optimize_ratio(value);
  }

  int get_stencil_ring_size() const {
    return m_minangle_remesh.get_remesher()->get_stencil_ring_size();
  }

  void set_stencil_ring_size(int value) {
    m_minangle_remesh.get_remesher()->set_stencil_ring_size(value);
  }

  OptimizeStrategy get_optimize_strategy() const {
    return m_minangle_remesh.get_remesher()->get_optimize_strategy();
  }

  void set_optimize_strategy(OptimizeStrategy value) {
    m_minangle_remesh.get_remesher()->set_optimize_strategy(value);
  }

  OptimizeType get_face_optimize_type() const {
    return m_minangle_remesh.get_remesher()->get_face_optimize_type();
  }

  void set_face_optimize_type(OptimizeType value) {
    m_minangle_remesh.get_remesher()->set_face_optimize_type(value);
  }

  OptimizeType get_edge_optimize_type() const {
    return m_minangle_remesh.get_remesher()->get_edge_optimize_type();
  }

  void set_edge_optimize_type(OptimizeType value) {
    m_minangle_remesh.get_remesher()->set_edge_optimize_type(value);
  }

  OptimizeType get_vertex_optimize_type() const {
    return m_minangle_remesh.get_remesher()->get_vertex_optimize_type();
  }

  void set_vertex_optimize_type(OptimizeType value) {
    m_minangle_remesh.get_remesher()->set_vertex_optimize_type(value);
  }

  bool get_optimize_after_local_operations() const {
    return m_minangle_remesh.get_remesher()->get_optimize_after_local_operations();
  }

  void set_optimize_after_local_operations(bool value) {
    m_minangle_remesh.get_remesher()->set_optimize_after_local_operations(value);
  }

  // 4) toggle view option
  void toggle_view_input();
  void toggle_view_remesh();
  void toggle_view_input_remesh();
  void toggle_view_minimal_angle();              // the minimal angle
  void toggle_view_mesh_edges();
  void toggle_view_mesh_plain_faces();
  void toggle_view_face_errors();
  void toggle_view_interpolated_feature_intensities();
  void toggle_view_element_classifications();
  void toggle_view_gaussian_curvatures();
  void toggle_view_maximal_normal_dihedrals();
  void toggle_view_normal_dihedrals();
  void toggle_view_face_in_start_points();           // face in links
  void toggle_view_face_in_end_points();
  void toggle_view_face_in_links();
  void toggle_view_face_out_start_points();          // face out links
  void toggle_view_face_out_end_points();
  void toggle_view_face_out_links();
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
  void toggle_view_face_feature_intensities();       // face samples
  void toggle_view_face_capacities();
  void toggle_view_face_weights();

  // 5) menu functions
  void eliminate_degenerations();               // input menu
  void split_input_long_edges();
  void input_properties();
  void split_borders();                         // isotropic remeshing menu
  void isotropic_remeshing();
  void reset_from_input();                      // min angle remeshing menu
  void generate_links();
  void remesh_properties();
  void minangle_remeshing();                           
  void initial_mesh_simplification();
  void split_local_longest_edge();  
  void increase_minimal_angle();
  void maximize_minimal_angle();
  void final_vertex_relocation();
  void test();                                  // test the minimal_angle_remeshing

  // 6) operations (may need update)
  void update_feature_intensities();

  // 7) status access
  bool get_link_initialized() const {
    return m_minangle_remesh.get_remesher()->get_links_initialized();
  }

  // 8) utilizations
  int get_optimize_type_index(OptimizeType ot) const;

  OptimizeType get_optimize_type(int index) const;

 private:
   // 1) for rendering
  void compile_shaders();
  void attrib_buffers(CGAL::QGLViewer*);
  void initialize_buffers();                      // initialize buffers
  void changed();                                 // compute elements
  void set_draw_render_types(DrawType draw_type, RenderType render_type);
  void reset_draw_render_types();
  void compute_elements();
  void compute_classified_edges(bool is_input,
      std::vector<float> *pos_normal_edges,
      std::vector<float> *pos_special_edges) const;
  void compute_faces(bool is_input, Color face_color,
      std::vector<float> *pos_faces, std::vector<float> *pos_face_normals,
      std::vector<float> *pos_face_colors, std::vector<float> *pos_boundaries,
      std::vector<float> *pos_samples) const;
  void compute_mesh_faces(bool is_input, Color face_color, FT sum_theta_value,
      FT dihedral_theta_value, FT max_error_threshold_value,
      std::vector<float> *pos_faces, std::vector<float> *pos_face_normals,
      std::vector<float> *pos_face_colors) const;
  void compute_all_voronois(bool is_input, FT sum_theta_value,
      FT dihedral_theta_value, std::vector<float> *pos_faces,
      std::vector<float> *pos_face_normals, 
      std::vector<float> *pos_face_colors,std::vector<float> *pos_boundaries,
      std::vector<float> *pos_samples) const;
  void compute_vertex_voronois(bool is_input, FT sum_theta_value,
      FT dihedral_theta_value, std::vector<float> *pos_faces,
      std::vector<float> *pos_face_normals, 
      std::vector<float> *pos_face_colors, std::vector<float> *pos_boundaries,
      std::vector<float> *pos_samples) const;
  void compute_edge_voronois(bool is_input, FT sum_theta_value,
      FT dihedral_theta_value, std::vector<float> *pos_faces, 
      std::vector<float> *pos_face_normals, 
      std::vector<float> *pos_face_colors, std::vector<float> *pos_boundaries,
      std::vector<float> *pos_samples) const;
  void compute_face_voronois(bool is_input, FT sum_theta_value,
      FT dihedral_theta_value, std::vector<float> *pos_faces,
      std::vector<float> *pos_face_normals, 
      std::vector<float> *pos_face_colors, std::vector<float> *pos_boundaries,
      std::vector<float> *pos_samples) const;
  void compute_edge_normal_dihedrals(bool is_input, FT dihedral_theta_value, 
      std::vector<float> *pos_faces, std::vector<float> *pos_face_normals,
      std::vector<float> *pos_face_colors, 
      std::vector<float> *pos_boundaries) const;
  void compute_edge_sample_properties(bool is_input, FT sum_theta_value,
      FT dihedral_theta_value, std::vector<float> *pos_faces,
      std::vector<float> *pos_face_normals, 
      std::vector<float> *pos_face_colors, std::vector<float> *pos_boundaries,
      std::vector<float> *pos_samples) const;
  Color get_vertex_sample_normalized_color(bool is_input, vertex_descriptor vd,
      FT min_value, FT max_value, FT h) const;
  Color get_vertex_classification_color(bool is_input,
                                        vertex_descriptor vd) const;
  void get_all_sample_normalized_colors(bool is_input, face_descriptor fd,
      FT min_value, FT max_value, FT h, Point_list *samples, 
      Color_list *colors) const;
  void get_edge_sample_normalized_colors(bool is_input, halfedge_descriptor hd,
      FT min_value, FT max_value, FT h, Color_list *colors) const;
  void get_face_sample_normalized_colors(bool is_input, face_descriptor fd,
      FT min_value, FT max_value, FT h, Point_list *samples,
      Color_list *colors) const;
  void get_disturbed_border_samples_with_weights(bool is_input, 
      face_descriptor fd, FT disturb_ratio, 
      std::map<Point, double, Point_Comp> *disturbed_border_samples) const;
  void compute_all_samples(bool is_input,
                           std::vector<float> *pos_samples) const;
  void compute_face_samples(bool is_input,
                            std::vector<float> *pos_samples) const;
  void compute_vertices(bool is_input, std::vector<float> *pos_samples) const;
  void compute_edges(bool is_input, std::vector<float> *pos_edges) const;
  void compute_min_radian_edges(bool is_input,
      std::vector<float> *pos_min_radian_edges) const;
  void compute_halfedge(const Mesh &mesh, halfedge_descriptor hd,
                        std::vector<float> *pos) const;
  void compute_vertex(const Mesh &mesh, vertex_descriptor vd,
                      std::vector<float> *pos) const;
  void compute_triangle(const Point &p1, const Point &p2, const Point &p3,
      const Normal &normal, const Color &color, std::vector<float> *pos_faces,
      std::vector<float> *pos_face_normals,
      std::vector<float> *pos_face_colors) const;
  void compute_triangle_point(const Point &p, const Normal &normal,
      const Color &color, std::vector<float> *pos_faces,
      std::vector<float> *pos_face_normals,
      std::vector<float> *pos_face_colors) const;
  void compute_face_start_points(bool is_input,
    std::vector<float> *pos_face_start_point) const;
  void compute_face_end_points(bool is_input,
    std::vector<float> *pos_face_end_point) const;
  void compute_face_links(bool is_input,
    std::vector<float> *pos_face_links) const;
  void compute_edge_start_points(bool is_input,
    std::vector<float> *pos_edge_start_points) const;
  void compute_edge_end_points(bool is_input,
    std::vector<float> *pos_edge_end_points) const;
  void compute_edge_links(bool is_input,
    std::vector<float> *pos_edge_links) const;
  void compute_vertex_start_points(bool is_input,
    std::vector<float> *pos_vertex_start_points) const;
  void compute_vertex_end_points(bool is_input,
    std::vector<float> *pos_vertex_end_points) const;
  void compute_vertex_links(bool is_input,
    std::vector<float> *pos_vertex_links) const;
  void inline compute_segment(const Point &p, const Point &q,
                              std::vector<float> *pos) const;
  void inline compute_point(const Point &point,std::vector<float> *pos) const;
  void inline compute_normal(const Normal &normal, 
                             std::vector<float> *pos) const;
  void inline compute_color(const Color &color, std::vector<float> *pos) const;
  Color get_rgb_color(FT h, FT s, FT v) const;
  
  // 2) opeartions
  void reset();         // reset from input
  void update_bbox();
  void normalize(FT radius, Mesh *mesh) const;
  double calculate_input_edge_length() const;
  bool open_surface_mesh(QString file_name, Mesh *mesh) const;
  bool read_ply(std::ifstream &in, Mesh *mesh) const;
  inline const Mesh_properties* get_mesh_properties(bool is_input) const {
    const Minangle_remesher *remesher = m_minangle_remesh.get_remesher();
    return is_input ? remesher->get_input() : remesher->get_remesh();
  }

 private:
  // 1) general data
  // 1.1£© rendeing data
  ManipulatedFrame *m_frame;
  QOpenGLFunctions_2_1 *gl;
  bool gl_init;
  Bbox m_bbox;
  bool are_buffers_initialized;
  // 1.2) member data
  Mesh *m_pInput, *m_pRemesh;
  Minangle_remesh m_minangle_remesh;

  // 1.3) parameter settings (isotropic remeshing parameters)
  double m_target_edge_length;
  int m_smooth_iteration_count;

  // 2) Shaders elements
  enum VAOs {
    ka_Input_faces,               //       input faces or cells
    ka_Input_boundaries,          //       cell boundaries
    ka_Input_samples,             //       voronoi samples
    ka_Input_normal_edges,        //       edges
    ka_Input_special_edges,       //       minimal radian edges or crease edges
    ka_Remesh_faces,              //       remesh
    ka_Remesh_boundaries,
    ka_Remesh_samples,
    ka_Remesh_normal_edges,
    ka_Remesh_special_edges,
    ka_Face_in_start,             // face in links
    ka_Face_in_end,
    ka_Face_in_links,
    ka_Face_out_start,            // face out links
    ka_Face_out_end,
    ka_Face_out_links,
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
    kb_Input_face_pos,             // input
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
    kb_Face_in_start_points,       // face in links
    kb_Face_in_end_points,
    kb_Face_in_link_lines,
    kb_Face_out_start_points,      // face out links
    kb_Face_out_end_points,
    kb_Face_out_link_lines,
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
  bool m_view_input;
  bool m_view_remesh;
  bool m_view_face_in_start_points;              // face in links
  bool m_view_face_in_end_points;
  bool m_view_face_in_links;
  bool m_view_face_out_start_points;             // face out links
  bool m_view_face_out_end_points;
  bool m_view_face_out_links;
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
  bool m_view_mesh_edges;
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
  std::vector<float> pos_face_in_start_points;   // face in links
  std::vector<float> pos_face_in_end_points;
  std::vector<float> pos_face_in_links;
  std::vector<float> pos_face_out_start_points;  // face out links
  std::vector<float> pos_face_out_end_points;
  std::vector<float> pos_face_out_links;
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

  const double TARGET_EDGE_LENGTH = 0.2;
  const int SMOOTH_ITERATION_COUNT = 3;

}; // end class Scene

#endif // CGAL_SCENE_H
