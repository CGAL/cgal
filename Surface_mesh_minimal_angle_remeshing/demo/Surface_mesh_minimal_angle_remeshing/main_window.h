#ifndef CGAL_MAINWINDOW_H
#define CGAL_MAINWINDOW_H

// Qt
//#include <QtOpenGL/qgl.h>
#include <CGAL/Qt/DemosMainWindow.h>
// local
#include "Scene.h"

//class Viewer;
namespace Ui {
  class MainWindow;
}

class MainWindow : public CGAL::Qt::DemosMainWindow {
  Q_OBJECT
 public:

  enum OpenType {
    k_open_both = 0,
    k_open_input,
    k_open_remesh
  };

  // life cycle
  explicit MainWindow(QWidget *parent = 0);
  ~MainWindow();

  // disable copy/move construction
  MainWindow(const MainWindow &) = delete;
  MainWindow(const MainWindow &&) = delete;
  MainWindow &operator = (const MainWindow &) = delete;
  MainWindow &operator = (const MainWindow &&) = delete;

 public slots:
  void updateViewerBBox();
  void open(QString file_name);
  void setAddKeyFrameKeyboardModifiers(Qt::KeyboardModifiers);

 protected slots:

  // settings
  void quit();
  void readSettings();
  void writeSettings();

  // drag & drop
  void dropEvent(QDropEvent *event);
  void closeEvent(QCloseEvent *event);
  void dragEnterEvent(QDragEnterEvent *event);

  // file menu
  void on_actionFile_open_triggered();
  void on_actionFile_open_input_triggered();
  void on_actionFile_open_remesh_triggered();
  void on_actionFile_save_remesh_as_triggered();

  // edit menu
  void on_actionEdit_copy_snapshot_triggered();
  void on_actionEdit_save_snapshot_triggered();

  // input menu
  void on_actionInput_eliminate_degenerations_triggered();
  void on_actionInput_split_long_edges_triggered();
  void on_actionInput_properties_triggered();

  // isotropic remeshing menu
  void on_actionIsotropic_split_borders_triggered();
  void on_actionIsotropic_remeshing_triggered();
  void on_actionIsotropic_parameter_settings_triggered();

  //min angle remeshing menu
  void on_actionMinAngle_remesh_reset_from_input_triggered();
  void on_actionMinAngle_remesh_generate_links_triggered();
  void on_actionMinAngle_remeshing_triggered();
  void on_actionMinAngle_initial_mesh_simplification_triggered();
  void on_actionMinAngle_split_local_longest_edge_triggered();
  void on_actionMinAngle_increase_minimal_angle_triggered();
  void on_actionMinAngle_maximize_minimal_angle_triggered();
  void on_actionMinAngle_final_Vertex_relocation_triggered();
  void on_actionMinAngle_parameter_settings_triggered();
  void on_actionMinAngle_remesh_properties_triggered();
  void on_actionTest_triggered();
  
  // view menu
  void on_actionView_input_triggered();                           // input
  void on_actionView_remesh_triggered();                          // remesh
  void on_actionView_toggle_input_remesh_triggered();
  void on_actionView_mesh_edges_triggered();                      // surface mesh edges
  void on_actionView_minimal_angle_triggered();                   // surface mesh face properties
  void on_actionView_mesh_faces_plain_triggered();
  void on_actionView_mesh_faces_ifi_triggered();
  void on_actionView_mesh_faces_errors_triggered();
  void on_actionView_all_sample_feature_intensities_triggered();  // all samples properties
  void on_actionView_all_sample_capacities_triggered();
  void on_actionView_all_sample_weights_triggered();
  void on_actionView_element_classifications_triggered();          // vertex sample properties
  void on_actionView_gaussian_curvatures_triggered();
  void on_actionView_maximal_normal_dihedrals_triggered();
  void on_actionView_vertex_feature_intensities_triggered();
  void on_actionView_vertex_capacities_triggered();
  void on_actionView_vertex_weights_triggered();
  void on_actionView_normal_dihedrals_triggered();                // edge samples properties
  void on_actionView_edge_feature_intensities_triggered();        
  void on_actionView_edge_capacities_triggered();
  void on_actionView_edge_weights_triggered();
  void on_actionView_face_feature_intensities_triggered();        // face samples properties
  void on_actionView_face_capacities_triggered();
  void on_actionView_face_weights_triggered();
  void on_actionView_face_in_start_points_triggered();            // face in links
  void on_actionView_face_in_end_points_triggered();
  void on_actionView_face_in_links_triggered();
  void on_actionView_face_out_start_points_triggered();           // face out links
  void on_actionView_face_out_end_points_triggered();
  void on_actionView_face_out_links_triggered();
  void on_actionView_edge_in_start_points_triggered();            // edge in links
  void on_actionView_edge_in_end_points_triggered();
  void on_actionView_edge_in_links_triggered();
  void on_actionView_edge_out_start_points_triggered();           // edge out links
  void on_actionView_edge_out_end_points_triggered();
  void on_actionView_edge_out_links_triggered();
  void on_actionView_vertex_in_start_points_triggered();          // vertex in links
  void on_actionView_vertex_in_end_points_triggered();
  void on_actionView_vertex_in_links_triggered();
  void on_actionView_vertex_out_start_points_triggered();         // vertex out links
  void on_actionView_vertex_out_end_points_triggered();
  void on_actionView_vertex_out_links_triggered();

private:
  void open(QString file_name, OpenType open_type);
  void update_menu_items();

private:
  // objects
  Scene *m_pScene;
  Viewer *m_pViewer;
  Ui::MainWindow *ui;
};

#endif // CGAL_MAINWINDOW_H
