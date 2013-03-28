#ifndef SCENE_EDIT_POLYHEDRON_ITEM_H
#define SCENE_EDIT_POLYHEDRON_ITEM_H

#include "Scene_edit_polyhedron_item_config.h"
#include "Scene_polyhedron_item.h"
#include "Polyhedron_type.h"
#include <CGAL/boost/graph/halfedge_graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/properties.h>

#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <queue>

#include <QColor>
#include <QList>

#include <QGLViewer/manipulatedFrame.h>
#include <QGLViewer/qglviewer.h>

#include <CGAL/glu.h>

#include "ui_Deform_mesh.h"

#include "Property_maps_for_edit_plugin.h"
#include <CGAL/Deform_mesh.h> 

typedef Polyhedron_vertex_deformation_index_map<Polyhedron> Vertex_index_map;
typedef Polyhedron_edge_deformation_index_map<Polyhedron> Edge_index_map;

#if defined(CGAL_EIGEN3_ENABLED)
#include <CGAL/Eigen_solver_traits.h>
  #if defined(CGAL_SUPERLU_ENABLED)
    #include <Eigen/SuperLUSupport>
    typedef CGAL::Eigen_solver_traits<Eigen::SuperLU<CGAL::Eigen_sparse_matrix<double>::EigenType> > DefaultSolver;
  #else
    #include <Eigen/SparseLU>
    typedef CGAL::Eigen_solver_traits<
                Eigen::SparseLU<
                  CGAL::Eigen_sparse_matrix<double, Eigen::ColMajor>::EigenType,
                  Eigen::COLAMDOrdering<int> >  > DefaultSolver;
  #endif
#elif defined(CGAL_TAUCS_ENABLED)
  #include <CGAL/Taucs_solver_traits.h>
  typedef CGAL::Taucs_solver_traits<double> DefaultSolver;
#else
  typedef CGAL::Eigen_solver_traits<Eigen::BiCGSTAB<CGAL::Eigen_sparse_matrix<double>::EigenType> > DefaultSolver;
#endif

typedef CGAL::Deform_mesh<Polyhedron, DefaultSolver, Vertex_index_map, Edge_index_map, CGAL::ORIGINAL_ARAP> Deform_mesh;

typedef boost::graph_traits<Polyhedron>::vertex_descriptor		vertex_descriptor;
typedef boost::graph_traits<Polyhedron>::vertex_iterator		  vertex_iterator;
typedef boost::graph_traits<Polyhedron>::edge_descriptor		  edge_descriptor;
typedef boost::graph_traits<Polyhedron>::edge_iterator		    edge_iterator;
typedef boost::graph_traits<Polyhedron>::in_edge_iterator		  in_edge_iterator;
typedef boost::graph_traits<Polyhedron>::out_edge_iterator		  out_edge_iterator;

typedef Deform_mesh::Point  Point;
typedef Deform_mesh::Vector Vector;

/// For storing associated data with handle group
struct Handle_group_data
{
  Deform_mesh::Handle_group handle_group;
  qglviewer::ManipulatedFrame* frame;
  qglviewer::Vec frame_initial_center; // initial center of frame
  Scene_interface::Bbox bbox;

  Handle_group_data(Deform_mesh::Handle_group handle_group, qglviewer::ManipulatedFrame* frame = 0)
    : handle_group(handle_group), frame(frame), bbox(0,0,0,0,0,0)
  { }
};

// To hold pressing states together
struct Mouse_keyboard_state
{
  bool ctrl_pressing;
  bool shift_pressing;
  bool left_button_pressing;
  bool right_button_pressing;

  Mouse_keyboard_state() 
    : ctrl_pressing(false), shift_pressing(false), left_button_pressing(false), right_button_pressing(false)
  { }
};

#include <QGLViewer/camera.h>
class ManipulatedFrameWithCustomGrabber 
  : public qglviewer::ManipulatedFrame 
{
public:
  ManipulatedFrameWithCustomGrabber(std::list<Handle_group_data>* all_frames, Mouse_keyboard_state* state) 
    : all_frames(all_frames), state(state)
  { }

  void checkIfGrabsMouse(int x, int y, const qglviewer::Camera* const camera)
  {
    if(state->ctrl_pressing && (state->left_button_pressing || state->right_button_pressing) ) 
    { // user is deforming currently don't change the state 
      return;  
    }

    if(! state->ctrl_pressing) 
    { // if ctrl is not pressed then deactivate all handle manipulators 

      // this part required because we are using our custom update function (when pressed/released to ctrl)
      QGLViewer* viewer = *QGLViewer::QGLViewerPool().begin(); 
      viewer->setMouseGrabber(NULL);  
      ////////////////////

      // this part is required because when user release ctrl but keep left/right button pressed,
      // previous state at frame is kept. 
      action_ = QGLViewer::NO_MOUSE_ACTION;
      stopSpinning();
      ////////////////////////////

      setGrabsMouse(false);    
      return;
    }
    
    // now find closest frame and grab mouse
    // all frames can't be empty (at least contains 'this')
    std::list<Handle_group_data>::iterator min_it = all_frames->begin();    
    const qglviewer::Vec& pos_it = camera->projectedCoordinatesOf(min_it->frame->position());
    float min_dist = std::pow(pos_it.x - x, 2) + std::pow(pos_it.y - y, 2);

    for(std::list<Handle_group_data>::iterator it = all_frames->begin(); it != all_frames->end(); ++it)
    {
      const qglviewer::Vec& pos_it = camera->projectedCoordinatesOf(it->frame->position());
      float dist = std::pow(pos_it.x - x, 2) + std::pow(pos_it.y - y, 2);
      if(dist < min_dist) {
        min_dist = dist;
        min_it = it;
      }
    }

    // This part required because we are using our custom update function (when pressed/released to ctrl)
    if(min_it->frame == this)
    {
      QGLViewer* viewer = *QGLViewer::QGLViewerPool().begin();
      viewer->setMouseGrabber(this);      
    }
    //This part required... ends//////////////////
    setGrabsMouse(min_it->frame == this);
  }
protected:
  std::list<Handle_group_data>* all_frames;
  Mouse_keyboard_state* state;
};


// This class represents a polyhedron in the OpenGL scene
class SCENE_EDIT_POLYHEDRON_ITEM_EXPORT Scene_edit_polyhedron_item 
  : public Scene_item {
  Q_OBJECT
public:  
  /// Create an Scene_edit_polyhedron_item from a Scene_polyhedron_item.
  /// The ownership of the polyhedron is moved to the new edit_polyhedron
  /// item.
  Scene_edit_polyhedron_item(Scene_polyhedron_item* poly_item);
  ~Scene_edit_polyhedron_item();

  /// Returns 0, so that one cannot clone an "edit polyhedron" item.
  Scene_edit_polyhedron_item* clone() const;

  // Function for displaying meta-data of the item
  QString toolTip() const;

  void setColor(QColor c);
  void setName(QString n);
  void setVisible(bool b);
  void setRenderingMode(RenderingMode m);
  
  // Indicate if rendering mode is supported
  bool supportsRenderingMode(RenderingMode m) const { return (m!=PointsPlusNormals); }
  // Points/Wireframe/Flat/Gouraud OpenGL drawing in a display list
  void draw() const;
  void draw_bbox(const Scene_interface::Bbox& bb ) const;
  void gl_draw_edge(double px, double py, double pz,
                          double qx, double qy, double qz) const;
  void gl_draw_point(const Point& p) const;

  bool manipulatable() const { return true; }
  qglviewer::ManipulatedFrame* manipulatedFrame();

  // Get wrapped polyhedron
  Polyhedron*       polyhedron();
  const Polyhedron* polyhedron() const;

  /// Returns a Scene_polyhedron_item from the edit polyhedron item, and
  /// transfer the ownership of the polyhedron to it.
  /// The item 'this' must be destroy just after a call to this function.
  Scene_polyhedron_item* to_polyhedron_item() const;

  // Get dimensions
  bool isFinite() const { return true; }
  bool isEmpty() const;
  Bbox bbox() const;

  // take mouse events from viewer, main-window does not work
  // take keyboard events from main-window, which is more stable
  bool eventFilter(QObject *target, QEvent *event);
  
protected:
  void timerEvent(QTimerEvent *event);

public slots:
  void changed();
  void select(double orig_x,
              double orig_y,
              double orig_z,
              double dir_x,
              double dir_y,
              double dir_z);
  
  void vertex_has_been_selected(void* void_ptr); // a vertex is selected by shift + left_click
  void deform();                                 // deform the mesh
  
signals:
  void mesh_deformed(Scene_edit_polyhedron_item* edit_item);      // emits when deformation is completed
  void mesh_repaint_needed(Scene_edit_polyhedron_item* edit_item);
public:

  Ui::DeformMesh* ui_widget;

  Scene_polyhedron_item* poly_item;
  qglviewer::ManipulatedFrame* frame;

  Deform_mesh deform_mesh;
  Deform_mesh::Handle_group active_group;
typedef std::list<Handle_group_data> Handle_group_data_list;
  Handle_group_data_list handle_frame_map; // keep list of handle_groups with assoc data

  bool show_roi; // draw roi points
  bool show_as_sphere;
  double length_of_axis; // for drawing axis at a handle group

  // by interleaving 'viewer's events (check constructor), keep followings:
  Mouse_keyboard_state state;

public:
  // Deformation related functions //
  void insert_handle(vertex_descriptor v)
  {
    if(!is_there_any_handle_group()) {
      print_message("There is no handle group, create one!");
      return; 
    } // no handle group to insert

    deform_mesh.insert_handle(active_group, v);

    refresh_manipulated_frame_center(active_group);
  }

  void insert_roi(vertex_descriptor v)
  {
    deform_mesh.insert_roi(v);
  }
  
  void erase_handle(vertex_descriptor v)
  {
    deform_mesh.erase_handle(v); // API should be safe enough to do that (without checking empty handle group etc.)

    if(!is_there_any_handle_group()) {
      print_message("There is no handle group!");
      return; 
    }
    
    refresh_manipulated_frame_center(active_group);
  }

  void erase_roi(vertex_descriptor v)
  {
    deform_mesh.erase_roi(v);
    erase_handle(v); // also erase from handles
  }

  void set_all_vertices_as_roi()
  {
    vertex_iterator vb, ve;
    for(boost::tie(vb, ve) = boost::vertices(*polyhedron()); vb != ve; ++vb)
    {
      insert_roi(*vb);
    }   
  }

  void clear_roi()
  {
    for(Handle_group_data_list::iterator it = handle_frame_map.begin(); it != handle_frame_map.end(); ++it)
    {
      delete it->frame;
    }
    handle_frame_map.clear();
    deform_mesh.clear();

    create_handle_group(); // create one new handle group
  } 

  void create_handle_group()
  {
    active_group = deform_mesh.create_handle_group();

    qglviewer::ManipulatedFrame* new_frame = 
      new ManipulatedFrameWithCustomGrabber(&handle_frame_map, &state);
    new_frame->setRotationSensitivity(2.0f);
    handle_frame_map.push_back(Handle_group_data(active_group, new_frame));
    refresh_manipulated_frame_center(active_group);

    connect(new_frame, SIGNAL(modified()), this, SLOT(deform()));  // OK we are deforming via timer,
    // but it makes demo more responsive if we also add this signal

    print_message("A new empty handle group is created.");
  }

  void delete_handle_group(bool create_new = true)
  {
    if(!is_there_any_handle_group()) { 
      print_message("There is no handle group to be deleted!");
      return; 
    } // no handle group

    // delete group representative    
    for(Handle_group_data_list::iterator it = handle_frame_map.begin(); it != handle_frame_map.end(); ++it)
    {
      if(it->handle_group == active_group)
      {
        delete it->frame;
        handle_frame_map.erase(it); 
        break;
      }
    }
    deform_mesh.erase_handle(active_group);

    // assign another handle_group to active_group
    Deform_mesh::Handle_group hgb, hge;
    if( is_there_any_handle_group(hgb, hge) )
    { 
      active_group = hgb; 
    } // no handle group 
    else if(create_new)
    { 
      create_handle_group(); 
    }
  }

  void prev_handle_group()
  {
    Deform_mesh::Handle_group hgb, hge;
    if( !is_there_any_handle_group(hgb, hge) ) {
      print_message("There is no handle group to iterate on!");
      return; 
    }
    // shift
    if(hgb == active_group) { active_group = --hge; }
    else                    {--active_group; }    
  }

  void next_handle_group()
  {
    Deform_mesh::Handle_group hgb, hge;
    if( !is_there_any_handle_group(hgb, hge) ) {
      print_message("There is no handle group to iterate on!");
      return; 
    }
    // shift
    if(--hge == active_group) { active_group = hgb; }
    else                      {++active_group; }    
  }

  void pivoting_end()
  {       
    for(Handle_group_data_list::iterator it = handle_frame_map.begin(); it != handle_frame_map.end(); ++it)
    {
      it->frame_initial_center = it->frame->position();
    }
    for(Handle_group_data_list::iterator it = handle_frame_map.begin(); it != handle_frame_map.end(); ++it)
    {
      it->frame->blockSignals(false);
    }
  }

  void pivoting_begin()
  {
    for(Handle_group_data_list::iterator it = handle_frame_map.begin(); it != handle_frame_map.end(); ++it)
    {
      it->frame->blockSignals(true);
    }
  }

protected:
  // Deformation related functions //
  void print_message(const QString& message)
  {
    // std::cout << message.toStdString() << std::endl;
    //if(ui_widget != NULL)
    //  ui_widget->MessageTextEdit->appendPlainText(message);
  }

  bool is_there_any_handle_group(Deform_mesh::Handle_group& hgb, Deform_mesh::Handle_group& hge)
  {
    boost::tie(hgb, hge) = deform_mesh.handle_groups();
    return hgb != hge;
  }

  bool is_there_any_handle_group()
  {
    Deform_mesh::Handle_group hgb, hge;
    return is_there_any_handle_group(hgb, hge);
  }

  bool is_there_any_handle()
  {
    Deform_mesh::Handle_group hgb, hge;
    if(!is_there_any_handle_group(hgb, hge)) { return false; }
    for(; hgb != hge; ++hgb)
    {
      Deform_mesh::Handle_iterator hb, he;
      boost::tie(hb, he) = deform_mesh.handles(hgb);
      if(hb != he) { return true; }      
    }
    return false;
  }

  void process_selection(vertex_descriptor v, int k_ring, bool is_roi, bool is_insert, bool use_euclidean)
  {
    std::cout << "Process k-ring: " << k_ring << " roi: " << is_roi << " insert: " << is_insert << std::endl;

    std::map<vertex_descriptor, int> neighs = use_euclidean ?
      extract_k_ring_with_distance(*polyhedron(), v, k_ring) :
      extract_k_ring(*polyhedron(), v, k_ring);

    for(std::map<vertex_descriptor, int>::iterator it = neighs.begin(); it != neighs.end(); ++it)
    {
      vertex_descriptor vh = it->first;
      if(is_roi) {
        if(is_insert) { insert_roi(vh); }
        else          { erase_roi(vh);  }
      }
      else {
        if(is_insert) { insert_handle(vh); }
        else          { erase_handle(vh);  }
      }
    }
  }

  void refresh_all_handle_centers()
  {
    Deform_mesh::Handle_group hgb, hge;
    for(boost::tie(hgb, hge) = deform_mesh.handle_groups(); hgb != hge; ++hgb)
    {
      refresh_manipulated_frame_center(hgb);
    }
  }

  void refresh_manipulated_frame_center(Deform_mesh::Handle_group hg)
  {
    Handle_group_data& hd = get_data(hg);

    hd.frame_initial_center = calculate_original_center(hg);
    hd.bbox = calculate_bbox(hg);

    qglviewer::ManipulatedFrame* hg_frame = hd.frame;
    hg_frame->blockSignals(true); // do not let it emit modified, which will cause a deformation
                                  // but we are just adjusting the center so it does not require a deformation
    hg_frame->setPosition(calculate_center(hg));
    hg_frame->blockSignals(false);
  }

  qglviewer::Vec calculate_center(Deform_mesh::Handle_group hg)
  {
    Point center_acc(0, 0, 0);
    Deform_mesh::Handle_iterator hb, he;
    std::size_t counter = 0;
    for(boost::tie(hb, he) = deform_mesh.handles(hg); hb != he; ++hb, ++counter)
    {
      center_acc = center_acc + ((*hb)->point() - CGAL::ORIGIN);
      // center_acc = center_acc + (deform_mesh.original_position(*hb) - CGAL::ORIGIN);
    }
    if(counter == 0) { return qglviewer::Vec(0,0,0); } 
    return qglviewer::Vec(center_acc.x() / counter, center_acc.y() / counter, center_acc.z() / counter);
  }

  qglviewer::Vec calculate_original_center(Deform_mesh::Handle_group hg)
  {
    Point center_acc(0, 0, 0);
    Deform_mesh::Handle_iterator hb, he;
    std::size_t counter = 0;
    for(boost::tie(hb, he) = deform_mesh.handles(hg); hb != he; ++hb, ++counter)
    {
      center_acc = center_acc + (deform_mesh.original_position(*hb) - CGAL::ORIGIN);
    }
    if(counter == 0) { return qglviewer::Vec(0,0,0); } 
    return qglviewer::Vec(center_acc.x() / counter, center_acc.y() / counter, center_acc.z() / counter);
  }

  Scene_interface::Bbox calculate_bbox(Deform_mesh::Handle_group hg)
  {    
    Deform_mesh::Handle_iterator hb, he;
    boost::tie(hb, he) = deform_mesh.handles(hg);
    if(hb == he) { return Scene_interface::Bbox(0,0,0,0,0,0); }
    const Deform_mesh::Point& p_i = deform_mesh.original_position(*hb);
    Scene_interface::Bbox bbox(p_i.x(), p_i.y(), p_i.z(),
                               p_i.x(), p_i.y(), p_i.z());
    for(; hb != he; ++hb)
    {
      const Deform_mesh::Point& p_i = deform_mesh.original_position(*hb);
      Scene_interface::Bbox bbox_it(p_i.x(), p_i.y(), p_i.z(),
                                    p_i.x(), p_i.y(), p_i.z());
      bbox = bbox + bbox_it;
    }
    return bbox;
  }

  
  Handle_group_data& get_data(Deform_mesh::Handle_group hg)
  {
    for(Handle_group_data_list::iterator it = handle_frame_map.begin(); it != handle_frame_map.end(); ++it)
    {
      if(it->handle_group == hg)
      {
        return *it;
      }
    }
    return *handle_frame_map.end(); // crash
  }

  Handle_group_data get_data(Deform_mesh::Const_handle_group hg) const
  {
    for(Handle_group_data_list::const_iterator it = handle_frame_map.begin(); it != handle_frame_map.end(); ++it)
    {
      if(it->handle_group == hg)
      {
        return *it;
      }
    }
    return *handle_frame_map.end(); // crash
  }

  std::map<vertex_descriptor, int> extract_k_ring(const Polyhedron &P, vertex_descriptor v, int k)
  {
    std::map<vertex_descriptor, int>  D;
    std::queue<vertex_descriptor>     Q;
    Q.push(v); D[v] = 0;

    int dist_v;
    while( !Q.empty() && (dist_v = D[Q.front()]) < k ) {
      v = Q.front();
      Q.pop();

      out_edge_iterator e, e_end;
      for(boost::tie(e, e_end) = boost::out_edges(v, P); e != e_end; e++)
      {
        vertex_descriptor new_v = boost::target(*e, P);
        if(D.insert(std::make_pair(new_v, dist_v + 1)).second) {
          Q.push(new_v);
        }
      }   
    }
    return D;
  }

  std::map<vertex_descriptor, int> extract_k_ring_with_distance(const Polyhedron &P, vertex_descriptor root, int k)
  {
    // tuple: <vertex, geodistance, euclidian distane>
    std::map<vertex_descriptor, int>  D;
    std::queue<vertex_descriptor>     Q;
    Q.push(root); D[root] = 0;

    int dist_v;
    double max_distance = 0.0;
    while( !Q.empty() && (dist_v = D[Q.front()]) < k ) {
      vertex_descriptor v = Q.front();
      Q.pop();

      out_edge_iterator e, e_end;
      for(boost::tie(e, e_end) = boost::out_edges(v, P); e != e_end; e++)
      {
        vertex_descriptor new_v = boost::target(*e, P);
        if(D.insert(std::make_pair(new_v, dist_v + 1)).second) {
          max_distance = (std::max)( (new_v->point() - root->point()).squared_length(), max_distance);
          Q.push(new_v);
        }
      }   
    }
    // now Q contains all nonprocessed 
    while( !Q.empty() ) {
      vertex_descriptor v = Q.front();
      Q.pop();

      out_edge_iterator e, e_end;
      for(boost::tie(e, e_end) = boost::out_edges(v, P); e != e_end; e++)
      {
        vertex_descriptor new_v = boost::target(*e, P);
        double distance = (new_v->point() - root->point()).squared_length();
        if(distance < max_distance)
        {
          if(D.insert(std::make_pair(new_v, dist_v + 1)).second) {
            Q.push(new_v);
          }
        }
      }   
    }

    return D;
  }

  void activate_closest_manipulated_frame(QPoint p)
  {
    QGLViewer* viewer = *QGLViewer::QGLViewerPool().begin();
    qglviewer::Camera* camera = viewer->camera();
    for(std::list<Handle_group_data>::iterator it = handle_frame_map.begin(); it != handle_frame_map.end(); ++it)
    {
      it->frame->checkIfGrabsMouse(p.x(), p.y(), camera);
    }
  }

protected:
  GLUquadric* quadric; // for drawing spheres

}; // end class Scene_edit_polyhedron_item

#endif // SCENE_EDIT_POLYHEDRON_ITEM_H
