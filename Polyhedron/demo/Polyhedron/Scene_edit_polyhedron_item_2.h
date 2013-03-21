#ifndef SCENE_EDIT_POLYHEDRON_ITEM_2_H
#define SCENE_EDIT_POLYHEDRON_ITEM_2_H

#include "Scene_edit_polyhedron_item_config_2.h"
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
#include "Custom_manipulated_frame.h"

#include "ui_Deform_mesh_2.h"

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

typedef CGAL::Deform_mesh<Polyhedron, DefaultSolver, Vertex_index_map, Edge_index_map> Deform_mesh;

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
  qglviewer::Vec initial_center;

  Handle_group_data(Deform_mesh::Handle_group handle_group, qglviewer::ManipulatedFrame* frame = 0)
    : handle_group(handle_group), frame(frame)
  { }
};

// This class represents a polyhedron in the OpenGL scene
class SCENE_EDIT_POLYHEDRON_ITEM_2_EXPORT Scene_edit_polyhedron_item_2 
  : public Scene_item {
  Q_OBJECT
public:  
  /// Create an Scene_edit_polyhedron_item_2 from a Scene_polyhedron_item.
  /// The ownership of the polyhedron is moved to the new edit_polyhedron
  /// item.
  Scene_edit_polyhedron_item_2(Scene_polyhedron_item* poly_item);
  ~Scene_edit_polyhedron_item_2();

  /// Returns 0, so that one cannot clone an "edit polyhedron" item.
  Scene_edit_polyhedron_item_2* clone() const;

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

  // bool keyPressEvent(QKeyEvent*);

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
  void mesh_deformed(Scene_edit_polyhedron_item_2* edit_item); // emits when deformation is completed

public:
  Ui::DeformMesh_2* ui_widget;

  Scene_polyhedron_item* poly_item;
  qglviewer::ManipulatedFrame* frame;

  std::set<vertex_descriptor> roi;
  std::set<vertex_descriptor> handles;

  Deform_mesh deform_mesh;
  Deform_mesh::Handle_group active_group;
  typedef std::list<Handle_group_data> Handle_group_data_list;
  Handle_group_data_list handle_frame_map;

  bool show_roi; // draw roi points

public:
  // Deformation related functions //
  void insert_handle(vertex_descriptor v)
  {
    if(!is_there_any_handle_group()) {
      ui_widget->MessageTextEdit->appendPlainText("There is no handle group, create one!");
      return; 
    } // no handle group to insert

    if(handles.insert(v).second)
    {
      deform_mesh.insert_handle(active_group, v);
      refresh_manipulated_frame_center(active_group);
    }
    insert_roi(v); // also insert as roi
  }

  void insert_roi(vertex_descriptor v)
  {
    if(roi.insert(v).second)
    {
      deform_mesh.insert_roi(v);
    }
  }
  
  void erase_handle(vertex_descriptor v)
  {
    // this part is not totally correct - should add checking whether active_group contains v
    if(handles.erase(v) == 1)
    {
      deform_mesh.erase_handle(active_group, v);
      refresh_manipulated_frame_center(active_group);
    }
  }

  void erase_roi(vertex_descriptor v)
  {
    if(roi.erase(v) == 1)
    {
      deform_mesh.erase_roi(v);
    }
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

  void create_handle_group()
  {
    active_group = deform_mesh.create_handle_group();

    qglviewer::ManipulatedFrame* new_frame = new CustomManipulatedFrame(0);//new qglviewer::ManipulatedFrame();
    handle_frame_map.push_back(Handle_group_data(active_group, new_frame));
    refresh_manipulated_frame_center(active_group);

    connect(new_frame, SIGNAL(modified()), this, SLOT(deform()));  

    ui_widget->MessageTextEdit->appendPlainText("A new empty handle group is created.");
  }

  void delete_handle_group()
  {
    if(!is_there_any_handle_group()) { 
      ui_widget->MessageTextEdit->appendPlainText("There is no handle group to be deleted!");
      return; 
    } // no handle group

    // first erase handles one by one
    Deform_mesh::Handle_iterator hb, he;
    boost::tie(hb, he) = deform_mesh.handles(active_group);
    while(hb != he)
    {
      Deform_mesh::Handle_iterator tmp = hb; ++tmp;      
      erase_handle(*hb);
      hb = tmp;
    }
    // now delete group representative    
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


    Deform_mesh::Handle_group hgb, hge;
    if( !is_there_any_handle_group(hgb, hge) )
    { return; } // no handle group 
    active_group = hgb;    
  }

  void prev_handle_group()
  {
    Deform_mesh::Handle_group hgb, hge;
    if( !is_there_any_handle_group(hgb, hge) ) {
      ui_widget->MessageTextEdit->appendPlainText("There is no handle group to iterate on!");
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
      ui_widget->MessageTextEdit->appendPlainText("There is no handle group to iterate on!");
      return; 
    }
    // shift
    if(--hge == active_group) { active_group = hgb; }
    else                      {++active_group; }    
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

  void process_selection(vertex_descriptor v, int k_ring, bool is_roi, bool is_insert)
  {
    std::cout << "Process k-ring: " << k_ring << " roi: " << is_roi << " insert: " << is_insert << std::endl;

    std::map<vertex_descriptor, int> neighs = extract_k_ring(*polyhedron(), v, k_ring);
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

  void refresh_manipulated_frame_center(Deform_mesh::Handle_group hg)
  {
    qglviewer::Vec center = calculate_center(hg);
    Handle_group_data& hd = get_data(hg);

    hd.initial_center = center;

    qglviewer::ManipulatedFrame* hg_frame = hd.frame;
    hg_frame->blockSignals(true); // do not let it emit modified, which will cause a deformation
                                  // but we are just adjusting the center so it does not require a deformation
    hg_frame->setPosition(center);
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
    }
    if(counter == 0) { return qglviewer::Vec(0,0,0); } 
    return qglviewer::Vec(center_acc.x() / counter, center_acc.y() / counter, center_acc.z() / counter);
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
}; // end class Scene_edit_polyhedron_item_2

#endif // SCENE_EDIT_POLYHEDRON_ITEM_2_H
