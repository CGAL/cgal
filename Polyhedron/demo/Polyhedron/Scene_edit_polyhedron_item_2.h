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

public slots:
  void changed();
  void select(double orig_x,
              double orig_y,
              double orig_z,
              double dir_x,
              double dir_y,
              double dir_z);
  
  void vertex_has_been_selected(void* void_ptr);

public:
  Ui::DeformMesh_2* ui_widget;

  Scene_polyhedron_item* poly_item;
  qglviewer::ManipulatedFrame* frame;

  std::set<vertex_descriptor> roi;
  std::set<vertex_descriptor> handles;

  Deform_mesh deform_mesh;
  Deform_mesh::Handle_group active_group;

public:
 void insert_handle(vertex_descriptor v)
 {
    if(handles.insert(v).second)
    {
      deform_mesh.insert_handle(active_group, v);
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

  void create_handle_group()
  {
    active_group = deform_mesh.create_handle_group();
  }
  
 void erase_handle(vertex_descriptor v)
 {
    if(handles.erase(v) == 1)
    {
      deform_mesh.erase_handle(active_group, v);
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

  void prev_handle_group()
  {
    Deform_mesh::Handle_group hgb, hge;
    boost::tie(hgb, hge) = deform_mesh.handle_groups();
    if(hgb == active_group) { active_group = hge; }
    else                    {--active_group; }    
  }

  void next_handle_group()
  {
    Deform_mesh::Handle_group hgb, hge;
    boost::tie(hgb, hge) = deform_mesh.handle_groups();
    if(hge == active_group) { active_group = hgb; }
    else                    {++active_group; }    
  }

  void process_selection(vertex_descriptor v, int k_ring, bool is_roi, bool is_insert)
  {
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
