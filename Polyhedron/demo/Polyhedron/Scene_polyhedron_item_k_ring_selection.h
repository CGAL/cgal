#ifndef SCENE_POLYHEDRON_ITEM_K_RING_SELECTION_H
#define SCENE_POLYHEDRON_ITEM_K_RING_SELECTION_H
#include "Scene_polyhedron_item_k_ring_selection_config.h"
#include "Polyhedron_type.h"
#include "Scene_polyhedron_item.h"

#include <QGLViewer/qglviewer.h>
#include <QKeyEvent>
#include <QMouseEvent>
#include <QMainWindow>
#include <QObject>

#include <map>
#include <queue>
#include "One_ring_iterators.h"

class SCENE_POLYHEDRON_ITEM_K_RING_SELECTION_EXPORT Scene_polyhedron_item_k_ring_selection 
  : public QObject
{
  Q_OBJECT
public:
  struct Active_handle { enum Type{ VERTEX = 0, FACET = 1, EDGE = 2 }; };

  // Hold mouse keyboard state together
  struct Mouse_keyboard_state
  {
    Mouse_keyboard_state() : shift_pressing(false), left_button_pressing(false) { }
    bool shift_pressing, left_button_pressing;
  };

  Mouse_keyboard_state  state;

  Active_handle::Type    active_handle_type;
  int                    k_ring;
  Scene_polyhedron_item* poly_item;

  Scene_polyhedron_item_k_ring_selection() {}

  Scene_polyhedron_item_k_ring_selection
    (Scene_polyhedron_item* poly_item, QMainWindow* mw, Active_handle::Type aht, int k_ring)
  {
    init(poly_item, mw, aht, k_ring);
  }

  void init(Scene_polyhedron_item* poly_item, QMainWindow* mw, Active_handle::Type aht, int k_ring) {
    this->poly_item = poly_item;
    this->active_handle_type = aht;
    this->k_ring = k_ring;

    poly_item->enable_facets_picking(true);
    poly_item->set_color_vector_read_only(true);

    QGLViewer* viewer = *QGLViewer::QGLViewerPool().begin();
    viewer->installEventFilter(this);
    mw->installEventFilter(this);

    connect(poly_item, SIGNAL(selected_vertex(void*)), this, SLOT(vertex_has_been_selected(void*)));
    connect(poly_item, SIGNAL(selected_facet(void*)), this, SLOT(facet_has_been_selected(void*)));
    connect(poly_item, SIGNAL(selected_edge(void*)), this, SLOT(edge_has_been_selected(void*)));
  }

public slots:
  // slots are called by signals of polyhedron_item
  void vertex_has_been_selected(void* void_ptr) 
  {
    if(active_handle_type != Active_handle::VERTEX) { return; }
    process_selection( static_cast<Polyhedron::Vertex*>(void_ptr)->halfedge()->vertex() );
  }
  void facet_has_been_selected(void* void_ptr)
  {
    if(active_handle_type != Active_handle::FACET) { return; }
    process_selection( static_cast<Polyhedron::Facet*>(void_ptr)->halfedge()->facet() );
  }
  void edge_has_been_selected(void* void_ptr) 
  {
    if(active_handle_type != Active_handle::EDGE) { return; }
    process_selection( static_cast<Polyhedron::Halfedge*>(void_ptr)->opposite()->opposite() );
  }

signals:
  void selected(const std::map<Polyhedron::Vertex_handle, int>&);
  void selected(const std::map<Polyhedron::Facet_handle, int>&);
  void selected(const std::map<Polyhedron::Halfedge_handle, int>&);

protected:
  template<class HandleType>
  void process_selection(HandleType clicked) {
    const std::map<HandleType, int>& selection = extract_k_ring(clicked, k_ring);
    emit selected(selection);
  }

  template<class HandleType>
  std::map<HandleType, int> extract_k_ring(HandleType v, int k)
  {
    std::map<HandleType, int>  D;
    std::queue<HandleType>     Q;
    Q.push(v); D[v] = 0;

    int dist_v;
    while( !Q.empty() && (dist_v = D[Q.front()]) < k ) {
      v = Q.front();
      Q.pop();

      for(One_ring_iterator<HandleType> circ(v); circ; ++circ)
      {
        HandleType new_v = circ;
        if(D.insert(std::make_pair(new_v, dist_v + 1)).second) {
          Q.push(new_v);
        }
      }
    }
    return D;
  }

  bool eventFilter(QObject* /*target*/, QEvent *event)
  {
    // This filter is both filtering events from 'viewer' and 'main window'
    // key events
    if(event->type() == QEvent::KeyPress || event->type() == QEvent::KeyRelease)  {
      QKeyEvent *keyEvent = static_cast<QKeyEvent*>(event);
      Qt::KeyboardModifiers modifiers = keyEvent->modifiers();

      state.shift_pressing = modifiers.testFlag(Qt::ShiftModifier);
    }
    // mouse events
    if(event->type() == QEvent::MouseButtonPress || event->type() == QEvent::MouseButtonRelease) {
      QMouseEvent* mouse_event = static_cast<QMouseEvent*>(event);
      if(mouse_event->button() == Qt::LeftButton) {
        state.left_button_pressing = event->type() == QEvent::MouseButtonPress;
      }   
    }

    // use mouse move event for paint-like selection
    if(event->type() == QEvent::MouseMove &&
      (state.shift_pressing && state.left_button_pressing) )    
    { // paint with mouse move event 
      QMouseEvent* mouse_event = static_cast<QMouseEvent*>(event);
      QGLViewer* viewer = *QGLViewer::QGLViewerPool().begin();
      qglviewer::Camera* camera = viewer->camera();

      bool found = false;
      const qglviewer::Vec& point = camera->pointUnderPixel(mouse_event->pos(), found);
      if(found)
      {
        const qglviewer::Vec& orig = camera->position();
        const qglviewer::Vec& dir = point - orig;
        poly_item->select(orig.x, orig.y, orig.z, dir.x, dir.y, dir.z);
      }
    }//end MouseMove
    return false;
  }
};

#endif
