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

#include <CGAL/boost/graph/selection.h>

class SCENE_POLYHEDRON_ITEM_K_RING_SELECTION_EXPORT Scene_polyhedron_item_k_ring_selection 
  : public QObject
{
  Q_OBJECT
public:
  struct Active_handle { enum Type{ VERTEX = 0, FACET = 1, EDGE = 2 }; };

  typedef boost::graph_traits<Polyhedron>::edge_descriptor edge_descriptor;

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
  bool is_active;

  Scene_polyhedron_item_k_ring_selection() {}

  Scene_polyhedron_item_k_ring_selection
    (Scene_polyhedron_item* poly_item, QMainWindow* mw, Active_handle::Type aht, int k_ring)
      :is_active(false)
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

public Q_SLOTS:
  // slots are called by signals of polyhedron_item
  void vertex_has_been_selected(void* void_ptr) 
  {
    is_active=true;
    if(active_handle_type != Active_handle::VERTEX) { return; }
    process_selection( static_cast<Polyhedron::Vertex*>(void_ptr)->halfedge()->vertex() );
  }
  void facet_has_been_selected(void* void_ptr)
  {
    is_active=true;
    if(active_handle_type != Active_handle::FACET) { return; }
    process_selection( static_cast<Polyhedron::Facet*>(void_ptr)->halfedge()->facet() );
  }
  void edge_has_been_selected(void* void_ptr) 
  {
    is_active=true;
    if(active_handle_type != Active_handle::EDGE) { return; }
    process_selection( edge(static_cast<Polyhedron::Halfedge*>(void_ptr)->opposite()->opposite(), *poly_item->polyhedron()) );
  }

Q_SIGNALS:
  void selected(const std::set<Polyhedron::Vertex_handle>&);
  void selected(const std::set<Polyhedron::Facet_handle>&);
  void selected(const std::set<edge_descriptor>&);
  void endSelection();

protected:

  template<class HandleType>
  void process_selection(HandleType clicked) {
    const std::set<HandleType>& selection = extract_k_ring(clicked, k_ring);
    Q_EMIT selected(selection);
  }

  template <class Handle>
  struct Is_selected_from_set{
    std::set<Handle>& selection;
    Is_selected_from_set(std::set<Handle>& selection)
      :selection(selection) {}
    friend bool get(Is_selected_from_set<Handle> map, Handle k)
    {
      return map.selection.count(k);
    }
    friend void put(Is_selected_from_set<Handle> map, Handle k, bool b)
    {
      if (b)
        map.selection.insert(k);
      else
        map.selection.erase(k);
    }
  };

  std::set<Polyhedron::Vertex_handle>
  extract_k_ring(Polyhedron::Vertex_handle clicked, unsigned int k)
  {
    std::set<Polyhedron::Vertex_handle> selection;
    selection.insert(clicked);
    if (k>0)
      CGAL::dilate_vertex_selection(CGAL::make_array(clicked),
                                    *poly_item->polyhedron(),
                                    k,
                                    Is_selected_from_set<Polyhedron::Vertex_handle>(selection),
                                    CGAL::Emptyset_iterator());

    return selection;
  }

  std::set<Polyhedron::Facet_handle>
  extract_k_ring(Polyhedron::Facet_handle clicked, unsigned int k)
  {
    std::set<Polyhedron::Facet_handle> selection;
    selection.insert(clicked);
    if (k>0)
      CGAL::dilate_face_selection(CGAL::make_array(clicked),
                                  *poly_item->polyhedron(),
                                  k,
                                  Is_selected_from_set<Polyhedron::Facet_handle>(selection),
                                  CGAL::Emptyset_iterator());

    return selection;
  }

  std::set<edge_descriptor>
  extract_k_ring(edge_descriptor clicked, unsigned int k)
  {
    std::set<edge_descriptor> selection;
    selection.insert(clicked);
    if (k>0)
      CGAL::dilate_edge_selection(CGAL::make_array(clicked),
                                  *poly_item->polyhedron(),
                                  k,
                                  Is_selected_from_set<edge_descriptor>(selection),
                                  CGAL::Emptyset_iterator());

    return selection;
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
        if (!state.left_button_pressing)
          if (is_active)
          {
            Q_EMIT endSelection();
            is_active=false;
          }
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
