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
  struct Active_handle {
    enum Type{ VERTEX = 0, FACET = 1, EDGE = 2 , CONNECTED_COMPONENT = 3, PATH = 4};
  };

  typedef boost::graph_traits<Polyhedron>::edge_descriptor edge_descriptor;

  // Hold mouse keyboard state together
  struct Mouse_keyboard_state
  {
    Mouse_keyboard_state() : shift_pressing(false), left_button_pressing(false) { }
    bool shift_pressing, left_button_pressing;
  };

  Mouse_keyboard_state  state;
  QMainWindow* mainwindow;
  Active_handle::Type    active_handle_type;
  int                    k_ring;
  Scene_polyhedron_item* poly_item;
  bool is_active;
  bool is_current_selection;
  bool is_highlighting;

  Scene_polyhedron_item_k_ring_selection() {}

  Scene_polyhedron_item_k_ring_selection
    (Scene_polyhedron_item* poly_item, QMainWindow* mw, Active_handle::Type aht, int k_ring)
      :is_active(false),is_current_selection(true), is_edit_mode(false)
  {
    init(poly_item, mw, aht, k_ring);
  }

  void setEditMode(bool b)
  {
    is_edit_mode = b;
    QGLViewer* viewer = *QGLViewer::QGLViewerPool().begin();
    //for highlighting
    viewer->setMouseTracking(b);
  }

  void init(Scene_polyhedron_item* poly_item, QMainWindow* mw, Active_handle::Type aht, int k_ring) {
    this->poly_item = poly_item;
    this->active_handle_type = aht;
    this->k_ring = k_ring;
    mainwindow = mw;
    is_highlighting = false;
    is_ready_to_highlight = true;
    is_ready_to_paint_select = true;
    poly_item->enable_facets_picking(true);
    poly_item->set_color_vector_read_only(true);

    QGLViewer* viewer = *QGLViewer::QGLViewerPool().begin();
    viewer->installEventFilter(this);
    mw->installEventFilter(this);
#if QGLVIEWER_VERSION >= 0x020501
    viewer->setMouseBindingDescription(Qt::Key_D, Qt::ShiftModifier, Qt::LeftButton, "(When in selection plugin) Removes the clicked primitive from the selection. ");
#else
    viewer->setMouseBindingDescription(Qt::SHIFT + Qt::LeftButton,  "(When in selection plugin) When D is pressed too, removes the clicked primitive from the selection. ");
#endif
    connect(poly_item, SIGNAL(selected_vertex(void*)), this, SLOT(vertex_has_been_selected(void*)));
    connect(poly_item, SIGNAL(selected_facet(void*)), this, SLOT(facet_has_been_selected(void*)));
    connect(poly_item, SIGNAL(selected_edge(void*)), this, SLOT(edge_has_been_selected(void*)));
  }
  void setCurrentlySelected(bool b)
  {
    is_current_selection = b;
  }
public Q_SLOTS:
  // slots are called by signals of polyhedron_item
  void vertex_has_been_selected(void* void_ptr) 
  {
    is_active=true;
    if(active_handle_type == Active_handle::VERTEX || active_handle_type == Active_handle::PATH)
      process_selection( static_cast<Polyhedron::Vertex*>(void_ptr)->halfedge()->vertex() );
    updateIsTreated();
  }
  void facet_has_been_selected(void* void_ptr)
  {
    is_active=true;
    if (active_handle_type == Active_handle::FACET
      || active_handle_type == Active_handle::CONNECTED_COMPONENT)
      process_selection(static_cast<Polyhedron::Facet*>(void_ptr)->halfedge()->facet());
    updateIsTreated();
  }
  void edge_has_been_selected(void* void_ptr) 
  {
    is_active=true;
    if(active_handle_type == Active_handle::EDGE)
      process_selection( edge(static_cast<Polyhedron::Halfedge*>(void_ptr)->opposite()->opposite(), *poly_item->polyhedron()) );
    updateIsTreated();
  }

  void paint_selection()
  {
    if(is_ready_to_paint_select)
    {
      // paint with mouse move event
      QGLViewer* viewer = *QGLViewer::QGLViewerPool().begin();
      qglviewer::Camera* camera = viewer->camera();

      bool found = false;
      const qglviewer::Vec& point = camera->pointUnderPixel(paint_pos, found);
      if(found)
      {
        const qglviewer::Vec& orig = camera->position();
        const qglviewer::Vec& dir = point - orig;
        poly_item->select(orig.x, orig.y, orig.z, dir.x, dir.y, dir.z);
      }
      is_ready_to_paint_select = false;
    }
  }

  void highlight()
  {
    if(is_ready_to_highlight)
    {
      // highlight with mouse move event
      QGLViewer* viewer = *QGLViewer::QGLViewerPool().begin();
      qglviewer::Camera* camera = viewer->camera();

      bool found = false;
      const qglviewer::Vec& point = camera->pointUnderPixel(hl_pos, found);
      if(found)
      {
        const qglviewer::Vec& orig = camera->position();
        const qglviewer::Vec& dir = point - orig;
        is_highlighting = true;
        poly_item->select(orig.x, orig.y, orig.z, dir.x, dir.y, dir.z);
        is_highlighting = false;
      }
      else
      {
        Q_EMIT clearHL();
      }
      is_ready_to_highlight = false;
    }
  }

Q_SIGNALS:
  void selected(const std::set<Polyhedron::Vertex_handle>&);
  void selected(const std::set<Polyhedron::Facet_handle>&);
  void selected(const std::set<edge_descriptor>&);
  void selected_HL(const std::set<Polyhedron::Vertex_handle>&);
  void selected_HL(const std::set<Polyhedron::Facet_handle>&);
  void selected_HL(const std::set<edge_descriptor>&);
  void toogle_insert(const bool);
  void endSelection();
  void resetIsTreated(); 
  void isCurrentlySelected(Scene_polyhedron_item_k_ring_selection*);
  void clearHL();

protected:

  void updateIsTreated()
  {
    static ushort i = 0;
    i++;
    if(i==3)
    {
      i = 0;
      Q_EMIT resetIsTreated();
    }
  }
  template<class HandleType>
  void process_selection(HandleType clicked) {
    const std::set<HandleType>& selection = extract_k_ring(clicked, k_ring);
    if(is_highlighting)
    {
      Q_EMIT selected_HL(selection);
    }
    else
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
      CGAL::expand_vertex_selection(CGAL::make_array(clicked),
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
      CGAL::expand_face_selection(CGAL::make_array(clicked),
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
      CGAL::expand_edge_selection(CGAL::make_array(clicked),
                                  *poly_item->polyhedron(),
                                  k,
                                  Is_selected_from_set<edge_descriptor>(selection),
                                  CGAL::Emptyset_iterator());
    return selection;
  }


  bool eventFilter(QObject* target, QEvent *event)
  {
    // This filter is both filtering events from 'viewer' and 'main window'

    // key events
      if(event->type() == QEvent::KeyPress || event->type() == QEvent::KeyRelease)  {
      QKeyEvent *keyEvent = static_cast<QKeyEvent*>(event);
      Qt::KeyboardModifiers modifiers = keyEvent->modifiers();

      state.shift_pressing = modifiers.testFlag(Qt::ShiftModifier);
    }

    if(event->type() == QEvent::KeyPress
            && state.shift_pressing
            && static_cast<QKeyEvent*>(event)->key()==Qt::Key_D)
    {
     Q_EMIT toogle_insert(false);
    }
    else if(event->type() == QEvent::KeyRelease
            && static_cast<QKeyEvent*>(event)->key()==Qt::Key_D)
    {
     Q_EMIT toogle_insert(true);
    }
    // mouse events
    if(event->type() == QEvent::MouseButtonPress || event->type() == QEvent::MouseButtonRelease) {
      QMouseEvent* mouse_event = static_cast<QMouseEvent*>(event);
      if(mouse_event->button() == Qt::LeftButton) {
        state.left_button_pressing = event->type() == QEvent::MouseButtonPress;
        if(is_edit_mode)
          Q_EMIT clearHL();
        if (!state.left_button_pressing)
          if (is_active)
          {
            Q_EMIT endSelection();
            is_active=false;
          }
      }
      //to avoid the contextual menu to mess up the states.
      else if(mouse_event->button() == Qt::RightButton) {
        state.left_button_pressing = false;
        state.shift_pressing = false;
      }
    }
    // use mouse move event for paint-like selection
    if( (event->type() == QEvent::MouseMove
         || (event->type() == QEvent::MouseButtonPress
             && static_cast<QMouseEvent*>(event)->button() == Qt::LeftButton))
        && (state.shift_pressing && state.left_button_pressing) )
    {
      Q_EMIT isCurrentlySelected(this);
      if(!is_current_selection)
        return false;
      if(target == mainwindow)
      {
        QGLViewer* viewer = *QGLViewer::QGLViewerPool().begin();
        viewer->setFocus();
        return false;
      }
      is_ready_to_paint_select = true;
      QMouseEvent* mouse_event = static_cast<QMouseEvent*>(event);
      paint_pos = mouse_event->pos();
      if(!is_edit_mode || event->type() == QEvent::MouseButtonPress)
        QTimer::singleShot(0,this,SLOT(paint_selection()));
    }
    //if in edit_mode and the mouse is moving without left button pressed :
    // highlight the primitive under cursor
    else if(is_edit_mode && event->type() == QEvent::MouseMove && !state.left_button_pressing)
    {
      if(target == mainwindow)
      {
        QGLViewer* viewer = *QGLViewer::QGLViewerPool().begin();
        viewer->setFocus();
        return false;
      }

      is_ready_to_highlight = true;
      QMouseEvent* mouse_event = static_cast<QMouseEvent*>(event);
      hl_pos = mouse_event->pos();
      QTimer::singleShot(0, this, SLOT(highlight()));
    }//end MouseMove
    return false;
  }
  bool is_edit_mode;
  bool is_ready_to_highlight;
  bool is_ready_to_paint_select;
  QPoint hl_pos;
  QPoint paint_pos;

};

#endif
