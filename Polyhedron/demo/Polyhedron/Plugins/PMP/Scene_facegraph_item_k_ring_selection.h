#ifndef SCENE_FACEGRAPH_ITEM_K_RING_SELECTION_H
#define SCENE_FACEGRAPH_ITEM_K_RING_SELECTION_H
#include "Scene_facegraph_item_k_ring_selection_config.h"
#ifdef USE_SURFACE_MESH
#include "Scene_surface_mesh_item.h"
#include <CGAL/iterator.h>
#else
#include "Scene_polyhedron_item.h"
#include "Polyhedron_type.h"
#endif
#include <set>
#include <QGLViewer/qglviewer.h>
#include <QKeyEvent>
#include <QMouseEvent>
#include <QMainWindow>
#include <QObject>

#include <map>
#include <queue>

#include <CGAL/boost/graph/selection.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/border.h>

#include <CGAL/Polygon_2.h>

#ifdef USE_SURFACE_MESH
typedef Scene_surface_mesh_item Scene_facegraph_item;
typedef Scene_surface_mesh_item::Kernel FG_Traits;
#else
typedef Scene_polyhedron_item Scene_facegraph_item;
typedef Kernel FG_Traits;
#endif

typedef Scene_facegraph_item::FaceGraph FaceGraph;
typedef boost::graph_traits<FaceGraph>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<FaceGraph>::edge_descriptor edge_descriptor;
typedef boost::graph_traits<FaceGraph>::face_descriptor face_descriptor;
typedef boost::graph_traits<FaceGraph>::halfedge_descriptor halfedge_descriptor;

struct FG_is_selected_edge_property_map{
  typedef boost::property_map<FaceGraph,boost::edge_index_t>::type EImap;

  std::vector<bool>* is_selected_ptr;
  EImap* edge_index_map;
  FG_is_selected_edge_property_map()
    : is_selected_ptr(NULL), edge_index_map(NULL) {}
  FG_is_selected_edge_property_map(std::vector<bool>& is_selected, EImap* map)
    : is_selected_ptr( &is_selected), edge_index_map(map)
  {}

  std::size_t id(edge_descriptor ed) {
    return get(*edge_index_map, ed);
  }

  friend bool get(FG_is_selected_edge_property_map map, edge_descriptor ed)
  {
    CGAL_assertion(map.is_selected_ptr!=NULL);
    return (*map.is_selected_ptr)[map.id(ed)];
  }

  friend void put(FG_is_selected_edge_property_map map, edge_descriptor ed, bool b)
  {
    CGAL_assertion(map.is_selected_ptr!=NULL);
    (*map.is_selected_ptr)[map.id(ed)]=b;
  }
};

class SCENE_FACEGRAPH_ITEM_K_RING_SELECTION_EXPORT Scene_facegraph_item_k_ring_selection
  : public QObject
{
  Q_OBJECT
public:
  struct Active_handle {
    enum Type{ VERTEX = 0, FACET = 1, EDGE = 2 , CONNECTED_COMPONENT = 3, PATH = 4};
  };

  typedef CGAL::Polygon_2<FG_Traits> Polygon_2;
  typedef std::vector<FG_Traits::Point_2> Polyline_2;
  typedef std::vector<Polyline_2> Polylines;

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
  Scene_facegraph_item* poly_item;
  bool is_active;
  bool is_current_selection;
  bool is_highlighting;

  Scene_facegraph_item_k_ring_selection() {}

  Scene_facegraph_item_k_ring_selection
    (Scene_facegraph_item* poly_item, QMainWindow* mw, Active_handle::Type aht, int k_ring)
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

  void init(Scene_facegraph_item* poly_item, QMainWindow* mw, Active_handle::Type aht, int k_ring) {
    this->poly_item = poly_item;
    this->active_handle_type = aht;
    this->k_ring = k_ring;
    polyline = new Polylines(0);
    polyline->push_back(Polyline_2());
    mainwindow = mw;
    is_highlighting = false;
    is_ready_to_highlight = true;
    is_ready_to_paint_select = true;
    is_lasso_active = false;

#ifndef USE_SURFACE_MESH
    poly_item->enable_facets_picking(true);
    poly_item->set_color_vector_read_only(true);
#endif
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
  void set_lasso_mode(bool b) { is_lasso_active = b; }

public Q_SLOTS:
  // slots are called by signals of polyhedron_item
  void vertex_has_been_selected(void* void_ptr) 
  {
    is_active=true;
    if(active_handle_type == Active_handle::VERTEX || active_handle_type == Active_handle::PATH)
    {
#ifdef USE_SURFACE_MESH
      std::size_t h = reinterpret_cast<std::size_t>(void_ptr);
      process_selection( static_cast<vertex_descriptor>(h) );
#else
      process_selection( static_cast<Polyhedron::Vertex*>(void_ptr)->halfedge()->vertex() );
#endif
    }
    updateIsTreated();
  }
  void facet_has_been_selected(void* void_ptr)
  {
    is_active=true;
    if (active_handle_type == Active_handle::FACET
      || active_handle_type == Active_handle::CONNECTED_COMPONENT)
    {
#ifdef USE_SURFACE_MESH
      std::size_t h = reinterpret_cast<std::size_t>(void_ptr);
      process_selection( static_cast<face_descriptor>(h) );
#else
      process_selection( static_cast<Polyhedron::Facet*>(void_ptr)->halfedge()->facet() );
#endif
    }
    updateIsTreated();
  }
  void edge_has_been_selected(void* void_ptr) 
  {
    is_active=true;
    if(active_handle_type == Active_handle::EDGE)
    {
#ifdef USE_SURFACE_MESH
      std::size_t h = reinterpret_cast<std::size_t>(void_ptr);
      process_selection( static_cast<edge_descriptor>(h) );
#else
      process_selection( edge(static_cast<Polyhedron::Halfedge*>(void_ptr)->opposite()->opposite(), *poly_item->polyhedron()) );
#endif
    }
    updateIsTreated();
  }

  void paint_selection()
  {
    if(is_ready_to_paint_select)
    {
      const qglviewer::Vec offset = static_cast<CGAL::Three::Viewer_interface*>(QGLViewer::QGLViewerPool().first())->offset();
      // paint with mouse move event
      QGLViewer* viewer = *QGLViewer::QGLViewerPool().begin();
      qglviewer::Camera* camera = viewer->camera();

      bool found = false;
      const qglviewer::Vec& point = camera->pointUnderPixel(paint_pos, found) - offset;
      if(found)
      {
        const qglviewer::Vec& orig = camera->position() - offset;
        const qglviewer::Vec& dir = point - orig;
        poly_item->select(orig.x, orig.y, orig.z, dir.x, dir.y, dir.z);
      }
      is_ready_to_paint_select = false;
    }
  }

  void lasso_selection()
  {
    QGLViewer* viewer = *QGLViewer::QGLViewerPool().begin();
    const qglviewer::Vec offset = static_cast<CGAL::Three::Viewer_interface*>(viewer)->offset();

    qglviewer::Camera* camera = viewer->camera();
    const FaceGraph& poly = *poly_item->polyhedron();

    std::set<face_descriptor> face_sel;
    boost::property_map<FaceGraph,CGAL::vertex_point_t>::const_type vpmap = get(boost::vertex_point, poly);
    //select all faces if their screen projection is inside the lasso
    BOOST_FOREACH(face_descriptor f, faces(poly))
    {
      BOOST_FOREACH(vertex_descriptor v, CGAL::vertices_around_face(halfedge(f, poly), poly))
      {
        FG_Traits::Point_3 p = get(vpmap, v);
        qglviewer::Vec vp(p.x(), p.y(), p.z());
        qglviewer::Vec vsp = camera->projectedCoordinatesOf(vp+offset);
        if(is_vertex_selected(vsp))
        {
          face_sel.insert(f);
          break;
        }
      }
    }
    if(face_sel.empty())
    {
      contour_2d.clear();
      return;
    }
    //get border edges of the selected patches
    std::vector<halfedge_descriptor> boundary_edges;
    CGAL::Polygon_mesh_processing::border_halfedges(face_sel, poly, std::back_inserter(boundary_edges));
    std::vector<bool> mark(edges(poly).size(), false);
    boost::property_map<FaceGraph, boost::edge_index_t>::type edge_index
      = get(boost::edge_index, poly);
    FG_is_selected_edge_property_map spmap(mark, &edge_index);
    BOOST_FOREACH(halfedge_descriptor h, boundary_edges)
      put(spmap, edge(h, poly), true);

    boost::vector_property_map<int,
      boost::property_map<FaceGraph, boost::face_index_t>::type>
      fccmap;

    //get connected componant from the picked face
    std::set<face_descriptor> final_sel;
    //std::vector<Polyhedron::Face_handle> cc;
    std::size_t nb_cc = CGAL::Polygon_mesh_processing::connected_components(poly
          , fccmap
          , CGAL::Polygon_mesh_processing::parameters::edge_is_constrained_map(spmap));
    std::vector<bool> is_cc_done(nb_cc, false);

    BOOST_FOREACH(face_descriptor f, face_sel)
    {

      int cc_id = get(fccmap, f);
      if(is_cc_done[cc_id])
      {
        continue;
      }
      CGAL::Halfedge_around_face_circulator<FaceGraph> hafc(halfedge(f, poly), poly);
      CGAL::Halfedge_around_face_circulator<FaceGraph> end = hafc;
      double x(0), y(0), z(0);
      int total(0);
      CGAL_For_all(hafc, end)
      {
        FG_Traits::Point_3 p = get(vpmap, target(*hafc, poly));
        x+=p.x(); y+=p.y(); z+=p.z();
        total++;
      }
      if(total == 0)
        continue;
      qglviewer::Vec center(x/(double)total, y/(double)total, z/(double)total);
      const qglviewer::Vec& orig = camera->position() - offset;
      qglviewer::Vec direction = center - orig;
      if(poly_item->intersect_face(orig.x,
                                   orig.y,
                                   orig.z,
                                   direction.x,
                                   direction.y,
                                   direction.z,
                                   f))
      {
        is_cc_done[cc_id] = true;
      }
    }
    BOOST_FOREACH(face_descriptor f, faces(poly))
    {
      if(is_cc_done[get(fccmap, f)])
        final_sel.insert(f);
    }
    switch(active_handle_type)
    {
    case Active_handle::FACET:
      selected(final_sel);
      break;
    case Active_handle::EDGE:
    {
      std::set<edge_descriptor> e_sel;
      BOOST_FOREACH(face_descriptor f, final_sel)
      {
        BOOST_FOREACH(halfedge_descriptor h, CGAL::halfedges_around_face(halfedge(f, poly), poly))
        {
          FG_Traits::Point_3 p = get(vpmap, target(h, poly));
          qglviewer::Vec vp1(p.x(), p.y(), p.z());
          qglviewer::Vec vsp1 = camera->projectedCoordinatesOf(vp1+offset);
          p = get(vpmap, target(opposite(h, poly), poly));
          qglviewer::Vec vp2(p.x(), p.y(), p.z());
          qglviewer::Vec vsp2 = camera->projectedCoordinatesOf(vp2+offset);
          if(is_vertex_selected(vsp1) || is_vertex_selected(vsp2))
            e_sel.insert(edge(h, poly));
        }
      }
      selected(e_sel);
      break;
    }
    case Active_handle::VERTEX:
    {
      std::set<vertex_descriptor> v_sel;
      BOOST_FOREACH(face_descriptor f, final_sel)
      {
        BOOST_FOREACH(vertex_descriptor v, CGAL::vertices_around_face(halfedge(f, poly), poly))
        {
          FG_Traits::Point_3 p = get(vpmap, v);
          qglviewer::Vec vp(p.x(), p.y(), p.z());
          qglviewer::Vec vsp = camera->projectedCoordinatesOf(vp+offset);
          if(is_vertex_selected(vsp))
            v_sel.insert(v);
        }
      }
      selected(v_sel);
      break;
    }
    default:
      break;
    }
    contour_2d.clear();
  }

  void highlight()
  {
    const qglviewer::Vec offset = static_cast<CGAL::Three::Viewer_interface*>(QGLViewer::QGLViewerPool().first())->offset();
    if(is_ready_to_highlight)
    {
      // highlight with mouse move event
      QGLViewer* viewer = *QGLViewer::QGLViewerPool().begin();
      qglviewer::Camera* camera = viewer->camera();
      bool found = false;
      const qglviewer::Vec& point = camera->pointUnderPixel(hl_pos, found) - offset;
      if(found)
      {
        const qglviewer::Vec& orig = camera->position() - offset;
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
  void selected(const std::set<vertex_descriptor>&);
  void selected(const std::set<face_descriptor>&);
  void selected(const std::set<edge_descriptor>&);
  void selected_HL(const std::set<vertex_descriptor>&);
  void selected_HL(const std::set<face_descriptor>&);
  void selected_HL(const std::set<edge_descriptor>&);
  void toogle_insert(const bool);
  void endSelection();
  void resetIsTreated(); 
  void isCurrentlySelected(Scene_facegraph_item_k_ring_selection*);
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
    //keeps the highlighting on track if the brush_size is not 0
    int current_ring = 0;
    if(active_handle_type != Active_handle::PATH && !is_edit_mode)
      current_ring = k_ring;
    const std::set<HandleType>& selection = extract_k_ring(clicked, current_ring);
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

  std::set<vertex_descriptor>
  extract_k_ring(vertex_descriptor clicked, unsigned int k)
  {
    std::set<vertex_descriptor> selection;
    selection.insert(clicked);
    if (k>0)
      CGAL::expand_vertex_selection(CGAL::make_array(clicked),
                                    *poly_item->polyhedron(),
                                    k,
                                    Is_selected_from_set<vertex_descriptor>(selection),
                                    CGAL::Emptyset_iterator());

    return selection;
  }

  std::set<face_descriptor>
  extract_k_ring(face_descriptor clicked, unsigned int k)
  {
    std::set<face_descriptor> selection;
    selection.insert(clicked);
    if (k>0)
      CGAL::expand_face_selection(CGAL::make_array(clicked),
                                  *poly_item->polyhedron(),
                                  k,
                                  Is_selected_from_set<face_descriptor>(selection),
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
        {
          if (is_active)
          {
            Q_EMIT endSelection();
            is_active=false;
          }
          apply_path();
          lasso_selection();
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
      if(!is_lasso_active)
      {
        is_ready_to_paint_select = true;
        QMouseEvent* mouse_event = static_cast<QMouseEvent*>(event);
        paint_pos = mouse_event->pos();
        if(!is_edit_mode || event->type() == QEvent::MouseButtonPress)
          QTimer::singleShot(0,this,SLOT(paint_selection()));
      }
      else
      {
        sample_mouse_path();
      }
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
  bool is_lasso_active;
  QPoint hl_pos;
  QPoint paint_pos;
  Polyline_2 contour_2d;
  Polylines* polyline;
  Polyline_2& poly() const  { return polyline->front(); }
  Polygon_2 lasso;
  CGAL::Bbox_2 domain_rectangle;
  bool update_polyline () const
  {
    if (contour_2d.size() < 2 ||
        (!(poly().empty()) && contour_2d.back () == poly().back()))
      return false;


    if (!(poly().empty()) && contour_2d.back () == poly().back())
      return false;

    poly().clear();

    for (unsigned int i = 0; i < contour_2d.size (); ++ i)
      poly().push_back (contour_2d[i]);

    return true;
  }

  void sample_mouse_path()
  {
    CGAL::Three::Viewer_interface* viewer = static_cast<CGAL::Three::Viewer_interface*>(*QGLViewer::QGLViewerPool().begin());
    const QPoint& p = viewer->mapFromGlobal(QCursor::pos());
    contour_2d.push_back (FG_Traits::Point_2 (p.x(), p.y()));

    if (update_polyline ())
    {
      //update draw
      QPainter *painter = viewer->getPainter();
      QPen pen;
      pen.setColor(QColor(Qt::green));
      pen.setWidth(3);
      //Create a QImage of the screen and paint the lasso on top of it
      QImage image = viewer->grabFrameBuffer();
      painter->begin(viewer);
      painter->drawImage(QPoint(0,0), image);
      painter->setPen(pen);
      for(std::size_t i=0; i<polyline->size(); ++i)
      {
        Polyline_2 poly = (*polyline)[i];
        if(!poly.empty())
          for(std::size_t j=0; j<poly.size()-1; ++j)
          {
            painter->drawLine(poly[j].x(), poly[j].y(), poly[j+1].x(), poly[j+1].y());
          }
      }
      painter->end();
    }
  }
  void apply_path()
  {
    update_polyline ();
    domain_rectangle = CGAL::bbox_2 (contour_2d.begin (), contour_2d.end ());
    lasso = Polygon_2 (contour_2d.begin (), contour_2d.end ());
  }

  bool is_vertex_selected (qglviewer::Vec& p)
  {
    if (domain_rectangle.xmin () < p.x &&
        p.x < domain_rectangle.xmax () &&
        domain_rectangle.ymin () < p.y &&
        p.y < domain_rectangle.ymax ())
      {
/*
 * domain_freeform.has_on_bounded_side() requires the polygon to be simple, which is never the case.
 * However, it works very well even if the polygon is not simple, so we use this instead to avoid
 * the cgal_assertion on is_simple().*/


        if (CGAL::bounded_side_2(lasso.container().begin(),
                                 lasso.container().end(),
                                 FG_Traits::Point_2(p.x, p.y),
                                 lasso.traits_member())  == CGAL::ON_BOUNDED_SIDE)
          return true;
      }
    return false;
  }
};

#endif
