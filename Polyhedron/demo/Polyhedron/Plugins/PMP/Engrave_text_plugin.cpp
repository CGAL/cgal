//General Plugin Data
#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include <CGAL/Three/Three.h>

#include "ui_Engrave_dock_widget.h"
//Items
#include "Scene_surface_mesh_item.h"
#include "Scene_polyhedron_selection_item.h"
#include "Scene_polylines_item.h"
#include "Messages_interface.h"

#include <CGAL/Surface_mesh_parameterization/Error_code.h>
#include <CGAL/surface_mesh_parameterization.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/Polygon_mesh_processing/extrude.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Surface_mesh_shortest_path.h>
#include <CGAL/double.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>

#include <CGAL/linear_least_squares_fitting_2.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/centroid.h>

#include <CGAL/boost/graph/Face_filtered_graph.h>

#include <QPainterPath>
#include <QGraphicsScene>
#include <QGraphicsItem>
#include <QDialog>
#include <QMessageBox>

#include <CGAL/Qt/GraphicsViewNavigation.h>

using namespace CGAL::Three;
namespace SMP = CGAL::Surface_mesh_parameterization;
typedef EPICK::Point_2                                   Point_2;
typedef EPICK::Point_3                                   Point_3;

typedef boost::graph_traits<SMesh>::
edge_descriptor          edge_descriptor;
typedef boost::graph_traits<SMesh>::
halfedge_descriptor      halfedge_descriptor;
typedef boost::graph_traits<SMesh>::
vertex_descriptor        vertex_descriptor;

typedef boost::unordered_set<boost::graph_traits<SMesh>::
face_descriptor>                                         Component;

struct FaceInfo2
{
  FaceInfo2(){}
  int nesting_level;
  bool in_domain(){
    return nesting_level%2 == 1;
  }
};

template<typename PMAP,
         typename NMAP>
struct Bot
{
  Bot(NMAP nmap,
      double d,
      PMAP pmap):d(d),
    pmap(pmap),
    nmap(nmap){}
  template<typename VD, typename T>
  void operator()(const T& v1,VD v2) const
  {
    put(pmap, v2, get(pmap, v2)-d*get(nmap, v1));
  }
  double d;
  PMAP pmap;
  NMAP nmap;

};

template<typename PMAP,
         typename NMAP>
struct Top
{
  Top(NMAP nmap,
      PMAP pmap,
      double d):d(d),
    nmap(nmap),
    pmap(pmap){}

  template<typename VD, typename T>
  void operator()(const T& v1, VD v2) const
  {
    put(pmap, v2, get(pmap, v2)+d*get(nmap, v1));
  }
  double d;
  NMAP nmap;
  PMAP pmap;
};

typedef EPICK                                                        Gt;
typedef CGAL::Triangulation_vertex_base_2<Gt>                        Vb;
typedef CGAL::Triangulation_face_base_with_info_2<FaceInfo2,Gt >     Fbb;
typedef CGAL::Constrained_triangulation_face_base_2<Gt,Fbb>          Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb>                 TDS;
typedef CGAL::No_constraint_intersection_requiring_constructions_tag Tag;
typedef CGAL::Constrained_Delaunay_triangulation_2<Gt, TDS, Tag>     CDT;

//Parameterization and text displaying
class ParamItem : public QGraphicsItem
{
public :
  ParamItem(Component* component,
            const std::vector<std::vector<EPICK::Point_2> > &polylines,
            EPICK::Aff_transformation_2 transfo,
            SMesh* graph,
            QRectF brect)
    :
      QGraphicsItem(),
      bounding_rect(brect),
      component(component),
      polylines(polylines),
      graph(graph),
      transfo(transfo){}

  ~ParamItem()
  {
    delete component;
  }

  QRectF boundingRect() const
  {
    return bounding_rect;
  }

  void set_transfo(EPICK::Aff_transformation_2 t){ transfo = t;}

  void paint(QPainter *painter, const QStyleOptionGraphicsItem *, QWidget *)
  {
    QPen pen;
    QBrush brush;
    brush.setColor(QColor(100, 100, 255));
    brush.setStyle(Qt::SolidPattern);
    pen.setColor(Qt::black);
    pen.setWidth(0);
    painter->setPen(pen);
    painter->setBrush(brush);
    SMesh::Property_map<halfedge_descriptor, float> u;
    SMesh::Property_map<halfedge_descriptor, float> v;

    u = graph->add_property_map<halfedge_descriptor, float>
        ("h:u", 0.0f).first;
    v = graph->add_property_map<halfedge_descriptor, float>
        ("h:v", 0.0f).first;

    for( Component::iterator
         fi = component->begin();
         fi != component->end();
         ++fi)
    {
      boost::graph_traits<SMesh>::face_descriptor f(*fi);
      QPointF points[3];
      boost::graph_traits<SMesh>::halfedge_descriptor h = halfedge(f, *graph);;
      points[0] = QPointF(get(u, h), -get(v, h));
      h = next(halfedge(f, *graph), *graph);
      points[1] = QPointF(get(u, h), -get(v, h));
      h = next(next(halfedge(f, *graph), *graph), *graph);
      points[2] = QPointF(get(u, h), -get(v, h));
      painter->drawPolygon(points,3);
    }

    pen.setColor(Qt::red);
    pen.setWidth(0);
    painter->setPen(pen);
    for(std::size_t i =0; i<polylines.size(); ++i)
    {
      std::vector<QPointF> points;
      points.reserve(polylines[i].size());
      for(std::size_t j =0; j<polylines[i].size(); ++j)
      {
        Point_2 transfo_point = transfo.transform(polylines[i][j]);
        points.push_back(QPointF(transfo_point.x(),
                                 -transfo_point.y()));
      }
      painter->drawPolyline(points.data(), static_cast<int>(points.size()));
    }
  }

private:
  QString texMesh_name;
  QRectF bounding_rect;
  Component* component;
  const std::vector<std::vector<EPICK::Point_2> >& polylines;

  SMesh* graph;
  EPICK::Aff_transformation_2 transfo;
};

class Navigation : public CGAL::Qt::GraphicsViewNavigation
{
public:
  Navigation()
    :CGAL::Qt::GraphicsViewNavigation(),
      prev_pos(QPoint(0,0))
  { }

protected:
  bool eventFilter(QObject *obj, QEvent *ev)
  {
    QGraphicsView* v = qobject_cast<QGraphicsView*>(obj);
    if(v == NULL) {
      QWidget* viewport = qobject_cast<QWidget*>(obj);
      if(viewport == NULL) {
        return false;
      }
      v = qobject_cast<QGraphicsView*>(viewport->parent());
      if(v == NULL) {
        return false;
      }
    }
    switch(ev->type())
    {
    case QEvent::MouseMove: {
      QMouseEvent* me = static_cast<QMouseEvent*>(ev);
      if(is_dragging)
      {
        qreal dir[2] = {v->mapToScene(me->pos()).x() - prev_pos.x(),
                        v->mapToScene(me->pos()).y() - prev_pos.y()};

        v->translate(dir[0],dir[1]);
        v->update();
      }
      prev_pos = v->mapToScene(me->pos());
      break;
    }

    case QEvent::MouseButtonPress: {
      is_dragging = true;
      break;
    }
    case QEvent::MouseButtonRelease: {
      is_dragging = false;
      break;
    }
    case QEvent::Wheel: {
      QWheelEvent* event = static_cast<QWheelEvent*>(ev);
#if QT_VERSION < QT_VERSION_CHECK(5, 14, 0)
      QPoint pos = event->pos();
#else
      QPointF pos = event->position();
#endif
      QPointF old_pos = v->mapToScene(pos.x(), pos.y());
      if(event->angleDelta().y() <0)
        v->scale(1.2, 1.2);
      else
        v->scale(0.8, 0.8);
      QPointF new_pos = v->mapToScene(pos.x(), pos.y());
      QPointF delta = new_pos - old_pos;
      v->translate(delta.x(), delta.y());
      v->update();
      break;
    }

    case QEvent::MouseButtonDblClick: {
      v->fitInView(v->scene()->itemsBoundingRect(), Qt::KeepAspectRatio);
      break;
    }
    default:
      CGAL::Qt::GraphicsViewNavigation::eventFilter(obj, ev);
    }
    return false;
  }
private:
  bool is_dragging;
  QPointF prev_pos;
};


class EngraveWidget :
    public QDockWidget,
    public Ui::EngraveWidget
{
public:
  EngraveWidget(QString name, QWidget *parent)
    :QDockWidget(name,parent)
  {
    setupUi(this);
  }
};

class Q_DECL_EXPORT Engrave_text_plugin :
    public QObject,
    public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

private:
  typedef CGAL::Surface_mesh_shortest_path_traits<EPICK, SMesh> SP_traits;
  typedef CGAL::Surface_mesh_shortest_path<SP_traits> Surface_mesh_shortest_path;
  typedef Surface_mesh_shortest_path::Face_location Face_location;
  typedef CGAL::AABB_face_graph_triangle_primitive<SMesh> Primitive;
  typedef CGAL::AABB_traits<EPICK, Primitive> Tree_traits;
  typedef CGAL::AABB_tree<Tree_traits> Tree;
  typedef EPICK::Point_3 Point_3;
  Messages_interface* messages;

public :
  ~Engrave_text_plugin()
  {
    delete graphics_scene;
    delete navigation;
  }
  void init(QMainWindow*,
            CGAL::Three::Scene_interface*,
            Messages_interface* m) Q_DECL_OVERRIDE{
    //get refs
    this->scene = Three::scene();
    this->mw = Three::mainWindow();
    messages = m;

    //action
    QAction* actionFitText= new QAction("Fit Text", mw);
    connect(actionFitText, SIGNAL(triggered()),
            this, SLOT(showWidget()));
    _actions << actionFitText;
    //widget
    dock_widget = new EngraveWidget("Engraving", mw);
    dock_widget->setVisible(false); // do not show at the beginning
    addDockWidget(dock_widget);
    connect(dock_widget->visualizeButton, &QPushButton::clicked,
            this, &Engrave_text_plugin::visualize);
    connect(dock_widget->engraveButton, &QPushButton::clicked,
            this, &Engrave_text_plugin::engrave);
    connect(dock_widget->text_meshButton, &QPushButton::clicked,
            this, &Engrave_text_plugin::generateTextItem);
    connect(dock_widget->letter_checkBox, &QCheckBox::toggled,
            this, &Engrave_text_plugin::generateTextItem);

    //items
    visu_item = nullptr;
    sel_item = nullptr;
    textMesh = nullptr;
    sm = nullptr;

    //transfo
    angle = 0.0;
    scalX=1.0;
    scalY=1.0;
    translation = EPICK::Vector_2(0,0);
    pointsize = 15;
    locked = false;
    connect(dock_widget->text_prec_slider, &QSlider::valueChanged,
            this, [this](){
      pointsize = dock_widget->text_prec_slider->value();
      scene->setSelectedItem(scene->item_id(sel_item));
      visualize();
    });
    connect(dock_widget->scalX_slider, &QSlider::valueChanged,
            this, [this](){
      scalX = dock_widget->scalX_slider->value()/1000.0;
      scene->setSelectedItem(scene->item_id(sel_item));
      visualize();
    });
    connect(dock_widget->scalY_slider, &QSlider::valueChanged,
            this, [this](){
      scalY = dock_widget->scalY_slider->value()/1000.0;
      scene->setSelectedItem(scene->item_id(sel_item));
      visualize();
    });

    connect(dock_widget->bot_slider, &QSlider::valueChanged,
            this, [this](){
      if(textMesh)
        generateTextItem();
    });

    connect(dock_widget->top_slider, &QSlider::valueChanged,
            this, [this](){
      if(textMesh)
        generateTextItem();
    });

    connect(dock_widget->reset_button, &QPushButton::clicked,
            this, [this](){
      cleanup();
    });

    connect(dock_widget->t_up_pushButton, &QPushButton::clicked,
            this, [this](){
      translation += EPICK::Vector_2(0,0.005);
      scene->setSelectedItem(scene->item_id(sel_item));
      visualize();
    });

    connect(dock_widget->t_down_pushButton, &QPushButton::clicked,
            this, [this](){
      translation -= EPICK::Vector_2(0,0.005);
      scene->setSelectedItem(scene->item_id(sel_item));
      visualize();
    });

    connect(dock_widget->t_right_pushButton, &QPushButton::clicked,
            this, [this](){
      translation += EPICK::Vector_2(0.005,0);
      scene->setSelectedItem(scene->item_id(sel_item));
      visualize();
    });

    connect(dock_widget->t_left_pushButton, &QPushButton::clicked,
            this, [this](){
      translation -= EPICK::Vector_2(0.005,0);
      scene->setSelectedItem(scene->item_id(sel_item));
      visualize();
    });
    connect(dock_widget->rot_slider, &QSlider::valueChanged,
            this, [this](){
      if(!locked)
      {
        angle = dock_widget->rot_slider->value() * CGAL_PI/180.0;
        scene->setSelectedItem(scene->item_id(sel_item));
        visualize();
      }
    });
    graphics_scene = new QGraphicsScene(dock_widget);
    dock_widget->graphicsView->setScene(graphics_scene);
    dock_widget->graphicsView->setRenderHints(QPainter::Antialiasing);
    navigation = new Navigation();
    dock_widget->graphicsView->installEventFilter(navigation);
    dock_widget->graphicsView->viewport()->installEventFilter(navigation);
  }
  bool applicable(QAction*) const Q_DECL_OVERRIDE
  {
    return qobject_cast<Scene_polyhedron_selection_item*>
        (scene->item(scene->mainSelectionIndex()));
  }
  QList<QAction*> actions() const Q_DECL_OVERRIDE{
    return _actions;
  }
public Q_SLOTS:
  void showWidget()
  {
    dock_widget->setVisible(!dock_widget->isVisible());
  }

  void visualize() {
    if(!sel_item)
      sel_item =
          qobject_cast<Scene_polyhedron_selection_item*>
          (scene->item(scene->mainSelectionIndex()));
    if(!sel_item)
    {
      QMessageBox::information(mw, "Error", "No selection found.");
      return;
    }
    if(sel_item->selected_facets.empty())
    {
      QMessageBox::information(mw, "Error", "No selected facets.");
      cleanup();
      return;
    }
    if(!CGAL::is_closed(*sel_item->polyhedron()))
    {
      QMessageBox::information(mw, "Error", "The surface mesh must be closed.");
      cleanup();
      return;
    }

    connect(sel_item, &Scene_polyhedron_selection_item::aboutToBeDestroyed, this,
            [this](){
      sel_item = nullptr;
    });

    if(visu_item)
      scene->erase(scene->item_id(visu_item));
    visu_item = nullptr;

    if(!sm)
    {
      sm = new SMesh();
      sel_item->export_selected_facets_as_polyhedron(sm);
      SMesh::Halfedge_index hd =
          CGAL::Polygon_mesh_processing::longest_border(*sm).first;
      SMesh::Property_map<SMesh::Vertex_index, EPICK::Point_2> uv_map =
          sm->add_property_map<SMesh::Vertex_index, EPICK::Point_2>("v:uv").first;

      // Parameterized bool pmap
      boost::unordered_set<SMesh::Vertex_index> vs;
      SMP::internal::Bool_property_map< boost::unordered_set<SMesh::Vertex_index> > vpm(vs);

      // Parameterizer
      SMP::ARAP_parameterizer_3<SMesh> parameterizer;

      SMP::Error_code status = parameterizer.parameterize(*sm, hd, uv_map,
                                                          get(boost::vertex_index, *sm), vpm);
      if(status != SMP::OK) {
        std::cout << "Encountered a problem: " << status << std::endl;
        cleanup();
        return ;
      }

      std::cout << "Parameterized with ARAP (SM) computed." << std::endl;
      xmin = (std::numeric_limits<double>::max)();
      xmax = (std::numeric_limits<double>::min)();
      ymin = (std::numeric_limits<double>::max)();
      ymax = (std::numeric_limits<double>::min)();
      uv_map_3 =
          sm->add_property_map<SMesh::Vertex_index, Point_3>("v:uv3").first;
      for(SMesh::Vertex_index v : sm->vertices())
      {
        uv_map_3[v] = Point_3(uv_map[v][0], uv_map[v]
            [1], 0);
        if(uv_map[v][0] > xmax)
          xmax = uv_map[v][0];
        if(uv_map[v][0] < xmin)
          xmin = uv_map[v][0];

        if(uv_map[v][1] > ymax)
          ymax = uv_map[v][1];
        if(uv_map[v][1] < ymin)
          ymin = uv_map[v][1];
      }

      CGAL::linear_least_squares_fitting_2(
            uv_map.begin(),
            uv_map.end(),
            bf_line,
            CGAL::Dimension_tag<0>());

      EPICK::Vector_2 A(bf_line.to_vector()),
          B(EPICK::Point_2(0,0),
            EPICK::Point_2(1,0));
      if (A.x()<0) A=-A;
          angle = std::acos(A.x()/CGAL::sqrt(A.squared_length()));
          if ( A.y()<0 ) angle+=3*CGAL_PI/2.;
          if (angle>2*CGAL_PI) angle-=2*CGAL_PI;

      locked = true;
      dock_widget->rot_slider->setSliderPosition(angle*180.0/CGAL_PI);
      locked = false;
    }
    //create Text Polyline
    QPainterPath path;
    QFont font;
    font.setPointSize(pointsize);
    path.addText(QPoint(xmin,ymin), font, dock_widget->lineEdit->text());
    QList<QPolygonF> polys = path.toSubpathPolygons();
    polylines.clear();
    float pxmin(8000),pxmax(-8000),
        pymin(8000), pymax(-8000);

    Q_FOREACH(QPolygonF poly, polys){
      Q_FOREACH(QPointF pf, poly)
      {
        EPICK::Point_2 v = EPICK::Point_2(pf.x(),-pf.y());
        if(v.x() < pxmin)
          pxmin = v.x();
        if(v.x() > pxmax)
          pxmax = v.x();
        if(v.y() < pymin)
          pymin = v.y();
        if(v.y() > pymax)
          pymax = v.y();
      }
    }
      Q_FOREACH(QPolygonF poly, polys){
        polylines.push_back(std::vector<EPICK::Point_2>());
        Q_FOREACH(QPointF pf, poly)
        {
          EPICK::Point_2 v = EPICK::Point_2(pf.x(),-pf.y());
          polylines.back().push_back(EPICK::Point_2(v.x()*(xmax-xmin)/(pxmax-pxmin) +xmin ,
                                                    v.y()*(ymax-ymin)/(pymax-pymin)+ymin
                                                    ));
        }
      }

    // build AABB-tree for face location queries
    Tree aabb_tree(faces(*sm).first, faces(*sm).second, *sm, uv_map_3);

    visu_item = new Scene_polylines_item;
    connect(visu_item, &Scene_polylines_item::aboutToBeDestroyed, this,
            [this](){
      visu_item = nullptr;
    });

    // compute 3D coordinates
    transfo =
        EPICK::Aff_transformation_2(CGAL::TRANSLATION,
                                    EPICK::Vector_2((xmax-xmin)/2+xmin,
                                                    (ymax-ymin)/2+ymin)+ translation)
        * EPICK::Aff_transformation_2(CGAL::ROTATION,sin(angle), cos(angle))
        * EPICK::Aff_transformation_2(scalX, 0.0,0.0,scalY)
        * EPICK::Aff_transformation_2(CGAL::TRANSLATION,
                                      EPICK::Vector_2(-(xmax-xmin)/2-xmin,
                                                      -(ymax-ymin)/2-ymin));
    for(const std::vector<EPICK::Point_2>& polyline : polylines)
    {
      visu_item->polylines.push_back(std::vector<Point_3>());
      for(const EPICK::Point_2& p : polyline)
      {
        EPICK::Point_2 p_2 = transfo.transform(p);

        Face_location loc = Surface_mesh_shortest_path::locate(
              Point_3(p_2.x(), p_2.y(), 0),
              aabb_tree, *sm, uv_map_3);
        visu_item->polylines.back().push_back(
              Surface_mesh_shortest_path::point(loc.first, loc.second,  *sm, sm->points()));
      }
    }
    visu_item->setName("Text");
    visu_item->setColor(QColor(Qt::red));
    scene->addItem(visu_item);
    dock_widget->engraveButton->setEnabled(true);
    dock_widget->text_meshButton->setEnabled(true);

    if(graphics_scene->items().empty())
    {
      Component* component = new Component();
      face_iterator bfit;
      for(bfit = faces(*sm).begin();
          bfit != faces(*sm).end();
          ++bfit)
      {
        component->insert(*bfit);
      }
      SMesh::Property_map<halfedge_descriptor, float> umap;
      SMesh::Property_map<halfedge_descriptor, float> vmap;
      umap = sm->add_property_map<halfedge_descriptor, float>
        ("h:u", 0.0f).first;
      vmap = sm->add_property_map<halfedge_descriptor, float>
        ("h:v", 0.0f).first;
      SMesh::Halfedge_iterator it;
      SMesh::Property_map<SMesh::Vertex_index, EPICK::Point_2> uv_map =
          sm->property_map<SMesh::Vertex_index, EPICK::Point_2>("v:uv").first;
      for(it = sm->halfedges_begin();
          it != sm->halfedges_end();
          ++it)
      {
        halfedge_descriptor hd(*it);
        EPICK::FT u = uv_map[target(hd, *sm)].x();
        EPICK::FT v = uv_map[target(hd, *sm)].y();
        put(umap, *it, static_cast<float>(u));
        put(vmap, *it, static_cast<float>(v));
      }

      //ParamItem does not take ownership of text_mesh_bottom
      ParamItem *param_item= new ParamItem(component, polylines, transfo, sm,
                                           QRectF(QPointF(xmin, -ymax), QPointF(xmax, -ymin)));
      graphics_scene->addItem(param_item);
      dock_widget->graphicsView->fitInView(param_item->boundingRect(), Qt::KeepAspectRatio);
    }
    else
    {
      ParamItem* param_item = static_cast<ParamItem*>(graphics_scene->items().first());
      param_item->set_transfo(transfo);
      dock_widget->graphicsView->activateWindow();
      graphics_scene->update();
    }
    // dock_widget->visualizeButton->setEnabled(false);
  }

  void create_text_mesh(SMesh& text_mesh)
  {
    if(!visu_item)
      return;
    if(!sel_item)
      return;
    if(sel_item->selected_facets.empty())
      return;
    if(!CGAL::is_closed(*sel_item->polyhedron()))
      return;
    CDT cdt;
    //polylines is duplicated so the transformation is only performed once
    std::vector<std::vector<EPICK::Point_2> > local_polylines
        = polylines;
    try{
      for(std::size_t i = 0;
          i < local_polylines.size(); ++i)
      {
        std::vector<Point_2>& points = local_polylines[i];
        for(std::size_t j = 0; j< points.size(); ++j)
        {
          Point_2 &p = points[j];
          p = transfo.transform(p);
        }
        cdt.insert_constraint(points.begin(),points.end());
      }
    }catch(std::runtime_error&)
    {
      QApplication::restoreOverrideCursor();
      throw;
    }
    if (cdt.dimension()!=2){
      QApplication::restoreOverrideCursor();
      std::cout << "Triangulation is not of dimension 2" << std::endl;
      return;
    }
    CGAL::Bbox_2 bbox= CGAL::bbox_2(local_polylines.front().begin(),
                                    local_polylines.front().end(),
                                    EPICK());
    Q_FOREACH(const std::vector<EPICK::Point_2>& points,
              local_polylines)
    {
      bbox += CGAL::bbox_2(points.begin(), points.end(), EPICK());
    }
    mark_nested_domains(cdt);

    SMesh text_mesh_bottom;
    cdt2_to_face_graph(cdt,
                       text_mesh_bottom);
    typedef boost::property_map<SMesh, CGAL::vertex_point_t>::type VPMap;
    typedef SMesh::Property_map<vertex_descriptor, EPICK::Vector_3> NPMAP;
    NPMAP vnormals =
        text_mesh_bottom.add_property_map<vertex_descriptor,
        EPICK::Vector_3 >("v:normal").first;

    if(!dock_widget->letter_checkBox->isChecked())
    {
      CGAL::Polygon_mesh_processing::compute_vertex_normals(text_mesh_bottom, vnormals);
    }
    else{
      // \todo Computing normals before the final
      // mesh would be better.

      //foreach CC
      SMesh::Property_map<face_descriptor, int> fcmap =
          text_mesh_bottom.add_property_map<face_descriptor, int>("f:cc", 0).first;
      std::size_t nb_cc = PMP::connected_components(text_mesh_bottom,
                                                    fcmap);
      for(std::size_t cc = 0; cc<nb_cc; ++cc)
      {
        //compute the average normal for the cc give it to every vertex
        EPICK::Vector_3 normal(0,0,0);
        CGAL::Face_filtered_graph<SMesh> fmesh(text_mesh_bottom,
                                               static_cast<int>(cc),
                                               fcmap);
        for(vertex_descriptor vd : vertices(fmesh))
        {
          normal += CGAL::Polygon_mesh_processing::compute_vertex_normal(vd, fmesh);
        }
        normal /= CGAL::sqrt(normal.squared_length());
        for(vertex_descriptor vd : vertices(fmesh))
        {
          put(vnormals, vd, normal);
        }
      }
    }
    Bot<VPMap, NPMAP> bot(vnormals, dock_widget->bot_slider->value()/100000.0,
                          get(CGAL::vertex_point, text_mesh));
    Top<VPMap, NPMAP> top(vnormals, get(CGAL::vertex_point, text_mesh),
                          dock_widget->top_slider->value()/100000.0);
    PMP::extrude_mesh(text_mesh_bottom, text_mesh, bot, top);
  }

  void engrave() {
    QApplication::setOverrideCursor(Qt::WaitCursor);
    SMesh text_mesh_complete;
    create_text_mesh(text_mesh_complete);

    if (PMP::does_self_intersect(text_mesh_complete))
    {
      QApplication::restoreOverrideCursor();
      CGAL::Three::Three::information("Error: text mesh self-intersects!");
      return;
    }

    SMesh result;
    if(!sel_item)
      return;
    CGAL::copy_face_graph(*sel_item->polyhedron(), result);
    bool OK = PMP::corefine_and_compute_difference(result, text_mesh_complete, result);

    if (!OK)
    {
      QApplication::restoreOverrideCursor();
      CGAL::Three::Three::information("Error: the output mesh is not manifold!");
      return;
    }

    CGAL::Polygon_mesh_processing::triangulate_faces(result);
    Scene_surface_mesh_item* result_item = new Scene_surface_mesh_item(
          result);
    scene->addItem(result_item);
    graphics_scene->clear();
    cleanup();
    if(textMesh)
    {
     textMesh->setVisible(false);
    }
    QApplication::restoreOverrideCursor();
    dock_widget->engraveButton->setEnabled(false);
    dock_widget->text_meshButton->setEnabled(false);
    dock_widget->visualizeButton->setEnabled(true);
  }

  void generateTextItem()
  {
    QApplication::setOverrideCursor(Qt::WaitCursor);

    if(textMesh)
    {
      textMesh->face_graph()->clear();
      create_text_mesh(*textMesh->face_graph());
      textMesh->invalidateOpenGLBuffers();
    }
    else
    {
      SMesh text_mesh;
      create_text_mesh(text_mesh);
      textMesh = new Scene_surface_mesh_item(text_mesh);
      connect(textMesh, &Scene_surface_mesh_item::aboutToBeDestroyed,
              this, [this](){
        textMesh = nullptr;
      });
      textMesh->setName("Extruded Text");
      scene->addItem(textMesh);
    }
    QApplication::restoreOverrideCursor();
  }

  void closure()Q_DECL_OVERRIDE
  {
    dock_widget->hide();
  }

private:
  template <class CDT>
  void
  mark_domains(CDT& ct,
               typename CDT::Face_handle start,
               int index,
               std::list<typename CDT::Edge>& border )
  {
    if(start->info().nesting_level != -1){
      return;
    }
    std::list<typename CDT::Face_handle> queue;
    queue.push_back(start);
    while(! queue.empty()){
      typename CDT::Face_handle fh = queue.front();
      queue.pop_front();
      if(fh->info().nesting_level == -1){
        fh->info().nesting_level = index;
        for(int i = 0; i < 3; i++){
          typename CDT::Edge e(fh,i);
          typename CDT::Face_handle n = fh->neighbor(i);
          if(n->info().nesting_level == -1){
            if(ct.is_constrained(e)) border.push_back(e);
            else queue.push_back(n);
          }
        }
      }
    }
  }


  template <class CDT>
  void
  mark_nested_domains(CDT& cdt)
  {
    for(typename CDT::All_faces_iterator it = cdt.all_faces_begin(); it !=
        cdt.all_faces_end(); ++it){
      it->info().nesting_level = -1;
    }
    std::list<typename CDT::Edge> border;
    mark_domains(cdt, cdt.infinite_face(), 0, border);
    while(! border.empty()){
      typename CDT::Edge e = border.front();
      border.pop_front();
      typename CDT::Face_handle n = e.first->neighbor(e.second);
      if(n->info().nesting_level == -1){
        mark_domains(cdt, n, e.first->info().nesting_level+1, border);
      }
    }
  }

  template <class CDT, class TriangleMesh>
  void cdt2_to_face_graph(const CDT& cdt,
                          TriangleMesh& tm)
  {

    Tree aabb_tree(faces(*sm).first, faces(*sm).second, *sm, uv_map_3);
    typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor vertex_descriptor;

    typedef std::map<typename CDT::Vertex_handle, vertex_descriptor> Map;
    Map descriptors;
    for (typename CDT::Finite_faces_iterator fit=cdt.finite_faces_begin(),
         fit_end=cdt.finite_faces_end();
         fit!=fit_end; ++fit)
    {
      if (!fit->info().in_domain()) continue;
      std::array<vertex_descriptor,3> vds;
      for(int i=0; i<3; ++i)
      {
        typename Map::iterator it;
        bool insert_ok;
        boost::tie(it,insert_ok) =
            descriptors.insert(std::make_pair(fit->vertex(i),vertex_descriptor()));
        if (insert_ok){
          const EPICK::Point_2& pt=fit->vertex(i)->point();
          Face_location loc = Surface_mesh_shortest_path::locate(
                Point_3(pt.x(), pt.y(), 0),
                aabb_tree, *sm, uv_map_3);
          it->second = add_vertex(Surface_mesh_shortest_path::point(loc.first, loc.second,
                                                                    *sm, sm->points()), tm);
        }
        vds[i]=it->second;
      }

      CGAL::Euler::add_face(vds, tm);
    }
  }

  void cleanup()
  {
    dock_widget->scalX_slider->setValue(1000);
    dock_widget->scalY_slider->setValue(1000);
    dock_widget->rot_slider->setValue(0);
    translation = EPICK::Vector_2(0,0);
    uv_map_3.reset();
    if(graphics_scene)
      graphics_scene->clear();
    if(sel_item)
    {
      scene->erase(scene->item_id(sel_item));
      sel_item =nullptr;
    }
    if(sm)
    {
      delete sm;
      sm = nullptr;
    }
    if(visu_item)
    {
      scene->erase(scene->item_id(visu_item));
      visu_item = nullptr;
    }
  }

  QList<QAction*> _actions;
  EngraveWidget* dock_widget;
  Scene_polylines_item* visu_item;
  Scene_polyhedron_selection_item* sel_item;
  Scene_surface_mesh_item* textMesh;
  double angle;
  double scalX;
  double scalY;
  EPICK::Vector_2 translation;
  EPICK::Aff_transformation_2 transfo;
  std::vector<std::vector<EPICK::Point_2> > polylines;
  SMesh::Property_map<SMesh::Vertex_index, Point_3> uv_map_3;
  SMesh* sm;
  float xmin, xmax, ymin, ymax;
  int pointsize;
  bool locked;
  EPICK::Line_2 bf_line;
  QGraphicsScene *graphics_scene;
  Navigation* navigation;
};
#include "Engrave_text_plugin.moc"

