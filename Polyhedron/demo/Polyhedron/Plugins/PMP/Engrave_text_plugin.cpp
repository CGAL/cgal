//General Plugin Data
#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include <CGAL/Three/Three.h>

#include "ui_Engrave_dock_widget.h"
//Items
#include "Scene_surface_mesh_item.h"
#include "Scene_polyhedron_selection_item.h"
#include "Scene_polylines_item.h"
#include "Nef_type.h"

//Actual code reqs
#include <CGAL/boost/graph/convert_nef_polyhedron_to_polygon_mesh.h>

#include <CGAL/Surface_mesh_parameterization/Error_code.h>
#include <CGAL/surface_mesh_parameterization.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/Polygon_mesh_processing/extrude.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Surface_mesh_shortest_path.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Delaunay_mesh_vertex_base_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <CGAL/Delaunay_mesher_2.h>

#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/centroid.h>

#include <QPainterPath>
#include <QGraphicsScene>
#include <QGraphicsItem>

#include <CGAL/Qt/GraphicsViewNavigation.h>

using namespace CGAL::Three;
namespace SMP = CGAL::Surface_mesh_parameterization;
typedef EPICK::Point_2                                            Point_2;

typedef boost::graph_traits<SMesh>::
edge_descriptor          edge_descriptor;
typedef boost::graph_traits<SMesh>::
halfedge_descriptor      halfedge_descriptor;
typedef boost::graph_traits<SMesh>::
vertex_descriptor        vertex_descriptor;

typedef boost::unordered_set<boost::graph_traits<SMesh>::
face_descriptor>                                                    Component;

struct FaceInfo2
{
  FaceInfo2(){}
  int nesting_level;
  bool in_domain(){
    return nesting_level%2 == 1;
  }
};

typedef EPICK                                                        Gt;
typedef CGAL::Delaunay_mesh_vertex_base_2<Gt>                        Vb;
typedef CGAL::Delaunay_mesh_face_base_2<Gt>                          Fm;
typedef CGAL::Triangulation_face_base_with_info_2<FaceInfo2,Gt,Fm>   Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb>                TDS;
typedef CGAL::No_intersection_tag                                   Tag;
typedef CGAL::Constrained_Delaunay_triangulation_2<Gt, TDS, Tag>    CDT;
typedef CGAL::Delaunay_mesh_size_criteria_2<CDT>               Criteria;
typedef CGAL::Delaunay_mesher_2<CDT, Criteria>                   Mesher;

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
    delete graph;
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
    SMesh::Property_map<halfedge_descriptor,std::pair<float, float> > uv;
    uv = graph->add_property_map<halfedge_descriptor,std::pair<float, float> >("h:uv",std::make_pair(0.0f,0.0f)).first;
    for( Component::iterator
         fi = component->begin();
         fi != component->end();
         ++fi)
    {
      boost::graph_traits<SMesh>::face_descriptor f(*fi);
      QPointF points[3];
      boost::graph_traits<SMesh>::halfedge_descriptor h = halfedge(f, *graph);;
      points[0] = QPointF(get(uv, h).first, get(uv, h).second);
      h = next(halfedge(f, *graph), *graph);
      points[1] = QPointF(get(uv, h).first, get(uv, h).second);
      h = next(next(halfedge(f, *graph), *graph), *graph);
      points[2] = QPointF(get(uv, h).first, get(uv, h).second);
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
                                 transfo_point.y()));
      }
      painter->drawPolyline(points.data(), points.size());
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
      QPointF old_pos = v->mapToScene(event->pos());
      if(event->delta() <0)
        v->scale(1.2, 1.2);
      else
        v->scale(0.8, 0.8);
      QPointF new_pos = v->mapToScene(event->pos());
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
  typedef CGAL::AABB_traits<Kernel, Primitive> Tree_traits;
  typedef CGAL::AABB_tree<Tree_traits> Tree;
  typedef EPICK::Point_3 Point_3;
  
public :
  
  void init(QMainWindow* , CGAL::Three::Scene_interface* , Messages_interface*) Q_DECL_OVERRIDE{
    //get refs
    this->scene = Three::scene();
    this->mw = Three::mainWindow();
    
    //action
    QAction* actionFitText= new QAction("Fit Text", mw);
    if(actionFitText) {
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
      
      //items
      visu_item = NULL;
      sel_item = NULL;
      
      //transfo
      angle = 0.0;
      translation = EPICK::Vector_2(0,0);
      connect(dock_widget->t_left_pushButton, &QPushButton::clicked,
              this, [this](){
        translation -= EPICK::Vector_2(0.05,0);
        scene->setSelectedItem(scene->item_id(sel_item));
        visualize();
      });
      
      connect(dock_widget->t_up_pushButton, &QPushButton::clicked,
              this, [this](){
        translation -= EPICK::Vector_2(0,0.05);
        scene->setSelectedItem(scene->item_id(sel_item));
        visualize();
      });
      
      connect(dock_widget->t_down_pushButton, &QPushButton::clicked,
              this, [this](){
        translation += EPICK::Vector_2(0,0.05);
        scene->setSelectedItem(scene->item_id(sel_item));
        visualize();
      });
      
      connect(dock_widget->t_right_pushButton, &QPushButton::clicked,
              this, [this](){
        translation += EPICK::Vector_2(0.05,0);
        scene->setSelectedItem(scene->item_id(sel_item));
        visualize();
      });
    }
    connect(dock_widget->r_right_pushButton, &QPushButton::clicked,
            this, [this](){
      angle += 0.15;
      scene->setSelectedItem(scene->item_id(sel_item));
      visualize();
    });
    connect(dock_widget->r_left_pushButton, &QPushButton::clicked,
            this, [this](){
      angle -= 0.15;
      scene->setSelectedItem(scene->item_id(sel_item));
      visualize();
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
      return;
    if(sel_item->selected_facets.empty())
      return;
    if(!CGAL::is_closed(*sel_item->polyhedron()))
      return;
    if(visu_item)
      scene->erase(scene->item_id(visu_item));
    visu_item = NULL;
    
    
    SMesh *sm = new SMesh();
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
    
    SMP::Error_code status = parameterizer.parameterize(*sm, hd, uv_map, get(boost::vertex_index, *sm), vpm);
    if(status != SMP::OK) {
      std::cout << "Encountered a problem: " << status << std::endl;
      return ;
    }
    
    std::cout << "Parameterized with ARAP (SM)!" << std::endl;
    uv_map_3 =
        sm->add_property_map<SMesh::Vertex_index, Point_3>("v:uv3").first;
    float xmin(std::numeric_limits<double>::max()), xmax(std::numeric_limits<double>::min()), 
        ymin(std::numeric_limits<double>::max()), ymax(std::numeric_limits<double>::min());
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
    std::cout<<"xmax = "<<xmax<<", ymax = "<<ymax<<std::endl;
    
    std::ofstream out("param_out.off");
    out << "OFF\n" << sm->number_of_vertices() << " " << sm->number_of_faces() << " 0\n";
    for(SMesh::Vertex_index v : sm->vertices())
      out << uv_map_3[v] << "\n";
    for(SMesh::Face_index f : faces(*sm))
    {
      SMesh::Halfedge_index h = sm->halfedge(f);
      out << "3 " << (unsigned) sm->target(h) << " "
          << (unsigned) sm->target(sm->next(h)) << " "
          << (unsigned) sm->source(h) << "\n";
    }
    
    
    //create Text Polyline
    typedef EPICK::Point_3 Point_3;
    Viewer_interface* viewer = static_cast<Viewer_interface*>(CGAL::QGLViewer::QGLViewerPool().first());
    QPainterPath path;
    QFont font;
    font.setPointSize(15);
    
    path.addText(QPoint(xmin,-ymin), font, dock_widget->lineEdit->text());
    viewer->getPainter()->begin(viewer);
    viewer->getPainter()->drawPath(path);
    viewer->getPainter()->end();
    QTransform trans;
    
    QList<QPolygonF> polys = path.toSubpathPolygons();
    QFontMetrics fm(font);
    int width=fm.width(dock_widget->lineEdit->text());
    int height=fm.height();
    std::cout<<"width = "<<width<<", height= "<<height<<std::endl;
    Q_FOREACH(QPolygonF poly, polys){
      polylines.push_back(std::vector<EPICK::Point_2>());
      Q_FOREACH(QPointF pf, poly)
      {
        EPICK::Point_2 v = EPICK::Point_2(pf.x(),-pf.y());
        polylines.back().push_back(EPICK::Point_2(v.x()*(xmax-xmin)/width +xmin ,
                                                  v.y()*(ymax-ymin)/height+ymin
                                                  ));
      }
    }
    
    
    
    // build AABB-tree for face location queries
    Tree aabb_tree(faces(*sm).first, faces(*sm).second, *sm, uv_map_3);
    
    visu_item = new Scene_polylines_item;
    // compute 3D coordinates
    EPICK::Aff_transformation_2 transfo = 
        EPICK::Aff_transformation_2(CGAL::TRANSLATION, EPICK::Vector_2((-(width)/2), ((height)/2))) 
        * EPICK::Aff_transformation_2(CGAL::ROTATION,sin(angle), cos(angle)) 
        * EPICK::Aff_transformation_2(CGAL::TRANSLATION, EPICK::Vector_2(((width)/2), -((height)/2)) + translation);
    BOOST_FOREACH(const std::vector<EPICK::Point_2>& polyline, polylines)
    {
      visu_item->polylines.push_back(std::vector<Point_3>());
      BOOST_FOREACH(const EPICK::Point_2& p, polyline)
      {
        EPICK::Point_2 p_2 = transfo.transform(p);
        
        Face_location loc = Surface_mesh_shortest_path::locate(
              Point_3(p_2.x(), p_2.y(), 0),
              aabb_tree, *sm, uv_map_3);
        visu_item->polylines.back().push_back(//p);
                                              Surface_mesh_shortest_path::point(loc.first, loc.second,  *sm, sm->points()));
      }
    }
    visu_item->setName("Text");
    visu_item->setColor(QColor(Qt::red));
    scene->addItem(visu_item);
    dock_widget->engraveButton->setEnabled(true);
    
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
      SMesh::Property_map<halfedge_descriptor,std::pair<float, float> > uv;
      uv = sm->add_property_map<halfedge_descriptor,std::pair<float, float> >(
            "h:uv",std::make_pair(0.0f,0.0f)).first;
      SMesh::Halfedge_iterator it;
      for(it = sm->halfedges_begin();
          it != sm->halfedges_end();
          ++it)
      {
        halfedge_descriptor hd(*it);
        EPICK::FT u = uv_map[target(hd, *sm)].x();
        EPICK::FT v = uv_map[target(hd, *sm)].y();
        put(uv, *it, std::make_pair(static_cast<float>(u),static_cast<float>(v)));
      }
      
      //ParamItem takes ownership of text_mesh_bottom
      ParamItem *param_item= new ParamItem(component, polylines, transfo, sm,
                                           QRectF(QPointF(xmin, ymin), QPointF(xmax, ymax))); 
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
  
  void engrave() {
    if(!visu_item)
      return;
    Scene_polyhedron_selection_item* sel_item  = 
        qobject_cast<Scene_polyhedron_selection_item*>
        (scene->item(scene->mainSelectionIndex()));
    if(!sel_item)
      return;
    if(sel_item->selected_facets.empty())
      return;
    if(!CGAL::is_closed(*sel_item->polyhedron()))
      return;
    
    /*   CDT cdt;
    try{
      Q_FOREACH(const std::vector<Kernel::Point_2>& points,
                polylines)
        cdt.insert_constraint(points.begin(),points.end());
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
    CGAL::Bbox_2 bbox= CGAL::bbox_2(polylines.front().begin(), polylines.front().end(), EPICK());
    Q_FOREACH(const std::vector<Kernel::Point_2>& points,
              polylines)
    {
      bbox += CGAL::bbox_2(points.begin(), points.end(), EPICK());
    }
    float diag = CGAL::sqrt(
          (bbox.xmax()-bbox.xmin())*(bbox.xmax()-bbox.xmin())
          +(bbox.ymax()-bbox.ymin())*(bbox.ymax()-bbox.ymin())
          );
    // start by marking the domain to mesh
    Criteria criteria(0.125, 0.05 * diag);
    Mesher mesher(cdt, criteria);
    
    mark_nested_domains(cdt);
    for(typename CDT::All_faces_iterator fit=cdt.all_faces_begin(),
        fit_end=cdt.all_faces_end();
        fit!=fit_end;++fit)
    {
      fit->set_in_domain(fit->info().in_domain());
    }
    mesher.init(true);
    
    
    mesher.refine_mesh();
    SMesh* text_mesh_bottom = new SMesh();
    cdt2_to_face_graph(cdt,
                       sm,
                       *text_mesh_bottom);*/
    //    PMP::extrude_mesh(text_mesh_bottom, text_mesh_complete,EPICK::Vector_3(0,0,0.3));
    
    /*CGAL::Surface_mesh<Exact_Kernel::Point_3> exact_text,
        exact_target;
    CGAL::copy_face_graph(text_mesh_complete, exact_text);
    CGAL::copy_face_graph(*sel_item->polyhedron(), exact_target);
    Nef_polyhedron nef_text(exact_text);
    Nef_polyhedron nef_target(exact_target);
    Nef_polyhedron new_nef = nef_target - nef_text;
    SMesh result;
    CGAL::convert_nef_polyhedron_to_polygon_mesh(new_nef, result);
    CGAL::Polygon_mesh_processing::triangulate_faces(result);*/
    //Scene_surface_mesh_item* result_item = new Scene_surface_mesh_item(text_mesh_bottom);
    //result);
    //scene->addItem(result_item);
    dock_widget->engraveButton->setEnabled(false);
    dock_widget->visualizeButton->setEnabled(true);
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
    for(typename CDT::All_faces_iterator it = cdt.all_faces_begin(); it != cdt.all_faces_end(); ++it){
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
                          SMesh sm,
                          TriangleMesh& tm)
  {
    
    Tree aabb_tree(faces(sm).first, faces(sm).second, sm, uv_map_3);
    typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor vertex_descriptor;
    
    typedef std::map<typename CDT::Vertex_handle, vertex_descriptor> Map;
    Map descriptors;
    for (typename CDT::Finite_faces_iterator fit=cdt.finite_faces_begin(),
         fit_end=cdt.finite_faces_end();
         fit!=fit_end; ++fit)
    {
      if (!fit->is_in_domain()) continue;
      CGAL::cpp11::array<vertex_descriptor,3> vds;
      for(int i=0; i<3; ++i)
      {
        typename Map::iterator it;
        bool insert_ok;
        boost::tie(it,insert_ok) =
            descriptors.insert(std::make_pair(fit->vertex(i),vertex_descriptor()));
        if (insert_ok){
          const Kernel::Point_2& pt=fit->vertex(i)->point();
          Face_location loc = Surface_mesh_shortest_path::locate(
                Point_3(pt.x(), pt.y(), 0),
                aabb_tree, sm, uv_map_3);
          it->second = add_vertex(Surface_mesh_shortest_path::point(loc.first, loc.second,  sm, sm.points()), tm);
        }
        vds[i]=it->second;
      }
      
      CGAL::Euler::add_face(vds, tm);
    }
  }
  
  QList<QAction*> _actions;
  EngraveWidget* dock_widget;
  Scene_polylines_item* visu_item;
  Scene_polyhedron_selection_item* sel_item;
  double angle;
  EPICK::Vector_2 translation;
  std::vector<std::vector<EPICK::Point_2> > polylines;
  SMesh::Property_map<SMesh::Vertex_index, Point_3> uv_map_3;
  QGraphicsScene *graphics_scene;
  Navigation* navigation;
}; 
#include "Engrave_text_plugin.moc"

