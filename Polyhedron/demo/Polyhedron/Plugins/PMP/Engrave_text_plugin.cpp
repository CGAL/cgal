//General Plugin Data
#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include <CGAL/Three/Three.h>

#include "ui_Engrave_dock_widget.h"
//Items
#include "Scene_surface_mesh_item.h"
#include "Scene_polyhedron_selection_item.h"
#include "Scene_polylines_item.h"

//Actual code reqs
#include <CGAL/Surface_mesh_parameterization/Error_code.h>
#include <CGAL/surface_mesh_parameterization.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/Surface_mesh_shortest_path.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/centroid.h>
#include <QPainterPath>

using namespace CGAL::Three;

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
      
      //transfo
      angle = 0.0;
      translation = EPICK::Vector_2(0,0);
      connect(dock_widget->t_left_pushButton, &QPushButton::clicked,
          this, [this](){
        translation -= EPICK::Vector_2(0.1,0);
        visualize();
      });
      
      connect(dock_widget->t_right_pushButton, &QPushButton::clicked,
          this, [this](){
        translation += EPICK::Vector_2(0.1,0);
        visualize();
      });
    }
    connect(dock_widget->r_right_pushButton, &QPushButton::clicked,
        this, [this](){
      angle += 0.15;
      visualize();
    });
    connect(dock_widget->r_left_pushButton, &QPushButton::clicked,
        this, [this](){
      angle -= 0.15;
      visualize();
    });
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
    Scene_polyhedron_selection_item* sel_item  = 
        qobject_cast<Scene_polyhedron_selection_item*>
        (scene->item(scene->mainSelectionIndex()));
    if(!sel_item)
      return;
    if(sel_item->selected_facets.empty())
      return;
    if(visu_item)
      scene->erase(scene->item_id(visu_item));
    visu_item = NULL;
    
    namespace SMP = CGAL::Surface_mesh_parameterization;
    typedef CGAL::Surface_mesh_shortest_path_traits<EPICK, SMesh> SP_traits;
    typedef CGAL::Surface_mesh_shortest_path<SP_traits> Surface_mesh_shortest_path;
    typedef Surface_mesh_shortest_path::Face_location Face_location;
    typedef CGAL::AABB_face_graph_triangle_primitive<SMesh> Primitive;
    typedef CGAL::AABB_traits<Kernel, Primitive> Tree_traits;
    typedef CGAL::AABB_tree<Tree_traits> Tree;
    typedef EPICK::Point_3 Point_3;
    SMesh sm;
    sel_item->export_selected_facets_as_polyhedron(&sm);
    SMesh::Halfedge_index hd =
        CGAL::Polygon_mesh_processing::longest_border(sm).first;
    SMesh::Property_map<SMesh::Vertex_index, EPICK::Point_2> uv_map =
        sm.add_property_map<SMesh::Vertex_index, EPICK::Point_2>("v:uv").first;
    
    // Parameterized bool pmap
    boost::unordered_set<SMesh::Vertex_index> vs;
    SMP::internal::Bool_property_map< boost::unordered_set<SMesh::Vertex_index> > vpm(vs);
    
    // Parameterizer
    SMP::ARAP_parameterizer_3<SMesh> parameterizer;
    
    SMP::Error_code status = parameterizer.parameterize(sm, hd, uv_map, get(boost::vertex_index, sm), vpm);
    if(status != SMP::OK) {
      std::cout << "Encountered a problem: " << status << std::endl;
      return ;
    }
    
    std::cout << "Parameterized with ARAP (SM)!" << std::endl;
    SMesh::Property_map<SMesh::Vertex_index, Point_3> uv_map_3 =
        sm.add_property_map<SMesh::Vertex_index, Point_3>("v:uv3").first;
    float xmin(std::numeric_limits<double>::max()), xmax(std::numeric_limits<double>::min()), 
          ymin(std::numeric_limits<double>::max()), ymax(std::numeric_limits<double>::min());
    for(SMesh::Vertex_index v : sm.vertices())
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
      out << "OFF\n" << sm.number_of_vertices() << " " << sm.number_of_faces() << " 0\n";
      for(SMesh::Vertex_index v : sm.vertices())
        out << uv_map_3[v] << "\n";
      for(SMesh::Face_index f : faces(sm))
      {
        SMesh::Halfedge_index h = sm.halfedge(f);
        out << "3 " << (unsigned) sm.target(h) << " "
                    << (unsigned) sm.target(sm.next(h)) << " "
                    << (unsigned) sm.source(h) << "\n";
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
    
    QList<QPolygonF> polys = path.toFillPolygons();
    //Scene_polylines_item* new_poly = new Scene_polylines_item();
    Scene_polylines_item::Polylines_container polylines;
    QFontMetrics fm(font);
    int width=fm.width(dock_widget->lineEdit->text());
    int height=fm.height();
    std::cout<<"width = "<<width<<", height= "<<height<<std::endl;
    Q_FOREACH(QPolygonF poly, polys){
      polylines.push_back(Scene_polylines_item::Polyline());
      Q_FOREACH(QPointF pf, poly)
      {
        EPICK::Point_2 v = EPICK::Point_2(pf.x(),-pf.y());
        polylines.back().push_back(Point_3(v.x()/2*(xmax-xmin)/width + (xmin) ,
                                           v.y()/2*(ymax-ymin)/height+ (ymin) ,
                                           0));
      }
    }
    
    
    
    
    
    // build AABB-tree for face location queries
    Tree aabb_tree(faces(sm).first, faces(sm).second, sm, uv_map_3);
    
    visu_item = new Scene_polylines_item;
    // compute 3D coordinates
    EPICK::Aff_transformation_2 transfo = 
    EPICK::Aff_transformation_2(CGAL::ROTATION,sin(angle), cos(angle)) * EPICK::Aff_transformation_2(CGAL::TRANSLATION, translation);
    BOOST_FOREACH(const std::vector<Point_3>& polyline, polylines)
    {
      visu_item->polylines.push_back(std::vector<Point_3>());
      BOOST_FOREACH(const Point_3& p, polyline)
      {
        EPICK::Point_2 p_2 = transfo.transform(EPICK::Point_2(p.x(), p.y()));
        
        Face_location loc = Surface_mesh_shortest_path::locate(
              Point_3(p_2.x(), p_2.y(), 0),
              aabb_tree, sm, uv_map_3);
        visu_item->polylines.back().push_back(//Point_3(p_2.x(), p_2.y(), 0));
              Surface_mesh_shortest_path::point(loc.first, loc.second,  sm, sm.points()));
      }
    }
    visu_item->setName("Text");
    visu_item->setColor(QColor(Qt::red));
    scene->addItem(visu_item);
    //dock_widget->engraveButton->setEnabled(true);
    //dock_widget->visualizeButton->setEnabled(false);
  }
  
  void engrave() {
    if(!visu_item)
      return;
    
    
    dock_widget->engraveButton->setEnabled(false);
    dock_widget->visualizeButton->setEnabled(true);
  }
  void closure()Q_DECL_OVERRIDE
  {
    dock_widget->hide();
  }
private:
  QList<QAction*> _actions;
  EngraveWidget* dock_widget;
  Scene_polylines_item* visu_item;
  double angle;
  EPICK::Vector_2 translation;
}; 
#include "Engrave_text_plugin.moc"

