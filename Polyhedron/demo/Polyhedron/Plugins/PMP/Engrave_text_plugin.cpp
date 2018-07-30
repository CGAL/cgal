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
          [this](){
        translation -= EPICK::Vector_2(1,0);
        visualize();
      });
      
      connect(dock_widget->t_right_pushButton, &QPushButton::clicked,
          [this](){
        translation += EPICK::Vector_2(1,0);
        visualize();
      });
    }
    connect(dock_widget->r_right_pushButton, &QPushButton::clicked,
        [this](){
      angle += 0.15;
      visualize();
    });
    connect(dock_widget->r_left_pushButton, &QPushButton::clicked,
        [this](){
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
    Viewer_interface* viewer = static_cast<Viewer_interface*>(CGAL::QGLViewer::QGLViewerPool().first());
    QPainterPath path;
    QFont font;
    font.setPointSize(15);
    path.addText(QPointF(viewer->width()/2,viewer->height()/2), font, dock_widget->lineEdit->text());
    viewer->getPainter()->begin(viewer);
    viewer->getPainter()->drawPath(path);
    viewer->getPainter()->end();
    QPolygonF poly = path.toFillPolygon();
    Scene_polylines_item* new_poly = new Scene_polylines_item();
    new_poly->polylines.push_back(Scene_polylines_item::Polyline());
    Q_FOREACH(QPointF pf, poly)
    {
      new_poly->polylines.front().push_back(Point_3(pf.x(), pf.y(), 0));
    }
    scene->addItem(new_poly);
   /* if(visu_item)
      scene->erase(scene->item_id(visu_item));
    visu_item = NULL;
    typedef EPICK::Point_3 Point_3;
    // read text file
    std::ifstream text_in("/home/gimeno/Data/Polylines/text.polylines.txt");
    Scene_polylines_item::Polylines_container polylines;
    int nb_pts = 0;
    while( text_in >> nb_pts)
    {
      polylines.push_back( std::vector<Point_3>() );
      polylines.back().reserve(nb_pts);
      double x, y, z;
      for (int i=0; i<nb_pts; ++i)
      {
        text_in >> x >> y >> z;
        if (z!=0)
          std::cerr << "z-coordinate of a point is not 0: " << z << "\n";
        polylines.back().push_back( Point_3(x,y,z) );
      }
    }
    
    namespace SMP = CGAL::Surface_mesh_parameterization;
    typedef CGAL::Surface_mesh_shortest_path_traits<EPICK, SMesh> SP_traits;
    typedef CGAL::Surface_mesh_shortest_path<SP_traits> Surface_mesh_shortest_path;
    typedef Surface_mesh_shortest_path::Face_location Face_location;
    typedef CGAL::AABB_face_graph_triangle_primitive<SMesh> Primitive;
    typedef CGAL::AABB_traits<Kernel, Primitive> Tree_traits;
    typedef CGAL::AABB_tree<Tree_traits> Tree;
    typedef EPICK::Point_3 Point_3;
    Scene_polyhedron_selection_item* sel_item  = 
        qobject_cast<Scene_polyhedron_selection_item*>
        (scene->item(scene->mainSelectionIndex()));
    Point_3 basis_3[3];
    int i=0;
    for(Scene_polyhedron_selection_item::Selection_set_vertex::iterator it 
          = sel_item->selected_vertices.begin();
        it != sel_item->selected_vertices.end();
        ++it)
    {
      basis_3[(i++)%3] = sel_item->polyhedron()->point(
            *it);
    }
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
    for(SMesh::Vertex_index v : sm.vertices())
      uv_map_3[v] = Point_3(uv_map[v][0], uv_map[v]
          [1], 0);
    // build AABB-tree for face location queries
    Tree aabb_tree(faces(sm).first, faces(sm).second, sm, uv_map_3);
    
    //get 2D basis
    EPICK::Vector_2 basis_2[3];
    for(int i=0; i<3; ++i)
    {
      Surface_mesh_shortest_path::Face_location coord = Surface_mesh_shortest_path::locate(basis_3[i], aabb_tree, sm, uv_map_3);
      EPICK::Point_2 face_points[3];
      SMesh::Face_index f = SMesh::Face_index(coord.first);
      face_points[0] = uv_map[source(halfedge(f,sm),sm)];
      face_points[1] = uv_map[target(halfedge(f,sm),sm)];
      face_points[2] = uv_map[target(next(halfedge(f,sm),sm),sm)];
      
      basis_2[i] = EPICK::Vector_2(
            coord.second[0]*face_points[0].x()
          +coord.second[1]*face_points[1].x()
          +coord.second[2]*face_points[2].x(),
          coord.second[0]*face_points[0].y()
        +coord.second[1]*face_points[1].y()
        +coord.second[2]*face_points[2].y());
    }
    
    Scene_polylines_item* visu_item = new Scene_polylines_item;
    // compute 3D coordinates
    EPICK::Aff_transformation_2 transfo = 
    EPICK::Aff_transformation_2(CGAL::ROTATION,sin(angle), cos(angle)) * EPICK::Aff_transformation_2(CGAL::TRANSLATION, translation);
    BOOST_FOREACH(const std::vector<Point_3>& polyline, polylines)
    {
      visu_item->polylines.push_back(std::vector<Point_3>());
      BOOST_FOREACH(const Point_3& p, polyline)
      {
        EPICK::Vector_2 p_2 = transfo.transform(basis_2[0]+basis_2[1]*p.x()+basis_2[2]*p.y());
        
        Face_location loc = Surface_mesh_shortest_path::locate(
              Point_3(p_2.x(), p_2.y(), 0),
              aabb_tree, sm, uv_map_3);
        visu_item->polylines.back().push_back(Surface_mesh_shortest_path::point(loc.first, loc.second,  sm, sm.points()));
      }
    }
    visu_item->setName("Text");
    visu_item->setColor(QColor(Qt::red));
    scene->addItem(visu_item);
    //dock_widget->engraveButton->setEnabled(true);
    //dock_widget->visualizeButton->setEnabled(false);
    */
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

