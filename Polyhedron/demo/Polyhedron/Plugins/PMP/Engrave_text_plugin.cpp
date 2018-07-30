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
#include <CGAL/Polygon_mesh_processing/extrude.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Delaunay_mesh_vertex_base_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <CGAL/Delaunay_mesher_2.h>

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
        translation -= EPICK::Vector_2(0.05,0);
        visualize();
      });
      
      connect(dock_widget->t_up_pushButton, &QPushButton::clicked,
              this, [this](){
        translation -= EPICK::Vector_2(0,-0.05);
        visualize();
      });
      
      connect(dock_widget->t_down_pushButton, &QPushButton::clicked,
              this, [this](){
        translation -= EPICK::Vector_2(0,0.05);
        visualize();
      });
      
      connect(dock_widget->t_right_pushButton, &QPushButton::clicked,
              this, [this](){
        translation += EPICK::Vector_2(0.05,0);
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
    
    QList<QPolygonF> polys = path.toSubpathPolygons();
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
        visu_item->polylines.back().push_back(//p);
        Surface_mesh_shortest_path::point(loc.first, loc.second,  sm, sm.points()));
      }
    }
    visu_item->setName("Text");
    visu_item->setColor(QColor(Qt::red));
    scene->addItem(visu_item);
    dock_widget->engraveButton->setEnabled(true);
    dock_widget->visualizeButton->setEnabled(false);
  }
  
  void engrave() {
    if(!visu_item)
      return;
    typedef CGAL::Projection_traits_xy_3<EPICK>                          Gt;
    typedef CGAL::Delaunay_mesh_vertex_base_2<Gt>                        Vb;
    typedef CGAL::Delaunay_mesh_face_base_2<Gt>                          Fm;
    typedef CGAL::Triangulation_face_base_with_info_2<FaceInfo2,Gt,Fm>   Fb;
    typedef CGAL::Triangulation_data_structure_2<Vb, Fb>                TDS;
    typedef CGAL::No_intersection_tag                                   Tag;
    typedef CGAL::Constrained_Delaunay_triangulation_2<Gt, TDS, Tag>    CDT;
    typedef CGAL::Delaunay_mesh_size_criteria_2<CDT>               Criteria;
    typedef CGAL::Delaunay_mesher_2<CDT, Criteria>                   Mesher;
    
    CDT cdt;
    try{
      Q_FOREACH(const std::vector<Kernel::Point_3>& points,
                visu_item->polylines)
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
    Scene_item::Bbox bbox= visu_item->bbox();
    float diag = CGAL::sqrt(
          (bbox.xmax()-bbox.xmin())*(bbox.xmax()-bbox.xmin())
          +(bbox.ymax()-bbox.ymin())*(bbox.ymax()-bbox.ymin())
          +(bbox.zmax()-bbox.zmax()) *(bbox.zmax()-bbox.zmax())
          );
    float zmax = bbox.zmax(),
        zmin = bbox.zmin();
    //project to z=zmin for mesh_2
    BOOST_FOREACH(Scene_polylines_item::Polyline poly, visu_item->polylines)
    {
      BOOST_FOREACH(EPICK::Point_3& p, poly)
      {
        p = EPICK::Point_3(p.x(), p.y(), zmin);
      }
    }
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
    SMesh text_mesh_bottom, text_mesh_complete;
    cdt2_to_face_graph(cdt,
                       text_mesh_bottom,
                       2,
                       zmin);
    PMP::extrude_mesh(text_mesh_bottom, text_mesh_complete,EPICK::Vector_3(0,0,zmax-zmin));
    
    Scene_surface_mesh_item* text_item = new Scene_surface_mesh_item(text_mesh_complete);
    text_item->setColor(QColor(Qt::green));
    scene->addItem(text_item);
    dock_widget->engraveButton->setEnabled(false);
    dock_widget->visualizeButton->setEnabled(true);
  }
  void closure()Q_DECL_OVERRIDE
  {
    dock_widget->hide();
  }
private:
  
  struct FaceInfo2
  {
    FaceInfo2(){}
    int nesting_level;
    bool in_domain(){
      return nesting_level%2 == 1;
    }
  };
  
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
  void cdt2_to_face_graph(const CDT& cdt, TriangleMesh& tm, int constant_coordinate_index, double constant_coordinate)
  {
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
          const Kernel::Point_3& pt=fit->vertex(i)->point();
          double coords[3] = {pt[0], pt[1], pt[2]};
          coords[2]=constant_coordinate;
          it->second = add_vertex(Kernel::Point_3(coords[0],coords[1],coords[2]), tm);
        }
        vds[i]=it->second;
      }
  
      CGAL::Euler::add_face(vds, tm);
    }
  }
  
  QList<QAction*> _actions;
  EngraveWidget* dock_widget;
  Scene_polylines_item* visu_item;
  double angle;
  EPICK::Vector_2 translation;
}; 
#include "Engrave_text_plugin.moc"

