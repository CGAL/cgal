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
      
    }
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
         
    // read text file
    std::vector< std::vector<Point_3> > polylines;
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
    dock_widget->engraveButton->setEnabled(true);
    dock_widget->visualizeButton->setEnabled(false);
  }
  
  void engrave() {
    if(!visu_item)
      return;
    
    namespace SMP = CGAL::Surface_mesh_parameterization;
    typedef CGAL::Surface_mesh_shortest_path_traits<EPICK, SMesh> SP_traits;
    typedef CGAL::Surface_mesh_shortest_path<SP_traits> Surface_mesh_shortest_path;
    typedef Surface_mesh_shortest_path::Face_location Face_location;
    typedef CGAL::AABB_face_graph_triangle_primitive<SMesh> Primitive;
    typedef CGAL::AABB_traits<Kernel, Primitive> Tree_traits;
    typedef CGAL::AABB_tree<Tree_traits> Tree;
    SMesh& sm = *qobject_cast<Scene_polyhedron_selection_item*>
        (scene->item(scene->mainSelectionIndex()))->polyhedron();
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
    
      SMesh::Property_map<SMesh::Vertex_index, EPICK::Point_3> uv_map_3 =
        sm.add_property_map<SMesh::Vertex_index, EPICK::Point_3>("v:uv3").first;
      for(SMesh::Vertex_index v : sm.vertices())
        uv_map_3[v] = Point_3(uv_map[v][0], uv_map[v]
            [1], 0);
      
/*      // debug to print 2d mesh
        std::ofstream out("param_out.off");
        out << "OFF\n" << sm.number_of_vertices() << " " << sm.number_of_faces() << " 0\n";
        for(Surface_mesh::Vertex_index v : sm.vertices())
          out << uv_map_3[v] << "\n";
        for(Surface_mesh::Face_index f : faces(sm))
        {
          Surface_mesh::Halfedge_index h = sm.halfedge(f);
          out << "3 " << (unsigned) sm.target(h) << " "
                      << (unsigned) sm.target(sm.next(h)) << " "
                      << (unsigned) sm.source(h) << "\n";
        }
      */
       
      
        // build AABB-tree for face location queries
        Tree aabb_tree(faces(sm).first, faces(sm).second, sm, uv_map_3);
      
        // compute 3D coordinates
        BOOST_FOREACH(const std::vector<Point_3>& polyline, visu_item->polylines().front());
        {
          BOOST_FOREACH(const Point_3& p, polyline)
          {
            Face_location loc = Surface_mesh_shortest_path::locate(p, aabb_tree, sm, uv_map_3);
          }
        }
    
    scene->erase(scene->item_id(visu_item));
    visu_item = NULL;
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
}; 
#include "Engrave_text_plugin.moc"
