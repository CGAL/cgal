#include <fstream>
#include <QApplication>
#include <QMainWindow>
#include <QMessageBox>
#include <QAction>
#include <QStringList>

#include "Scene_polyhedron_item.h"
#include "Scene_surface_mesh_item.h"
#include "Polyhedron_type.h"

#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>
#include <CGAL/Exact_rational.h>
#include <CGAL/Polyhedron_kernel.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Dualizer.h>
#include <CGAL/translate.h>
#include <CGAL/boost/graph/Dual.h>
#include <CGAL/boost/graph/copy_face_graph.h>

#include "Kernel_type.h"
typedef Kernel::Triangle_3 Triangle;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;
typedef Kernel::Plane_3 Plane;
typedef Kernel::FT FT;
using namespace CGAL::Three;
class Polyhedron_demo_kernel_plugin : 
  public QObject,
  public Polyhedron_demo_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

public:
  QList<QAction*> actions() const {
    return _actions;
  }

  bool applicable(QAction*) const { 
    return qobject_cast<Scene_polyhedron_item*>(scene->item(scene->mainSelectionIndex())) ||
      qobject_cast<Scene_surface_mesh_item*>(scene->item(scene->mainSelectionIndex()));
  }

   void init(QMainWindow* mainWindow, CGAL::Three::Scene_interface* scene_interface, Messages_interface*)
   {
     scene = scene_interface;
     mw = mainWindow;
     QAction* actionKernel = new QAction(tr("Kernel"), mainWindow);
     connect(actionKernel, SIGNAL(triggered()),
             this, SLOT(on_actionKernel_triggered()));
     _actions << actionKernel;
   }
public Q_SLOTS:
  void on_actionKernel_triggered();
private:
  QList<QAction*> _actions;
  CGAL::Three::Scene_interface* scene;
  QMainWindow* mw;

}; // end Polyhedron_demo_kernel_plugin


void Polyhedron_demo_kernel_plugin::on_actionKernel_triggered()
{
  const CGAL::Three::Scene_interface::Item_id index = scene->mainSelectionIndex();
  
  Scene_polyhedron_item* item = 
    qobject_cast<Scene_polyhedron_item*>(scene->item(index));

  Scene_surface_mesh_item* sm_item = 
    qobject_cast<Scene_surface_mesh_item*>(scene->item(index));

  if(item || sm_item)
  {
    QApplication::setOverrideCursor(Qt::WaitCursor);
    Polyhedron* pMesh = (item)?item->polyhedron():NULL;
    SMesh* sMesh = (sm_item)?sm_item->polyhedron():NULL;

    typedef CGAL::Exact_rational ET;
    typedef Polyhedron_kernel<Kernel,ET> Polyhedron_kernel;

    // get triangles from polyhedron
    std::list<Triangle> triangles;
    if(pMesh){
      get_triangles(*pMesh,std::back_inserter(triangles));
    } else {
      get_triangles(*sMesh,std::back_inserter(triangles));
    }
    // solve LP 
    std::cout << "Solve linear program..." << triangles.size();
    Polyhedron_kernel kernel;
    if(!kernel.solve(triangles.begin(),triangles.end()))
    {
      std::cout << "done (empty kernel)" << std::endl;
      QApplication::restoreOverrideCursor();
      QMessageBox::information(mw, tr("Empty Kernel"),
                               tr("The kernel of the polyhedron \"%1\" is empty.").
                               arg((item)?item->name():sm_item->name()));

      return;
    }
    std::cout << "done" << std::endl;

    // get inside point
    Point inside_point = kernel.inside_point();
    Vector translate = inside_point - CGAL::ORIGIN;

    // compute dual of translated polyhedron w.r.t. inside point.
    std::cout << "Compute dual of translated polyhedron...";
    std::list<Point> dual_points;
    std::list<Triangle>::iterator it;
    for(it = triangles.begin();
      it != triangles.end();
      it++)
    {
      const Triangle& triangle = *it;
      const Point p0 = triangle[0] - translate;
      const Point p1 = triangle[1] - translate;
      const Point p2 = triangle[2] - translate;
      Plane plane(p0,p1,p2); 
      Vector normal = plane.orthogonal_vector();
      normal = normal / std::sqrt(normal*normal);
      // compute distance to origin (do not use plane.d())
      FT distance_to_origin = std::sqrt(CGAL::squared_distance(Point(CGAL::ORIGIN),plane));
      Point dual_point = CGAL::ORIGIN + normal / distance_to_origin;
      dual_points.push_back(dual_point);
    }
    std::cout << "ok" << std::endl;

    // compute convex hull in dual space
    std::cout << "convex hull in dual space...";
    Polyhedron convex_hull;
    CGAL::convex_hull_3(dual_points.begin(),dual_points.end(),convex_hull);
    std::cout << "ok" << std::endl;

    if(pMesh){
      // add kernel as new polyhedron
      Polyhedron *pKernel = new Polyhedron;
      
      // dualize and translate back to get final kernel
      Dualizer<Polyhedron,Kernel> dualizer;
      dualizer.run(convex_hull,*pKernel);
      ::translate<Polyhedron,Kernel>(*pKernel,translate);
      pKernel->inside_out();
      
      Scene_polyhedron_item* new_item = new Scene_polyhedron_item(pKernel);
      new_item->setName(tr("%1 (kernel)").arg(item->name()));
      new_item->setColor(Qt::magenta);
      new_item->setRenderingMode(item->renderingMode());
      item->setRenderingMode(Wireframe);
      
      scene->addItem(new_item);
      scene->itemChanged(item);
    } else { // sMesh
      // add kernel as new polyhedron
      typedef SMesh SMesh;
      SMesh* pKernel = new SMesh;
      
      typedef CGAL::Dual<Polyhedron> Dual;
      typedef boost::graph_traits<Dual>::vertex_descriptor dual_vertex_descriptor;

      std::ofstream out("primal.off");
      out << convex_hull << std::endl;
      out.close();
      Dual dual(convex_hull);
      
      std::map<dual_vertex_descriptor,Point> vpm; 
      BOOST_FOREACH(boost::graph_traits<Polyhedron>::face_descriptor fd, faces(convex_hull)){  
        Point points[3];
        int i = 0;
        BOOST_FOREACH(boost::graph_traits<Polyhedron>::vertex_descriptor vd, 
                      vertices_around_face(halfedge(fd,convex_hull),convex_hull)){
          points[i++] = get(CGAL::vertex_point,convex_hull,vd);
        }
        Plane plane(points[0],points[1],points[2]); 
        FT sqd = CGAL::squared_distance(Point(CGAL::ORIGIN),plane);
        FT distance_to_origin = std::sqrt(sqd);
        Vector normal = plane.orthogonal_vector();
        normal = normal / std::sqrt(normal * normal);
        vpm[fd] = CGAL::ORIGIN + normal / distance_to_origin  + translate;
      }
      CGAL::copy_face_graph(dual, *pKernel,
                            CGAL::Emptyset_iterator(),
                            CGAL::Emptyset_iterator(),
                            CGAL::Emptyset_iterator(),
                            boost::make_assoc_property_map(vpm));

    
      // pKernel->inside_out();
      
      Scene_surface_mesh_item* new_item = new Scene_surface_mesh_item(pKernel);
      new_item->setName(tr("%1 (kernel)").arg(sm_item->name()));
      new_item->setColor(Qt::magenta);
      new_item->setRenderingMode(sm_item->renderingMode());
      sm_item->setRenderingMode(Wireframe);
      
      scene->addItem(new_item);
      scene->itemChanged(sm_item);
    }
    QApplication::restoreOverrideCursor();
  }
}

#include "Kernel_plugin.moc"
