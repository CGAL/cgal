#include <QApplication>
#include <QMainWindow>
#include <QMessageBox>
#include <QAction>
#include <QStringList>

#include "Scene_polyhedron_item.h"
#include "Polyhedron_type.h"

#include "Polyhedron_demo_plugin_helper.h"
#include "Polyhedron_demo_plugin_interface.h"

#include <CGAL/Polyhedron_kernel.h>
#include <CGAL/Gmpzf.h>
#include <CGAL/convex_hull_3.h>

#include <CGAL/Dualizer.h>
#include <CGAL/translate.h>

#include "Kernel_type.h"
typedef Kernel::Triangle_3 Triangle;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;
typedef Kernel::Plane_3 Plane;
typedef Kernel::FT FT;

class Polyhedron_demo_kernel_plugin : 
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_plugin_interface)

public:
  // used by Polyhedron_demo_plugin_helper
  QStringList actionsNames() const {
    return QStringList() << "actionKernel";
  }

  bool applicable() const { 
    return qobject_cast<Scene_polyhedron_item*>(scene->item(scene->mainSelectionIndex()));
  }

public slots:
  void on_actionKernel_triggered();

}; // end Polyhedron_demo_kernel_plugin


void Polyhedron_demo_kernel_plugin::on_actionKernel_triggered()
{
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
  
  Scene_polyhedron_item* item = 
    qobject_cast<Scene_polyhedron_item*>(scene->item(index));

  if(item)
  {
    Polyhedron* pMesh = item->polyhedron();

    typedef CGAL::Gmpzf ET; // choose exact integral type
    typedef Polyhedron_kernel<Kernel,ET> Polyhedron_kernel;

    // get triangles from polyhedron
    std::list<Triangle> triangles;
    get_triangles(*pMesh,std::back_inserter(triangles));

    // solve LP 
    std::cout << "Solve linear program...";
    Polyhedron_kernel kernel;
    if(!kernel.solve(triangles.begin(),triangles.end()))
    {
      std::cout << "done (empty kernel)" << std::endl;
      QMessageBox::information(mw, tr("Empty kernel"),
                               tr("The kernel of the polyhedron \"%1\" is empty.").
                               arg(item->name()));
      QApplication::restoreOverrideCursor();
      return;
    }
    std::cout << "done" << std::endl;

    // add kernel as new polyhedron
    Polyhedron *pKernel = new Polyhedron;

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

    // dualize and translate back to get final kernel
    Dualizer<Polyhedron,Kernel> dualizer;
    dualizer.run(convex_hull,*pKernel);
    ::translate<Polyhedron,Kernel>(*pKernel,translate);
    pKernel->inside_out();

    Scene_polyhedron_item* new_item = new Scene_polyhedron_item(pKernel);
    new_item->setName(tr("%1 (kernel)").arg(item->name()));
    new_item->setColor(Qt::magenta);
    new_item->setRenderingMode(item->renderingMode());
    scene->addItem(new_item);

    item->setRenderingMode(Wireframe);
    scene->itemChanged(item);

    QApplication::restoreOverrideCursor();
  }
}

Q_EXPORT_PLUGIN2(Polyhedron_demo_kernel_plugin, Polyhedron_demo_kernel_plugin)

#include "Polyhedron_demo_kernel_plugin.moc"
