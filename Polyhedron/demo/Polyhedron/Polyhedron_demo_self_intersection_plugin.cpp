#include <QApplication>
#include <QMessageBox>
#include <QMainWindow>
#include "Kernel_type.h"
#include "Polyhedron_type.h"
#include "Scene_polyhedron_item.h"

#include "Polyhedron_demo_plugin_helper.h"
#include "Polyhedron_demo_plugin_interface.h"

#include <CGAL/intersections.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/box_intersection_d.h>

#include <CGAL/self_intersect.h>
#include <CGAL/Make_triangle_soup.h>

typedef Kernel::Triangle_3 Triangle;

class Polyhedron_demo_self_intersection_plugin : 
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_plugin_interface)

public:
  // used by Polyhedron_demo_plugin_helper
  QStringList actionsNames() const {
    return QStringList() << "actionSelfIntersection";
  }

public slots:
  void on_actionSelfIntersection_triggered();

}; // end Polyhedron_demo_self_intersection_plugin

void Polyhedron_demo_self_intersection_plugin::on_actionSelfIntersection_triggered()
{
  Scene_interface::Item_id index = scene->mainSelectionIndex();
  
  Scene_polyhedron_item* item = 
    qobject_cast<Scene_polyhedron_item*>(scene->item(index));

  if(item)
  {
    Polyhedron* pMesh = item->polyhedron();

    // compute self-intersections
    typedef std::list<Triangle>::iterator Iterator;
    typedef CGAL::Box_intersection_d::Box_with_handle_d<double,3,Iterator> Box;
    std::list<Triangle> triangles; // intersecting triangles
    typedef std::back_insert_iterator<std::list<Triangle> > OutputIterator;
    std::cout << "Self-intersect...";
    ::self_intersect<Polyhedron,Kernel,OutputIterator>(*pMesh,std::back_inserter(triangles));
    std::cout << "ok (" << triangles.size() << " triangle(s))" << std::endl;

    // add intersecting triangles as a new polyhedron, i.e., a triangle soup.
    if(triangles.size() != 0)
    {
      Polyhedron *pSoup = new Polyhedron;
      Make_triangle_soup<Polyhedron,Kernel,Iterator> soup_builder;
      soup_builder.run(triangles.begin(),triangles.end(),*pSoup);

      Scene_polyhedron_item* new_item = new Scene_polyhedron_item(pSoup);
      new_item->setName(tr("%1 (intersecting triangles)").arg(item->name()));
      new_item->setColor(Qt::magenta);
      new_item->setRenderingMode(item->renderingMode());
      scene->addItem(new_item);
      item->setRenderingMode(Wireframe);
      scene->itemChanged(item);
    }
    else 
      QMessageBox::information(mw, tr("No self intersection"),
                               tr("The polyhedron \"%1\" does not self-intersect.").
                               arg(item->name()));
  }
}

Q_EXPORT_PLUGIN2(Polyhedron_demo_self_intersection_plugin, Polyhedron_demo_self_intersection_plugin)

#include "Polyhedron_demo_self_intersection_plugin.moc"
