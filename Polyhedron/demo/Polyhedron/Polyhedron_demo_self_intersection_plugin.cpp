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

#include <CGAL/Self_intersection_polyhedron_3.h>
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

  bool applicable(QAction*) const { 
    return qobject_cast<Scene_polyhedron_item*>(scene->item(scene->mainSelectionIndex()));
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

    typedef Polyhedron::Facet_const_handle Facet_const_handle;
    std::vector<std::pair<Facet_const_handle, Facet_const_handle> > facets;
    CGAL::self_intersect<Kernel>(*pMesh, back_inserter(facets));

    std::cout << "ok (" << facets.size() << " triangle pair(s))" << std::endl;

    // add intersecting triangles as a new polyhedron, i.e., a triangle soup.
    if(!facets.empty())
    {
      Polyhedron *pSoup = new Polyhedron;

      std::vector<Facet_const_handle> facets_flat;
      facets_flat.reserve(2 * facets.size());
      for(std::size_t i = 0; i < facets.size(); ++i) 
      { facets_flat.push_back(facets[i].first); 
        facets_flat.push_back(facets[i].second); 
      }

      Make_triangle_soup<Polyhedron,Kernel,std::vector<Facet_const_handle>::iterator> soup_builder;
      soup_builder.run(facets_flat.begin(), facets_flat.end(),*pSoup);

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
