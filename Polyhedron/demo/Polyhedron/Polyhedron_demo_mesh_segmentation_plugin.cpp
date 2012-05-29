#include "Polyhedron_demo_plugin_helper.h"
#include "Polyhedron_demo_plugin_interface.h"

#include "Scene_polyhedron_item.h"
#include "Scene_polyhedron_with_color_item.h"
#include "Polyhedron_type.h"
#include <CGAL/Surface_mesh_segmentation.h>
#include <QApplication>
#include <QMainWindow>
#include <QInputDialog>
#include <QTime>
#include <QAction>
#include <QDebug>

class Polyhedron_demo_mesh_segmentation_plugin : 
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_plugin_interface)

public:

  QList<QAction*> actions() const {
    return QList<QAction*>() << actionSegmentation;
  }

  void init(QMainWindow* mainWindow, Scene_interface* scene_interface) {
    this->scene = scene_interface;
    this->mw = mainWindow;
    actionSegmentation = new QAction("Mesh Segmentation", mw);
    if(actionSegmentation) {
      connect(actionSegmentation, SIGNAL(triggered()),this, SLOT(on_actionSegmentation_triggered()));
    }
  }


public slots:
  void on_actionSegmentation_triggered();
private:
  QAction*  actionSegmentation;
};

void Polyhedron_demo_mesh_segmentation_plugin::on_actionSegmentation_triggered()
{
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
  qDebug("afadsf");
  Scene_polyhedron_item* item = 
    qobject_cast<Scene_polyhedron_item*>(scene->item(index));

  if(item)
  {
    Polyhedron* pMesh = item->polyhedron();
	
	Scene_polyhedron_with_color_item* new_item = new Scene_polyhedron_with_color_item(pMesh);
	std::vector<QColor> color_vector;

	//int counter = 0;
 //   for(Polyhedron::Facet_iterator it = pMesh->facets_begin(); it != pMesh->facets_end(); ++it)
	//{
	//	int color = (counter++) % 255;
	//	color_vector.push_back(QColor(color,color,color));
	//	std::cout << color << std::endl;
	//}
	
	CGAL::Surface_mesh_segmentation<Polyhedron> segmentation(pMesh);	
  int findex=0;
	for(CGAL::Surface_mesh_segmentation<Polyhedron>::Facet_iterator facet_it = segmentation.mesh->facets_begin(); 
        facet_it != segmentation.mesh->facets_end(); ++facet_it,++findex)   
    {
      facet_it->set_patch_id(findex);
		int color = (int) (255 * segmentation.sdf_values[facet_it]);
		color_vector.push_back(QColor(color,color,color));
	}
	new_item->set_color_vector(color_vector);
	scene->itemChanged(index);
	scene->addItem(new_item);
	
  }
}

Q_EXPORT_PLUGIN2(Polyhedron_demo_mesh_segmentation_plugin, Polyhedron_demo_mesh_segmentation_plugin)

#include "Polyhedron_demo_mesh_segmentation_plugin.moc"
