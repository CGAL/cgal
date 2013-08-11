#include <QApplication>
#include <QMainWindow>
#include <QAction>
#include <QStringList>
#include <QInputDialog>

#include "Scene_polyhedron_item.h"
#include "Scene_polygon_soup_item.h"
#include "Polyhedron_type.h"

#include "Polyhedron_demo_plugin_helper.h"
#include "Polyhedron_demo_plugin_interface.h"

#include <CGAL/Mesh_3/Detect_features_in_polyhedra.h>

class Polyhedron_demo_detect_sharp_edges_plugin : 
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_plugin_interface)

public:
  void init(QMainWindow* mainWindow, Scene_interface* scene_interface) {
    this->scene = scene_interface;
    this->mw = mainWindow;
    actionSharEdges = new QAction("Detect sharp features", mw);
    actionSharEdges->setObjectName("detectSharpFeaturesAction");
    if(actionSharEdges) {
      connect(actionSharEdges, SIGNAL(triggered()),
              this, SLOT(detectSharpEdgesWithInputDialog()));
    }
  }

  bool applicable() const {
    Q_FOREACH(int index, scene->selectionIndices())
    {
      Scene_polyhedron_item* item =
        qobject_cast<Scene_polyhedron_item*>(scene->item(index));
      if (item) return true;
    }
    return false;
  }
  
  // used by Polyhedron_demo_plugin_helper
  QList<QAction*> actions() const {
    return QList<QAction*>() << actionSharEdges;
  }

public slots:
void detectSharpEdges(bool input_dialog = false, double angle = 60);
  void detectSharpEdgesWithInputDialog();

protected:
  Kernel::Vector_3 facet_normal(Polyhedron::Facet_handle f);
  bool is_sharp(Polyhedron::Halfedge_handle he);

private:
  QAction* actionSharEdges;
}; // end Polyhedron_demo_detect_sharp_edges_plugin

void Polyhedron_demo_detect_sharp_edges_plugin::detectSharpEdgesWithInputDialog()
{
  detectSharpEdges(true);
}

void Polyhedron_demo_detect_sharp_edges_plugin::detectSharpEdges(bool input_dialog,
                                                                 double angle)
{
  typedef std::pair<int,Polyhedron*> Poly_tuple;
  
  // Get selected items
  QList<Poly_tuple> polyhedrons;
  Q_FOREACH(int index, scene->selectionIndices())
  {
    Scene_polyhedron_item* item = 
      qobject_cast<Scene_polyhedron_item*>(scene->item(index));
    if(!item)
      return;
    
    Polyhedron* pMesh = item->polyhedron();
    if(!pMesh)
      return;
    
    polyhedrons << make_pair(index, pMesh);
  }

  if(input_dialog) {
    bool ok = true;
    angle = QInputDialog::getDouble(NULL, 
                                    tr("Sharp edges max angle"),
                                    tr("Angle in degrees between 0 and 180:"),
                                    angle, // value
                                    0.,          // min 
                                    180., // max
                                    2,          // decimals
                                    &ok);
    if(!ok) return;
  }
  // Detect edges
  QApplication::setOverrideCursor(Qt::WaitCursor);
  CGAL::Mesh_3::Detect_features_in_polyhedra<Polyhedron> detect_features;
  Q_FOREACH(Poly_tuple tuple, polyhedrons)
  {
    Polyhedron* pMesh = tuple.second;
    if(!pMesh) continue;

    for(Polyhedron::Edge_iterator
          eit = pMesh->edges_begin(),
          end = pMesh->edges_end(); eit != end; ++eit)
    {
      eit->set_feature_edge(false);
    }

    // Detect edges in current polyhedron
    detect_features.detect_sharp_edges(*pMesh, angle);
    
    // update scene
    scene->itemChanged(tuple.first);
  }
  
  // default cursor
  QApplication::restoreOverrideCursor();
}

Q_EXPORT_PLUGIN2(Polyhedron_demo_detect_sharp_edges_plugin, Polyhedron_demo_detect_sharp_edges_plugin)

#include "Polyhedron_demo_detect_sharp_edges_plugin.moc"
