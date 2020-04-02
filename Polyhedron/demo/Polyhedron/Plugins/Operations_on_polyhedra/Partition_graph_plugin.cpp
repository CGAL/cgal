#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>

#include "ui_PartitionDialog.h"
#include "Color_map.h"

#include "Scene_surface_mesh_item.h"

#include <CGAL/boost/graph/METIS/partition_graph.h>
#include <CGAL/boost/graph/METIS/partition_dual_graph.h>

#include <QString>
#include <QAction>
#include <QMenu>
#include <QMainWindow>
#include <QApplication>
#include <QElapsedTimer>
#include <QMessageBox>

typedef Scene_surface_mesh_item Scene_facegraph_item;

typedef Scene_facegraph_item::Face_graph FaceGraph;
class PartitionDialog :
    public QDialog,
    public Ui::PartitionDialog
{
  Q_OBJECT
public:
  PartitionDialog(QWidget* =0)
  {
    setupUi(this);
  }
};

using namespace CGAL::Three;
class Partition_graph_plugin :
    public QObject,
    public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

public:

  QList<QAction*> actions() const {
    return QList<QAction*>() << actionNodalPartition
                             << actionDualPartition;;
  }

  bool applicable(QAction*) const {
    return qobject_cast<Scene_facegraph_item*>(scene->item(scene->mainSelectionIndex()));
  }

  void init(QMainWindow* _mw, CGAL::Three::Scene_interface* scene_interface, Messages_interface*) {
    mw = _mw;
    this->scene = scene_interface;


    actionNodalPartition = new QAction(
                tr("Create a Nodal Graph Based Partition")
          , mw);
    if(actionNodalPartition) {
      connect(actionNodalPartition, SIGNAL(triggered()),this, SLOT(create_nodal_partition()));
    }

    actionDualPartition = new QAction(
                tr("Create a Dual Graph Based Partition")
          , mw);
    if(actionDualPartition) {
      connect(actionDualPartition, SIGNAL(triggered()),this, SLOT(create_dual_partition()));
    }
  }

private:
  QAction*  actionNodalPartition;
  QAction*  actionDualPartition;
  CGAL::Three::Scene_interface* scene;
  enum PARTITION_TYPE{
    NODAL=0,
    DUAL};

  void create_partition(PARTITION_TYPE type)
  {
    Scene_facegraph_item* item =
      qobject_cast<Scene_facegraph_item*>(scene->item(scene->mainSelectionIndex()));
    if(!item)
      return;
    if(!(CGAL::is_triangle_mesh(*item->face_graph())
         && is_valid(*item->face_graph())))
      return;
    PartitionDialog *dialog = new PartitionDialog();
       //opens the dialog
       if(!dialog->exec())
         return;
       int nparts = dialog->nparts_spinBox->value();
    QApplication::setOverrideCursor(Qt::WaitCursor);
    item->face_graph()->collect_garbage();
    item->color_vector().clear();
    if(!item->hasPatchIds()){
      item->setItemIsMulticolor(true);
      item->computeItemColorVectorAutomatically(true);
    }

    typedef boost::property_map<FaceGraph,CGAL::face_patch_id_t<int> >::type PatchIDMap;
    FaceGraph* fg =item->face_graph();
    boost::property_map<FaceGraph, boost::vertex_index_t>::type
      vimap = get(boost::vertex_index, *fg);

    PatchIDMap pidmap = get(CGAL::face_patch_id_t<int>(), *fg);
    std::map<boost::graph_traits<FaceGraph>::vertex_descriptor,
        int> vpm;
    if(type == DUAL)
      CGAL::METIS::partition_dual_graph(*fg,
                                        nparts,
                                        CGAL::parameters::vertex_partition_id_map(boost::make_assoc_property_map(vpm)).face_partition_id_map(pidmap).vertex_index_map(vimap)
                                        );
    else if(type == NODAL)
      CGAL::METIS::partition_graph(*fg,
                                   nparts,
                                   CGAL::parameters::vertex_partition_id_map(boost::make_assoc_property_map(vpm)).face_partition_id_map(pidmap).vertex_index_map(vimap)
                                   );
    item->setProperty("NbPatchIds", nparts);
    item->invalidateOpenGLBuffers();
    QApplication::restoreOverrideCursor();
    item->redraw();
  }
  public Q_SLOTS:
  void create_nodal_partition()
  {
    create_partition(NODAL);
  }

  void create_dual_partition()
  {
    create_partition(DUAL);
  }


}; // end class Polyhedron_demo_affine_transform_plugin

#include "Partition_graph_plugin.moc"
