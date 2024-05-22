#include <unordered_set>

#include <QApplication>
#include <QAction>
#include <QList>
#include <QMainWindow>
#include <QMessageBox>
#include <QInputDialog>
#include <QtDebug>

#include "Scene_polygon_soup_item.h"
#include "Scene_polylines_item.h"
#include "Scene_surface_mesh_item.h"
#include "Scene_points_with_normal_item.h"

#include <CGAL/Three/CGAL_Lab_plugin_interface.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>

#include <CGAL/array.h>
#include <CGAL/Three/Three.h>
#include "Messages_interface.h"
#include "ui_Repair_soup.h"
using namespace CGAL::Three;
class CGAL_Lab_orient_soup_plugin :
  public QObject,
  public CGAL_Lab_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::CGAL_Lab_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.CGALLab.PluginInterface/1.0" FILE "orient_soup_plugin.json")

public:
  void init(QMainWindow* mainWindow,
            CGAL::Three::Scene_interface* scene_interface,
            Messages_interface* m);
  bool applicable(QAction* action) const {
    for(CGAL::Three::Scene_interface::Item_id index : scene->selectionIndices()) {
      if(qobject_cast<Scene_polygon_soup_item*>(scene->item(index)))
        return true;
      else
        if (action==actionShuffle &&
            qobject_cast<Scene_surface_mesh_item*>(scene->item(index))
            )
          return true;
    }
    return false;
  }

  QList<QAction*> actions() const;

public Q_SLOTS:
  void orientSM();
  void shuffle();
  void displayNonManifoldEdges();
  void createPointsAndPolyline(std::vector<std::size_t> &nm_points, bool warn);
  void cleanSoup();

private:
  template<class Item>
  void apply_shuffle(Item* item,
                     const CGAL::Three::Scene_interface::Item_id& index);


  CGAL::Three::Scene_interface* scene;
  Messages_interface* messages;
  QMainWindow* mw;
  QAction* actionOrientSM;
  QAction* actionShuffle;
  QAction* actionDisplayNonManifoldEdges;
  QAction* actionClean;

}; // end CGAL_Lab_orient_soup_plugin

void CGAL_Lab_orient_soup_plugin::init(QMainWindow* mainWindow,
                                              CGAL::Three::Scene_interface* scene_interface,
                                              Messages_interface* m)
{
  scene = scene_interface;
  mw = mainWindow;
  messages = m;
  actionOrientSM = new QAction(tr("&Orient Polygon Soup"), mainWindow);
  actionOrientSM->setObjectName("actionOrientSM");
  actionOrientSM->setProperty("subMenuName", "Polygon Mesh Processing");
  connect(actionOrientSM, SIGNAL(triggered()),
          this, SLOT(orientSM()));

  actionShuffle = new QAction(tr("&Shuffle Polygon Soup"), mainWindow);
  actionShuffle->setProperty("subMenuName", "Polygon Mesh Processing");
  connect(actionShuffle, SIGNAL(triggered()),
          this, SLOT(shuffle()));

  actionDisplayNonManifoldEdges = new QAction(tr("Display Non Manifold Edges"),
                                              mainWindow);
  actionDisplayNonManifoldEdges->setProperty("subMenuName", "View");
  connect(actionDisplayNonManifoldEdges, SIGNAL(triggered()),
          this, SLOT(displayNonManifoldEdges()));

  actionClean = new QAction(tr("Clean Polygon Soup"), mainWindow);
  actionClean->setProperty("subMenuName", "Polygon Mesh Processing");
  connect(actionClean, &QAction::triggered,
          this, &CGAL_Lab_orient_soup_plugin::cleanSoup);
}

QList<QAction*> CGAL_Lab_orient_soup_plugin::actions() const {
  return QList<QAction*>()
      << actionOrientSM
      << actionShuffle
      << actionDisplayNonManifoldEdges
      << actionClean;
}

void set_vcolors(SMesh* smesh, std::vector<CGAL::IO::Color> colors)
{
  typedef SMesh SMesh;
  typedef boost::graph_traits<SMesh>::vertex_descriptor vertex_descriptor;
  SMesh::Property_map<vertex_descriptor, CGAL::IO::Color> vcolors =
    smesh->property_map<vertex_descriptor, CGAL::IO::Color >("v:color").value();
  bool created;
  boost::tie(vcolors, created) = smesh->add_property_map<SMesh::Vertex_index,CGAL::IO::Color>("v:color",CGAL::IO::Color(0,0,0));
  assert(colors.size()==smesh->number_of_vertices());
  int color_id = 0;
  for(vertex_descriptor vd : vertices(*smesh))
      vcolors[vd] = colors[color_id++];
}

void set_fcolors(SMesh* smesh, std::vector<CGAL::IO::Color> colors)
{
  typedef SMesh SMesh;
  typedef boost::graph_traits<SMesh>::face_descriptor face_descriptor;
  SMesh::Property_map<face_descriptor, CGAL::IO::Color> fcolors =
    smesh->property_map<face_descriptor, CGAL::IO::Color >("f:color").value();
  bool created;
   boost::tie(fcolors, created) = smesh->add_property_map<SMesh::Face_index,CGAL::IO::Color>("f:color",CGAL::IO::Color(0,0,0));
  assert(colors.size()==smesh->number_of_faces());
  int color_id = 0;
  for(face_descriptor fd : faces(*smesh))
      fcolors[fd] = colors[color_id++];
}

void CGAL_Lab_orient_soup_plugin::orientSM()
{
  for(CGAL::Three::Scene_interface::Item_id index : scene->selectionIndices())
  {
    Scene_polygon_soup_item* item =
      qobject_cast<Scene_polygon_soup_item*>(scene->item(index));

    if(item)
    {
      int create_items = QMessageBox::question(mw, "Orient Mesh", "Do you wish to extract the potential non manifold simplicies ?");
      std::vector<std::size_t> nm_points;
      if(!item->orient(nm_points)) {
         QMessageBox::information(mw, tr("Not orientable without self-intersections"),
                                      tr("The polygon soup \"%1\" is not directly orientable."
                                         " Some vertices have been duplicated and some self-intersections"
                                         " have been created.")
                                      .arg(item->name()));
      }
      if(create_items == QMessageBox::Yes)
      {
        createPointsAndPolyline(nm_points, true);
      }
      QApplication::setOverrideCursor(Qt::WaitCursor);
        SMesh* smesh = new SMesh();
        if(item->exportAsSurfaceMesh(smesh)) {
          if(!item->getVColors().empty())
            set_vcolors(smesh,item->getVColors());
          if(!item->getFColors().empty())
            set_fcolors(smesh,item->getFColors());
          Scene_surface_mesh_item* sm_item = new Scene_surface_mesh_item(smesh);
          sm_item->setName(item->name());
          sm_item->setRenderingMode(item->renderingMode());
          sm_item->setVisible(item->visible());
          sm_item->setProperty("source filename", item->property("source filename"));
          sm_item->setProperty("loader_name", item->property("loader_name"));
          scene->replaceItem(index, sm_item);
          item->deleteLater();
        } else {
          item->invalidateOpenGLBuffers();
          scene->itemChanged(item);
        }
      QApplication::restoreOverrideCursor();
    }
    else{
      CGAL::Three::Three::warning(tr("This function is only applicable on polygon soups."));
    }
  }
}

void CGAL_Lab_orient_soup_plugin::shuffle()
{
  for(CGAL::Three::Scene_interface::Item_id index : scene->selectionIndices())
  {
    Scene_polygon_soup_item* soup_item =
        qobject_cast<Scene_polygon_soup_item*>(scene->item(index));

    if(soup_item) {
      QApplication::setOverrideCursor(Qt::WaitCursor);
      soup_item->shuffle_orientations();
      QApplication::restoreOverrideCursor();
      continue;
    }
    Scene_surface_mesh_item* sm_item =
        qobject_cast<Scene_surface_mesh_item*>(scene->item(index));
    if(sm_item) {
      apply_shuffle(sm_item, index);
    }
  }
}

template<class Item>
void CGAL_Lab_orient_soup_plugin::apply_shuffle( Item* root_item,
                                                        const CGAL::Three::Scene_interface::Item_id& index)
{
  QApplication::setOverrideCursor(Qt::WaitCursor);
  Scene_polygon_soup_item* item = new Scene_polygon_soup_item();
  item->setName(root_item->name());
  item->setRenderingMode(root_item->renderingMode());
  item->setVisible(root_item->visible());
  item->setProperty("source filename", root_item->property("source filename"));
  item->load(root_item);
  item->shuffle_orientations();
  item->setColor(root_item->color());
  scene->replaceItem(index, item, true);
  root_item->deleteLater();
  QApplication::restoreOverrideCursor();
}

void CGAL_Lab_orient_soup_plugin::displayNonManifoldEdges()
{
  const CGAL::Three::Scene_interface::Item_id index = scene->mainSelectionIndex();

  Scene_polygon_soup_item* item =
    qobject_cast<Scene_polygon_soup_item*>(scene->item(index));

  if(item)
  {
    QApplication::setOverrideCursor(Qt::WaitCursor);
    item->setDisplayNonManifoldEdges(!item->displayNonManifoldEdges());
    scene->itemChanged(item);
    QApplication::restoreOverrideCursor();
  }
}
//todo: nm-points should probably be a pair, and the removal check be on both members
void CGAL_Lab_orient_soup_plugin::createPointsAndPolyline(std::vector<std::size_t>& nm_points, bool warn)
{

  const CGAL::Three::Scene_interface::Item_id index = scene->mainSelectionIndex();

  Scene_polygon_soup_item* item =
    qobject_cast<Scene_polygon_soup_item*>(scene->item(index));

  if(item)
  {
    QApplication::setOverrideCursor(Qt::WaitCursor);
    Scene_points_with_normal_item* points = nullptr;

    bool items_created = false;
    if(nm_points.empty() && warn)
    {
      delete points;
      CGAL::Three::Three::information(tr("There is no non-manifold vertex in this soup."));
    }
    else
    {
      points = new Scene_points_with_normal_item();
      items_created = true;
      for(std::size_t id : nm_points)
      {
        points->point_set()->insert(item->points()[id]);
      }
      points->setName(QString("Non Manifold Vertices of %1").arg(item->name()));
          points->setColor(QColor(Qt::red));

          scene->addItem(points);
    }
    Polygon_soup::Edges nm_edges = item->non_manifold_edges();
    if(!nm_edges.empty())
    {
      Scene_polylines_item* poly =
          new Scene_polylines_item();
      for(Polygon_soup::Edge edge : nm_edges)
      {
        Point_3 a(item->points()[edge[0]]), b(item->points()[edge[1]]);
        Scene_polylines_item::Polyline new_edge;
        new_edge.push_back(a);
        new_edge.push_back(b);
        poly->polylines.push_back(new_edge);
      }
      poly->setName(QString("Non Manifold Edges of %1").arg(item->name()));
      poly->setColor(QColor(Qt::red));

      scene->addItem(poly);
      items_created = true;
    }
    else
    {
      CGAL::Three::Three::information(tr("There is no non-manifold edge in this soup."));
    }
    QApplication::restoreOverrideCursor();
    if(!items_created)
      QMessageBox::information(mw, "Nothing Non-manifold", "No non-manifold edge nor vertex was found.");
  }
}



class RepairDialog :
    public QDialog,
    public Ui::Dialog
{
  Q_OBJECT
public:
  RepairDialog(QWidget* =nullptr)
  {
    setupUi(this);
  }
};

void CGAL_Lab_orient_soup_plugin::cleanSoup()
{
  const CGAL::Three::Scene_interface::Item_id index = scene->mainSelectionIndex();

  Scene_polygon_soup_item* item =
      qobject_cast<Scene_polygon_soup_item*>(scene->item(index));

  if(!item)
    return;
  RepairDialog dlg;
  if(!dlg.exec())
    return;
  bool b1 = dlg.eadCheckbox->isChecked(), b2 = dlg.rsoCheckBox->isChecked();
  item->repair(b1, b2);

  item->invalidateOpenGLBuffers();
  item->itemChanged();
}
#include "Orient_soup_plugin.moc"
