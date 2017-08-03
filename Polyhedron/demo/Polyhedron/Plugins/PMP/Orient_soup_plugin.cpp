#include <QApplication>
#include <QAction>
#include <QList>
#include <QMainWindow>
#include <QMessageBox>
#include <QtDebug>

#include "Scene_polygon_soup_item.h"
#include "Scene_polyhedron_item.h"
#include "Scene_surface_mesh_item.h"

#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>
#include "Messages_interface.h"
using namespace CGAL::Three;
class Polyhedron_demo_orient_soup_plugin : 
  public QObject,
  public Polyhedron_demo_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

public:
  void init(QMainWindow* mainWindow,
            CGAL::Three::Scene_interface* scene_interface,
            Messages_interface* m);
  bool applicable(QAction* action) const {
    Q_FOREACH(CGAL::Three::Scene_interface::Item_id index, scene->selectionIndices()) {
      if(qobject_cast<Scene_polygon_soup_item*>(scene->item(index)))
        return true;
      else
        if (action==actionShuffle &&
            (qobject_cast<Scene_polyhedron_item*>(scene->item(index))
            ||qobject_cast<Scene_surface_mesh_item*>(scene->item(index)))
            )
          return true;
    }
    return false;
  }

  QList<QAction*> actions() const;

public Q_SLOTS:
  void orientPoly();
  void orientSM();
  void shuffle();
  void displayNonManifoldEdges();

private:
  template<class Item>
  void apply_shuffle(Item* item,
                     const CGAL::Three::Scene_interface::Item_id& index);
  CGAL::Three::Scene_interface* scene;
  Messages_interface* messages;
  QMainWindow* mw;
  QAction* actionOrientPoly;
  QAction* actionOrientSM;
  QAction* actionShuffle;
  QAction* actionDisplayNonManifoldEdges;

}; // end Polyhedron_demo_orient_soup_plugin

void Polyhedron_demo_orient_soup_plugin::init(QMainWindow* mainWindow,
                                              CGAL::Three::Scene_interface* scene_interface,
                                              Messages_interface* m)
{
  scene = scene_interface;
  mw = mainWindow;
  messages = m;
  actionOrientPoly = new QAction(tr("&Orient Polygon Soup (as a polyhedron)"), mainWindow);
  actionOrientPoly->setObjectName("actionOrientPoly");
  actionOrientPoly->setProperty("subMenuName", "Polygon Mesh Processing");
  connect(actionOrientPoly, SIGNAL(triggered()),
          this, SLOT(orientPoly()));
  actionOrientSM = new QAction(tr("&Orient Polygon Soup (as a surface_mesh)"), mainWindow);
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
}

QList<QAction*> Polyhedron_demo_orient_soup_plugin::actions() const {
  return QList<QAction*>()
      << actionOrientPoly
      << actionOrientSM
      << actionShuffle
      << actionDisplayNonManifoldEdges;
}

void set_vcolors(SMesh* smesh, std::vector<CGAL::Color> colors)
{
  typedef SMesh SMesh;
  typedef boost::graph_traits<SMesh>::vertex_descriptor vertex_descriptor;
  SMesh::Property_map<vertex_descriptor, CGAL::Color> vcolors =
    smesh->property_map<vertex_descriptor, CGAL::Color >("v:color").first;
  bool created;
  boost::tie(vcolors, created) = smesh->add_property_map<SMesh::Vertex_index,CGAL::Color>("v:color",CGAL::Color(0,0,0));
  assert(colors.size()==smesh->number_of_vertices());
  int color_id = 0;
  BOOST_FOREACH(vertex_descriptor vd, vertices(*smesh))
      vcolors[vd] = colors[color_id++];
}

void set_fcolors(SMesh* smesh, std::vector<CGAL::Color> colors)
{
  typedef SMesh SMesh;
  typedef boost::graph_traits<SMesh>::face_descriptor face_descriptor;
  SMesh::Property_map<face_descriptor, CGAL::Color> fcolors =
    smesh->property_map<face_descriptor, CGAL::Color >("f:color").first;
  bool created;
   boost::tie(fcolors, created) = smesh->add_property_map<SMesh::Face_index,CGAL::Color>("f:color",CGAL::Color(0,0,0));
  assert(colors.size()==smesh->number_of_faces());
  int color_id = 0;
  BOOST_FOREACH(face_descriptor fd, faces(*smesh))
      fcolors[fd] = colors[color_id++];
}

void Polyhedron_demo_orient_soup_plugin::orientPoly()
{
  Q_FOREACH(CGAL::Three::Scene_interface::Item_id index, scene->selectionIndices())
  {
    Scene_polygon_soup_item* item =
      qobject_cast<Scene_polygon_soup_item*>(scene->item(index));

    if(item)
    {
      if(!item->orient()) {
         QMessageBox::information(mw, tr("Not orientable without self-intersections"),
                                      tr("The polygon soup \"%1\" is not directly orientable."
                                         " Some vertices have been duplicated and some self-intersections"
                                         " have been created.")
                                      .arg(item->name()));
      }


        Scene_polyhedron_item* poly_item = new Scene_polyhedron_item();
        if(item->exportAsPolyhedron(poly_item->polyhedron())) {
          poly_item->setName(item->name());
          poly_item->setColor(item->color());
          poly_item->setRenderingMode(item->renderingMode());
          poly_item->setVisible(item->visible());
          poly_item->invalidateOpenGLBuffers();
          poly_item->setProperty("source filename", item->property("source filename"));
          poly_item->setProperty("loader_name", item->property("loader_name"));
          scene->replaceItem(index, poly_item);
          item->deleteLater();
        } else {
          item->invalidateOpenGLBuffers();
          scene->itemChanged(item);
        }
    }
    else{
      messages->warning(tr("This function is only applicable on polygon soups."));
    }
  }
}

void Polyhedron_demo_orient_soup_plugin::orientSM()
{
  Q_FOREACH(CGAL::Three::Scene_interface::Item_id index, scene->selectionIndices())
  {
    Scene_polygon_soup_item* item =
      qobject_cast<Scene_polygon_soup_item*>(scene->item(index));

    if(item)
    {
      //     qDebug()  << tr("I have the item %1\n").arg(item->name());
      if(!item->orient()) {
         QMessageBox::information(mw, tr("Not orientable without self-intersections"),
                                      tr("The polygon soup \"%1\" is not directly orientable."
                                         " Some vertices have been duplicated and some self-intersections"
                                         " have been created.")
                                      .arg(item->name()));
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
      messages->warning(tr("This function is only applicable on polygon soups."));
    }
  }
}

void Polyhedron_demo_orient_soup_plugin::shuffle()
{
  BOOST_FOREACH(CGAL::Three::Scene_interface::Item_id index, scene->selectionIndices())
  {
    Scene_polygon_soup_item* soup_item =
        qobject_cast<Scene_polygon_soup_item*>(scene->item(index));

    if(soup_item) {
      QApplication::setOverrideCursor(Qt::WaitCursor);
      soup_item->shuffle_orientations();
      QApplication::restoreOverrideCursor();
      continue;
    }

    Scene_polyhedron_item* poly_item =
        qobject_cast<Scene_polyhedron_item*>(scene->item(index));
    if(poly_item) {
      apply_shuffle(poly_item, index);
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
void Polyhedron_demo_orient_soup_plugin::apply_shuffle( Item* root_item,
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
  scene->replaceItem(index, item);
  delete root_item;
  QApplication::restoreOverrideCursor();
}

void Polyhedron_demo_orient_soup_plugin::displayNonManifoldEdges()
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

#include "Orient_soup_plugin.moc"
