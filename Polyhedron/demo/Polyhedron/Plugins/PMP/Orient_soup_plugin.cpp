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

#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>

#include <CGAL/array.h>
#include <CGAL/Three/Three.h>
#include "Messages_interface.h"
#include "ui_Repair_soup.h"
using namespace CGAL::Three;
class Polyhedron_demo_orient_soup_plugin :
  public QObject,
  public Polyhedron_demo_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0" FILE "orient_soup_plugin.json")

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
  void createPointsAndPolyline();
  void cleanSoup();

private:
  template<class Item>
  void apply_shuffle(Item* item,
                     const CGAL::Three::Scene_interface::Item_id& index);
  void getNMPoints(std::set<std::size_t> &vertices_to_duplicate,
      Scene_polygon_soup_item* item);

  CGAL::Three::Scene_interface* scene;
  Messages_interface* messages;
  QMainWindow* mw;
  QAction* actionOrientSM;
  QAction* actionShuffle;
  QAction* actionNMToPolyline;
  QAction* actionDisplayNonManifoldEdges;
  QAction* actionClean;

}; // end Polyhedron_demo_orient_soup_plugin

void Polyhedron_demo_orient_soup_plugin::init(QMainWindow* mainWindow,
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
  actionNMToPolyline = new QAction(tr("Extract Non Manifold Simplices"), mainWindow);
  actionNMToPolyline->setProperty("subMenuName", "Polygon Mesh Processing");
  connect(actionNMToPolyline, &QAction::triggered,
          this, &Polyhedron_demo_orient_soup_plugin::createPointsAndPolyline);
  actionClean = new QAction(tr("Clean Polygon Soup"), mainWindow);
  actionClean->setProperty("subMenuName", "Polygon Mesh Processing");
  connect(actionClean, &QAction::triggered,
          this, &Polyhedron_demo_orient_soup_plugin::cleanSoup);
}

QList<QAction*> Polyhedron_demo_orient_soup_plugin::actions() const {
  return QList<QAction*>()
      << actionOrientSM
      << actionShuffle
      << actionNMToPolyline
      << actionDisplayNonManifoldEdges
      << actionClean;
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
  for(vertex_descriptor vd : vertices(*smesh))
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
  for(face_descriptor fd : faces(*smesh))
      fcolors[fd] = colors[color_id++];
}

void Polyhedron_demo_orient_soup_plugin::orientSM()
{
  Q_FOREACH(CGAL::Three::Scene_interface::Item_id index, scene->selectionIndices())
  {
    Scene_polygon_soup_item* item =
      qobject_cast<Scene_polygon_soup_item*>(scene->item(index));

    if(item)
    {
      int create_items = QMessageBox::question(mw, "Orient Mesh", "Do you wish to extract the potential non manifold simplicies ?");
      if(create_items == QMessageBox::Yes)
      {
        createPointsAndPolyline();
      }
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
      CGAL::Three::Three::warning(tr("This function is only applicable on polygon soups."));
    }
  }
}

void Polyhedron_demo_orient_soup_plugin::shuffle()
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
  scene->replaceItem(index, item, true);
  root_item->deleteLater();
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
void Polyhedron_demo_orient_soup_plugin::createPointsAndPolyline()
{

  const CGAL::Three::Scene_interface::Item_id index = scene->mainSelectionIndex();

  Scene_polygon_soup_item* item =
    qobject_cast<Scene_polygon_soup_item*>(scene->item(index));

  if(item)
  {
    QApplication::setOverrideCursor(Qt::WaitCursor);
    Scene_points_with_normal_item* points = new Scene_points_with_normal_item();
    std::set<std::size_t> nm_vertices;
    getNMPoints(nm_vertices, item);
    bool items_created = false;
    if(nm_vertices.empty())
    {
      delete points;
      CGAL::Three::Three::information(tr("There is no non-manifold vertex in this soup."));
    }
    else
    {
        items_created = true;
      for(std::size_t id : nm_vertices)
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

void Polyhedron_demo_orient_soup_plugin::getNMPoints(
    std::set<std::size_t > &vertices_to_duplicate,
    Scene_polygon_soup_item* item)
{
  typedef std::pair<std::size_t, std::size_t>                              V_ID_pair;
  typedef CGAL::Polygon_mesh_processing::internal::Polygon_soup_orienter<Polygon_soup::Points,
      Polygon_soup::Polygons> PSO;
  typedef PSO::Edge_map Edge_map;
  typedef std::set<V_ID_pair>                                              Marked_edges;

  Edge_map edges;
  edges.resize(item->points().size());
  Marked_edges m_edges;
  PSO::fill_edge_map(edges, m_edges, item->polygons());

  // for each vertex, indicates the list of polygon containing it
  std::vector< std::vector<std::size_t> > incident_polygons_per_vertex(item->points().size());
  std::size_t nb_polygons=item->polygons().size();
  for(std::size_t ip=0; ip<nb_polygons; ++ip)
  {
    for(std::size_t iv : item->polygons()[ip])
      incident_polygons_per_vertex[iv].push_back(ip);
  }

  std::size_t nbv = item->points().size();

  for (std::size_t v_id = 0; v_id < nbv; ++v_id)
  {
    const std::vector< std::size_t >& incident_polygons = incident_polygons_per_vertex[v_id];

    if ( incident_polygons.empty() ) continue; //isolated vertex
    std::set<std::size_t> visited_polygons;

    bool first_pass = true;
    for(std::size_t p_id : incident_polygons)
    {
      if ( !visited_polygons.insert(p_id).second ) continue; // already visited

      if (!first_pass)
      {
        vertices_to_duplicate.insert(v_id);
      }


      std::size_t nbv = item->polygons()[p_id].size(), pvid=0;
      for (; pvid!=nbv; ++pvid)
        if (v_id==item->polygons()[p_id][pvid]) break;
      CGAL_assertion( pvid!=nbv );
      std::size_t p = item->polygons()[p_id][ (pvid+nbv-1)%nbv ];
      std::size_t n = item->polygons()[p_id][ (pvid+1)%nbv ];
      const std::array<std::size_t,3>& neighbors = CGAL::make_array(p,v_id,n);

      std::size_t next = neighbors[2];

      do{
        std::size_t other_p_id;
        std::tie(next, other_p_id) = PSO::next_cw_vertex_around_source(v_id, next, item->polygons(), edges, m_edges);
        if (next==v_id) break;
        visited_polygons.insert(other_p_id);
      }
      while(next!=neighbors[0]);

      if (next==v_id){
        /// turn the otherway round
        next = neighbors[0];
        do{
          std::size_t other_p_id;
          std::tie(next, other_p_id) = PSO::next_ccw_vertex_around_target(next, v_id, item->polygons(), edges, m_edges);
          if (next==v_id) break;
          visited_polygons.insert(other_p_id);
        }
        while(true);
      }
      first_pass=false;
    }
  }

  //remove vertices already in NM edges
  //check edges of p_id.
  for(Scene_polygon_soup_item::Edge edge : item->non_manifold_edges())
  {
    vertices_to_duplicate.erase(edge[0]);
    vertices_to_duplicate.erase(edge[1]);
  }
}


class RepairDialog :
    public QDialog,
    public Ui::Dialog
{
  Q_OBJECT
public:
  RepairDialog(QWidget* =0)
  {
    setupUi(this);
  }
};

void Polyhedron_demo_orient_soup_plugin::cleanSoup()
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

