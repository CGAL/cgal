
// Needed for lloyd_optimize_mesh_2 which does it too late
// (and we don't want to spend the time on finding out who
// includes the header file that sets it too a value too low
#define  BOOST_PARAMETER_MAX_ARITY 8

#include <stdexcept>

#define CGAL_CT2_WANTS_TO_HAVE_EXTRA_ACTION_FOR_INTERSECTING_CONSTRAINTS
#define CGAL_CDT2_EXTRA_ACTION_FOR_INTERSECTING_CONSTRAINTS \
  throw std::runtime_error("This plugin does not deal with intersecting edges");

#include <QtCore/qglobal.h>

#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>

#include "Scene_polyhedron_item.h"
#include "Scene_surface_mesh_item.h"
#include "Scene_polylines_item.h"
#include "Scene_points_with_normal_item.h"
#include "Polyhedron_type.h"

#include <CGAL/iterator.h>

#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Projection_traits_yz_3.h>
#include <CGAL/Projection_traits_xz_3.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesh_vertex_base_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/lloyd_optimize_mesh_2.h>

#include <CGAL/boost/graph/Euler_operations.h>

#include <QTime>
#include <QAction>
#include <QMainWindow>
#include <QApplication>
#include <QString>
#include <QDialog>
#include <QMessageBox>
#include <QtPlugin>

#include <vector>
#include <algorithm>
#include <queue>
#include <sstream>

#include "ui_mesh_2_dialog.h"

struct FaceInfo2
{
  FaceInfo2(){}
  int nesting_level;
  bool in_domain(){
    return nesting_level%2 == 1;
  }
};

template <class CDT>
void
mark_domains(CDT& ct,
             typename CDT::Face_handle start,
             int index,
             std::list<typename CDT::Edge>& border )
{
  if(start->info().nesting_level != -1){
    return;
  }
  std::list<typename CDT::Face_handle> queue;
  queue.push_back(start);
  while(! queue.empty()){
    typename CDT::Face_handle fh = queue.front();
    queue.pop_front();
    if(fh->info().nesting_level == -1){
      fh->info().nesting_level = index;
      for(int i = 0; i < 3; i++){
        typename CDT::Edge e(fh,i);
        typename CDT::Face_handle n = fh->neighbor(i);
        if(n->info().nesting_level == -1){
          if(ct.is_constrained(e)) border.push_back(e);
          else queue.push_back(n);
        }
      }
    }
  }
}

template <class CDT>
void
mark_nested_domains(CDT& cdt)
{
  for(typename CDT::All_faces_iterator it = cdt.all_faces_begin(); it != cdt.all_faces_end(); ++it){
    it->info().nesting_level = -1;
  }
  std::list<typename CDT::Edge> border;
  mark_domains(cdt, cdt.infinite_face(), 0, border);
  while(! border.empty()){
    typename CDT::Edge e = border.front();
    border.pop_front();
    typename CDT::Face_handle n = e.first->neighbor(e.second);
    if(n->info().nesting_level == -1){
      mark_domains(cdt, n, e.first->info().nesting_level+1, border);
    }
  }
}

template <class CDT, class TriangleMesh>
void cdt2_to_face_graph(const CDT& cdt, TriangleMesh& tm, int constant_coordinate_index, double constant_coordinate)
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor vertex_descriptor;

  typedef std::map<typename CDT::Vertex_handle, vertex_descriptor> Map;
  Map descriptors;
  for (typename CDT::Finite_faces_iterator fit=cdt.finite_faces_begin(),
                                           fit_end=cdt.finite_faces_end();
                                           fit!=fit_end; ++fit)
  {
    if (!fit->is_in_domain()) continue;
    CGAL::cpp11::array<vertex_descriptor,3> vds;
    for(int i=0; i<3; ++i)
    {
      typename Map::iterator it;
      bool insert_ok;
      boost::tie(it,insert_ok) =
        descriptors.insert(std::make_pair(fit->vertex(i),vertex_descriptor()));
      if (insert_ok){
        const Kernel::Point_3& pt=fit->vertex(i)->point();
        double coords[3] = {pt[0], pt[1], pt[2]};
        coords[constant_coordinate_index]=constant_coordinate;
        it->second = add_vertex(Kernel::Point_3(coords[0],coords[1],coords[2]), tm);
      }
      vds[i]=it->second;
    }

    CGAL::Euler::add_face(vds, tm);
  }
}

class Polyhedron_demo_mesh_2_plugin :
  public QObject,
  public CGAL::Three::Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

public:
  void init(QMainWindow* mainWindow,
            CGAL::Three::Scene_interface* scene_interface,
            Messages_interface*)
  {
    this->scene = scene_interface;
    this->mw = mainWindow;

    actionMesh_2_ = new QAction("Mesh_2", mw);
    // actionMesh_2_->setProperty("subMenuName", "Polygon Mesh Processing");
    if (actionMesh_2_) {
      connect(actionMesh_2_, SIGNAL(triggered()),
        this, SLOT(run()));
    }
  }

  QList<QAction*> actions() const {
    return QList<QAction*>() << actionMesh_2_;
  }

  bool applicable(QAction*) const
  {
    Q_FOREACH(int index, scene->selectionIndices())
    {
      //if one polyhedron is found in the selection, it's fine
      if (qobject_cast<Scene_polylines_item*>(scene->item(index)))
        return true;
    }
    return false;
  }

private:
  Ui::mesh_2_dialog
  create_dialog(QDialog* dialog, double diag_length, bool no_seeds)
  {
    Ui::mesh_2_dialog ui;
    ui.setupUi(dialog);
    connect(ui.buttonBox, SIGNAL(accepted()), dialog, SLOT(accept()));
    connect(ui.buttonBox, SIGNAL(rejected()), dialog, SLOT(reject()));

    //connect checkbox to spinbox
    connect(ui.runLloyd_checkbox, SIGNAL(toggled(bool)),
            ui.nbIterations_spinbox, SLOT(setEnabled(bool)));
    connect(ui.runMesh2_checkbox, SIGNAL(toggled(bool)),
            ui.edgeLength_dspinbox, SLOT(setEnabled(bool)));

    //Set default parameter edge length
    ui.edgeLength_dspinbox->setDecimals(3);
    ui.edgeLength_dspinbox->setSingleStep(0.001);
    ui.edgeLength_dspinbox->setRange(1e-6 * diag_length, //min
                                     2.   * diag_length);//max
    ui.edgeLength_dspinbox->setValue(0.05 * diag_length);
    ui.edgeLength_dspinbox->setToolTip(
      "Diagonal length of the Bbox of the selection to mesh is "+
      QString::number(diag_length)+" - default is 5% of it");

    //Set default for nb iterations
    ui.nbIterations_spinbox->setSingleStep(1);
    ui.nbIterations_spinbox->setRange(1/*min*/, 1000/*max*/);
    ui.nbIterations_spinbox->setValue(1);

    // Run Mesh_2 by default, not Lloyd
    ui.runLloyd_checkbox->setChecked(false);
    ui.runMesh2_checkbox->setChecked(true);

    // Domain definition and disabling all options if no Mesh_2 run
    if (no_seeds){
      ui.radioSeedsIn->setDisabled(true);
      ui.radioSeedsOut->setDisabled(true);
      ui.radioNesting->setChecked(true);
    }
    else
      ui.radioSeedsOut->setChecked(true);
    return ui;
  }

  template <class ProjectionTraits>
  void mesh(const std::vector<Scene_polylines_item*>& polylines_items,
            const std::vector<Scene_points_with_normal_item*>& points_items,
            double diagonal_length,
            int constant_coordinate_index)
  {
    // extract seeds
    std::vector<Kernel::Point_3> seeds;
    Q_FOREACH(Scene_points_with_normal_item* points_item, points_items)
      Q_FOREACH(Point_set_3<Kernel>::Index it, *points_item->point_set())
        seeds.push_back(points_item->point_set()->point(it));

    // Create dialog box
    QDialog dialog(mw);
    Ui::mesh_2_dialog ui =
      create_dialog(&dialog, diagonal_length, seeds.empty());
    dialog.setWindowFlags(Qt::Dialog|Qt::CustomizeWindowHint|Qt::WindowCloseButtonHint);

    // Get values
    int i = dialog.exec();
    if (i == QDialog::Rejected)
    {
      std::cout << "2D Meshing aborted" << std::endl;
      return;
    }
    bool runMesh2 = ui.runMesh2_checkbox->isChecked();
    double target_length = ui.edgeLength_dspinbox->value();
    unsigned int nb_iter = ui.nbIterations_spinbox->value();
    bool runLloyd = ui.runLloyd_checkbox->isChecked();


    // wait cursor
    QApplication::setOverrideCursor(Qt::WaitCursor);
    typedef ProjectionTraits                                             Gt;
    typedef CGAL::Delaunay_mesh_vertex_base_2<Gt>                        Vb;
    typedef CGAL::Delaunay_mesh_face_base_2<Gt>                          Fm;
    typedef CGAL::Triangulation_face_base_with_info_2<FaceInfo2,Gt,Fm>   Fb;
    typedef CGAL::Triangulation_data_structure_2<Vb, Fb>                TDS;
    typedef CGAL::No_intersection_tag                                   Tag;
    typedef CGAL::Constrained_Delaunay_triangulation_2<Gt, TDS, Tag>    CDT;
    typedef CGAL::Delaunay_mesh_size_criteria_2<CDT>               Criteria;
    typedef CGAL::Delaunay_mesher_2<CDT, Criteria>                   Mesher;

    QTime time; // global timer
    time.start();

    std::cout << " Building Constrained_Delaunay_triangulation_2..."
              << std::flush;
    CDT cdt;

    QTime ltime; //local timer
    ltime.start();
    double constant_coordinate =
      polylines_items.back()->polylines.back().back()[constant_coordinate_index];
    try{
      Q_FOREACH(Scene_polylines_item* polylines_item, polylines_items)
        Q_FOREACH(const std::vector<Kernel::Point_3>& points,
                    polylines_item->polylines)
          cdt.insert_constraint(points.begin(),points.end());
    }catch(std::runtime_error)
    {
      QApplication::restoreOverrideCursor();
      throw;
    }
    std::cout << " done (" << ltime.elapsed() << " ms)" << std::endl;

    if (cdt.dimension()!=2){
      QApplication::restoreOverrideCursor();
      std::cout << "Triangulation is not of dimension 2" << std::endl;
      return;
    }

    // start by marking the domain to mesh
    Criteria criteria(0.125, target_length);
    Mesher mesher(cdt, criteria);
    bool use_seeds=ui.radioSeedsOut->isChecked() ||
                   ui.radioSeedsIn->isChecked();
    bool use_nesting = ui.radioNesting->isChecked();

    if (!use_seeds){
      if (use_nesting){
        mark_nested_domains(cdt);
        for(typename CDT::All_faces_iterator fit=cdt.all_faces_begin(),
                                             fit_end=cdt.all_faces_end();
                                             fit!=fit_end;++fit)
        {
          fit->set_in_domain(fit->info().in_domain());
        }
      }
      mesher.init(use_nesting);
    }
    else
      mesher.set_seeds(seeds.begin(), seeds.end(), ui.radioSeedsIn->isChecked(), true);

    if (runMesh2){
      time.restart();
      std::cout << " Running refine_Delaunay_mesh_2 ..." << std::flush;
      mesher.refine_mesh();
      std::cout << " done (" << ltime.elapsed() << " ms)" << std::endl;
    }

    if (runLloyd){
      ltime.restart();
      std::cout << " Running lloyd_optimize_mesh_2..." << std::flush;
      CGAL::lloyd_optimize_mesh_2(cdt,
        CGAL::parameters::max_iteration_number = nb_iter);
      std::cout << " done (" << ltime.elapsed() << " ms)" << std::endl;
    }

    // export result as a polyhedron item
    QString iname =
      polylines_items.size()==1?
      polylines_items.front()->name()+QString("_meshed_"):
      QString("2dmesh_");
    iname+=QString::number(target_length);
    if (runLloyd) iname+=QString("_Lloyd_")+QString::number(nb_iter);
    
    if(mw->property("is_polyhedron_mode").toBool()){
      Scene_polyhedron_item* poly_item = new Scene_polyhedron_item();
      poly_item->setName(iname);
      cdt2_to_face_graph(cdt,
                         *poly_item->polyhedron(),
                         constant_coordinate_index,
                         constant_coordinate);
      scene->addItem(poly_item);
      poly_item->invalidateOpenGLBuffers();
    }else{
      Scene_surface_mesh_item* poly_item = new Scene_surface_mesh_item();
      poly_item->setName(iname);
      cdt2_to_face_graph(cdt,
                         *poly_item->polyhedron(),
                         constant_coordinate_index,
                         constant_coordinate);
      scene->addItem(poly_item);
      poly_item->invalidateOpenGLBuffers();
    }
    std::cout << "ok (" << time.elapsed() << " ms)" << std::endl;
    // default cursor
    QApplication::restoreOverrideCursor();
  }

  int detect_constant_coordinate(
    const std::vector<Scene_polylines_item*>& polylines_items,
    const std::vector<Scene_points_with_normal_item*>& points_items)
  {
    int res=-1;
    Kernel::Point_3 ref = polylines_items.front()->polylines.front().front();
    Q_FOREACH(Scene_polylines_item* polylines_item, polylines_items)
      Q_FOREACH(const std::vector<Kernel::Point_3>& points,
                  polylines_item->polylines)
        Q_FOREACH(const Kernel::Point_3& pt, points)
        {
          int nbe=0, candidate=-1;
          for (int i=0; i<3; ++i)
            if (ref[i]==pt[i]){
              ++nbe;
              candidate=i;
            }
          if (nbe==0) return -1;
          if (nbe==1){
            if (res==-1)
              res=candidate;
            else
              if(res!=candidate) return -1;
          }
        }
    if (res==-1) return res;
    Q_FOREACH(Scene_points_with_normal_item* points_item, points_items)
      Q_FOREACH(Point_set_3<Kernel>::Index pt, *points_item->point_set())
      if (points_item->point_set()->point(pt)[res]!=ref[res])
          return -1;
    return res;
  }

public Q_SLOTS:

  void run()
  {
    QApplication::setOverrideCursor(Qt::WaitCursor);
    //collect input polylines
    std::vector<Scene_polylines_item*> polylines_items;
    std::vector<Scene_points_with_normal_item*> points_items;
    double inf = std::numeric_limits<double>::infinity();
    CGAL::Three::Scene_interface::Bbox bbox(inf,inf,inf,-inf,-inf,-inf);
    Q_FOREACH(int index, scene->selectionIndices())
    {
      Scene_polylines_item* polylines_item =
        qobject_cast<Scene_polylines_item*>(scene->item(index));
      if (polylines_item){
        polylines_items.push_back(polylines_item);
        bbox=bbox+polylines_item->bbox();
      }
      else{
        Scene_points_with_normal_item* points_item =
          qobject_cast<Scene_points_with_normal_item*>(scene->item(index));
        if (points_item)
          points_items.push_back(points_item);
      }
    }

    double diag = CGAL::sqrt(
          (bbox.xmax()-bbox.xmin())*(bbox.xmax()-bbox.xmin())
          +(bbox.ymax()-bbox.ymin())*(bbox.ymax()-bbox.ymin())
          +(bbox.zmax()-bbox.zmax()) *(bbox.zmax()-bbox.zmax())
          );
    QApplication::restoreOverrideCursor();
    switch( detect_constant_coordinate(polylines_items, points_items) )
    {
      using namespace CGAL;
      typedef Kernel K;
      case 0:
        mesh<Projection_traits_yz_3<K> >(polylines_items, points_items, diag, 0);
      break;
      case 1:
        mesh<Projection_traits_xz_3<K> >(polylines_items, points_items, diag, 1);
      break;
      case 2:
        mesh<Projection_traits_xy_3<K> >(polylines_items, points_items, diag, 2);
      break;
      default:
        QMessageBox::critical(mw,
                              "Invalid Input Data",
                              "Polylines and seed points must all be in "
                              "the same xy, yz or xz plane");
    }
  }

private:
  QAction* actionMesh_2_;

}; // end Mesh_2_plugin

#include "Mesh_2_plugin.moc"
