#define CGAL_TETRAHEDRAL_REMESHING_VERBOSE
#define CGAL_DUMP_REMESHING_STEPS
#define CGAL_TETRAHEDRAL_REMESHING_DEBUG

#include <QtCore/qglobal.h>

#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>

#include "Scene_c3t3_item.h"
#include "C3t3_type.h"

#include <CGAL/Tetrahedral_remeshing/Remeshing_triangulation_3.h>
#include <CGAL/tetrahedral_remeshing.h>

#include <QTime>
#include <QAction>
#include <QMainWindow>
#include <QApplication>
#include <QString>
#include <QInputDialog>
#include <QtPlugin>
#include <QMessageBox>

#include "ui_Tetrahedral_remeshing_dialog.h"

namespace CGAL {

  namespace internal {

    template<typename TDS_src, typename TDS_tgt>
    struct Vertex_converter
    {
      typename TDS_tgt::Vertex operator()(const typename TDS_src::Vertex& v_src) const
      {
        typedef typename CGAL::Kernel_traits<
          typename TDS_src::Vertex::Point>::Kernel GT_src;
        typedef typename CGAL::Kernel_traits<
          typename TDS_tgt::Vertex::Point>::Kernel GT_tgt;
        CGAL::Cartesian_converter<GT_src, GT_tgt> conv;

        typedef typename TDS_tgt::Vertex::Point Tgt_point;

        typename TDS_tgt::Vertex v_tgt;
        v_tgt.set_point(Tgt_point(conv(v_src.point())));
        v_tgt.set_time_stamp(-1);
//        v_tgt.set_dimension(v_src.dimension());
        return v_tgt;
      }
      void operator()(const typename TDS_src::Vertex& v_src,
                      typename TDS_tgt::Vertex& v_tgt) const
      {
        typedef typename CGAL::Kernel_traits<
          typename TDS_src::Vertex::Point>::Kernel GT_src;
        typedef typename CGAL::Kernel_traits<
          typename TDS_tgt::Vertex::Point>::Kernel GT_tgt;
        CGAL::Cartesian_converter<GT_src, GT_tgt> conv;

        typedef typename TDS_tgt::Vertex::Point Tgt_point;

        v_tgt.set_point(Tgt_point(conv(v_src.point())));
//        v_tgt.set_dimension(v_src.dimension());
      }
    };

    template<typename TDS_src, typename TDS_tgt>
    struct Cell_converter
    {
      typename TDS_tgt::Cell operator()(const typename TDS_src::Cell& c_src) const
      {
        typename TDS_tgt::Cell c_tgt;
        c_tgt.set_subdomain_index(c_src.subdomain_index());
        c_tgt.set_time_stamp(-1);
        return c_tgt;
      }
      void operator()(const typename TDS_src::Cell& c_src,
                      typename TDS_tgt::Cell& c_tgt) const
      {
        c_tgt.set_subdomain_index(c_src.subdomain_index());
      }
    };

    template<typename T3, typename Remeshing_tr>
    void build_remeshing_triangulation(const T3& tr,
                                       Remeshing_tr& remeshing_tr)
    {
      typedef typename T3::Triangulation_data_structure Tds;
      typedef typename Remeshing_tr::Tds                RTds;

      remeshing_tr.clear();
      remeshing_tr.set_infinite_vertex(
        remeshing_tr.tds().copy_tds(
          tr.tds(),
          tr.infinite_vertex(),
          Vertex_converter<Tds, RTds>(),
          Cell_converter<Tds, RTds>()));
    }

    template<typename T3, typename Remeshing_tr>
    void build_from_remeshing_triangulation(const Remeshing_tr& remeshing_tr,
                                            T3& tr)
    {
      typedef typename T3::Triangulation_data_structure Tds;
      typedef typename Remeshing_tr::Tds                RTds;

      tr.clear();
      tr.set_infinite_vertex(
        tr.tds().copy_tds(
          remeshing_tr.tds(),
          remeshing_tr.infinite_vertex(),
          Vertex_converter<RTds, Tds>(),
          Cell_converter<RTds, Tds>()));
    }

    void update_c3t3(C3t3& c3t3)
    {
      for (typename C3t3::Triangulation::Finite_facets_iterator
        fit = c3t3.triangulation().finite_facets_begin();
        fit != c3t3.triangulation().finite_facets_end();
        ++fit)
      {
        typename C3t3::Triangulation::Facet f = *fit;
        typename C3t3::Triangulation::Cell::Subdomain_index
          s1 = f.first->subdomain_index(),
          s2 = f.first->neighbor(f.second)->subdomain_index();
        if (s1 != s2)
        {
          if (s1 > s2)
            std::swap(s1, s2);
          c3t3.add_to_complex(f, s1 + 100 * s2);// std::make_pair(s1, s2));
        }
      }
      for (typename C3t3::Triangulation::Finite_cells_iterator
        cit = c3t3.triangulation().finite_cells_begin();
        cit != c3t3.triangulation().finite_cells_end();
        ++cit)
      {
        typename C3t3::Triangulation::Cell::Subdomain_index
          si = cit->subdomain_index();
        if (si != 0)
        {
          cit->set_subdomain_index(0);//o.w. add_to_complex() does nothing
          c3t3.add_to_complex(cit, si);
        }
      }
    }
  }
}

using namespace CGAL::Three;
class Polyhedron_demo_tetrahedral_remeshing_plugin :
  public QObject,
  public Polyhedron_demo_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0" FILE "tetrahedral_remeshing_plugin.json")

public:
  void init(QMainWindow* mainWindow, Scene_interface* scene_interface, Messages_interface*)
  {
    this->scene = scene_interface;
    this->mw = mainWindow;

    actionTetrahedralRemeshing_ = new QAction("Tetrahedral Remeshing", mw);
    if (actionTetrahedralRemeshing_) {
      connect(actionTetrahedralRemeshing_, SIGNAL(triggered()),
        this, SLOT(tetrahedral_remeshing()));
    }
  }

  QList<QAction*> actions() const {
    return QList<QAction*>() << actionTetrahedralRemeshing_;
  }

  bool applicable(QAction*) const
  {
    return qobject_cast<Scene_c3t3_item*>(scene->item(scene->mainSelectionIndex()));
  }

public Q_SLOTS:
  void tetrahedral_remeshing()
  {
    typedef CGAL::Tetrahedral_remeshing::Remeshing_triangulation_3<EPICK> Remeshing_triangulation;

    const Scene_interface::Item_id index = scene->mainSelectionIndex();

    Scene_c3t3_item* c3t3_item =
      qobject_cast<Scene_c3t3_item*>(scene->item(index));

    if (c3t3_item)
    {
      // Create dialog box
      QDialog dialog(mw);
      Ui::Tetrahedral_remeshing_dialog ui
        = tet_remeshing_dialog(&dialog, c3t3_item);

      // Get values
      int i = dialog.exec();
      if (i == QDialog::Rejected)
      {
        std::cout << "Remeshing aborted" << std::endl;
        return;
      }
      double target_length = ui.edgeLength_dspinbox->value();
      unsigned int nb_iter = ui.nbIterations_spinbox->value();
      bool protect = ui.protect_checkbox->isChecked();

      // wait cursor
      QApplication::setOverrideCursor(Qt::WaitCursor);

      QTime time;
      time.start();

      Remeshing_triangulation tr;
      CGAL::internal::build_remeshing_triangulation(c3t3_item->c3t3().triangulation(), tr);

      std::cout << "Remeshing triangulation built (" << time.elapsed() << " ms)" << std::endl;
      time.restart();

      CGAL::tetrahedral_adaptive_remeshing(tr, target_length,
        CGAL::parameters::protect_boundaries(protect)
        .number_of_iterations(nb_iter));

      std::cout << "Remeshing done (" << time.elapsed() << " ms)" << std::endl;
      time.restart();

      c3t3_item->c3t3().clear();
      CGAL::internal::build_from_remeshing_triangulation(tr, c3t3_item->c3t3().triangulation());
      CGAL::internal::update_c3t3(c3t3_item->c3t3());

      std::cout << "Back conversion done (" << time.elapsed() << " ms)" << std::endl;

      c3t3_item->c3t3_changed();
      this->scene->itemChanged(index);

      // default cursor
      QApplication::restoreOverrideCursor();
    }
    else
    {
      std::cout << "Can't remesh that type of thing" << std::endl;
    }
  }

private:
  Scene_interface *scene;
  QMainWindow* mw;

  Ui::Tetrahedral_remeshing_dialog
  tet_remeshing_dialog(QDialog* dialog,
                       Scene_c3t3_item* c3t3_item)
  {
    Ui::Tetrahedral_remeshing_dialog ui;
    ui.setupUi(dialog);
    connect(ui.buttonBox, SIGNAL(accepted()), dialog, SLOT(accept()));
    connect(ui.buttonBox, SIGNAL(rejected()), dialog, SLOT(reject()));

    //Set default parameters
    Scene_interface::Bbox bbox = c3t3_item->bbox();
    ui.objectName->setText(c3t3_item->name());

    ui.objectNameSize->setText(
      tr("Object bbox size (w,h,d):  <b>%1</b>,  <b>%2</b>,  <b>%3</b>")
      .arg(bbox.xmax()-bbox.xmin(), 0, 'g', 3)
      .arg(bbox.ymax()-bbox.ymin(), 0, 'g', 3)
      .arg(bbox.zmax()-bbox.zmin(), 0, 'g', 3));

    double diago_length = CGAL::sqrt((bbox.xmax()-bbox.xmin())*(bbox.xmax()-bbox.xmin())
                                   + (bbox.ymax()-bbox.ymin())*(bbox.ymax()-bbox.ymin())
                                   + (bbox.zmax()-bbox.zmin())*(bbox.zmax()-bbox.zmin()));
    double log = std::log10(diago_length);
    unsigned int nb_decimals = (log > 0) ? 5 : (std::ceil(-log)+3);

    ui.edgeLength_dspinbox->setDecimals(nb_decimals);
    ui.edgeLength_dspinbox->setSingleStep(1e-3);
    ui.edgeLength_dspinbox->setRange(1e-6 * diago_length, //min
                                     2.   * diago_length);//max
    ui.edgeLength_dspinbox->setValue(0.05 * diago_length);

    std::ostringstream oss;
    oss << "Diagonal length of the Bbox of the triangulation to remesh is ";
    oss << diago_length << "." << std::endl;
    oss << "Default is 5% of it" << std::endl;
    ui.edgeLength_dspinbox->setToolTip(QString::fromStdString(oss.str()));

    ui.nbIterations_spinbox->setSingleStep(1);
    ui.nbIterations_spinbox->setRange(1/*min*/, 1000/*max*/);
    ui.nbIterations_spinbox->setValue(1);

    ui.protect_checkbox->setChecked(false);

    return ui;
  }


private:
  QAction* actionTetrahedralRemeshing_;

}; // end Polyhedron_demo_tetrahedral_remeshing_plugin

#include "Tetrahedral_remeshing_plugin.moc"
