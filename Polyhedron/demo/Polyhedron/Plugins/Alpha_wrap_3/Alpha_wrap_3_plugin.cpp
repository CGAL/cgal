#include <QtCore/qglobal.h>

#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>

#include "Scene_surface_mesh_item.h"
#include "Scene_polygon_soup_item.h"
#include "Scene_polyhedron_selection_item.h"
#include "Scene_polylines_item.h"
#include "Scene_points_with_normal_item.h"

// Since we want to do visualization and interruption, it's better to use the sorted priority queue,
// even if it is slower
#define CGAL_AW3_USE_SORTED_PRIORITY_QUEUE

#include <CGAL/alpha_wrap_3.h>

#include <QAction>
#include <QApplication>
#include <QDialog>
#include <QElapsedTimer>
#include <QMainWindow>
#include <QMessageBox>
#include <QTextStream>
#include <QString>
#include <QThread>
#include <QTimer>
#include <QTranslator>
#include <QtPlugin>

#include <vector>
#include <algorithm>
#include <queue>
#include <sstream>
#include <stdexcept>

#include "ui_alpha_wrap_3_dialog.h"

using TS_Oracle = CGAL::Alpha_wraps_3::internal::Triangle_soup_oracle<Kernel>;
using SS_Oracle = CGAL::Alpha_wraps_3::internal::Segment_soup_oracle<Kernel, TS_Oracle>;
using Oracle = CGAL::Alpha_wraps_3::internal::Point_set_oracle<Kernel, SS_Oracle>;
using Wrapper = CGAL::Alpha_wraps_3::internal::Alpha_wrapper_3<Oracle>;

// Here is the pipeline for the interruption box:
// - The main window is connected to a wrapping thread, which performs the wrapping.
// - The wrapping has a visitor, AW3_interrupter_visitor, which has a shared_ptr to a Boolean
// - When the user clicks the box, the Boolean is switched to *false*
// - The wrapping thread creates the wip mesh

// Here is the pipeline for the iterative visualization:
// - The main window is connected to a wrapping thread, which performs the wrapping.
// - The wrapping has a visitor, Iterative_AW3_visualization_visitor
// - The visitor has a shared pointer to an emiter (can't emit directly from the visitor)
// - The visitor has shared pointers to a polygon soup (+ colors), which it updates itself
//   before emitting a signal
// - There is a pause in the emition because it needs to wait for the main thread to draw the
//   polygon soup before the visitor updates the polygon soup.

struct Iterative_update_emiter
  : public QObject
{
  Q_OBJECT

public:
  void emit_new_iteration(int sid)
  {
    Q_EMIT new_iteration(sid);
    CGAL::Three::Three::getMutex()->lock();
    Three::getWaitCondition()->wait(CGAL::Three::Three::getMutex());
    CGAL::Three::Three::getMutex()->unlock();
  }

  void emit_last_iteration(int sid)
  {
    // Last iteration only updates the (existing) soup item's properties, but there is no change
    // in geometry, so there is no need to wait for the main thread to update the main window.
    Q_EMIT last_iteration(sid);
  }

Q_SIGNALS:
  void new_iteration(int);
  void last_iteration(int sid);
};

struct Iterative_AW3_visualization_visitor
  : public CGAL::Alpha_wraps_3::internal::Wrapping_default_visitor
{
private:
  const bool visualize_iterations;

  std::shared_ptr<std::vector<Kernel::Point_3> > points;
  std::shared_ptr<std::vector<std::vector<std::size_t> > > faces;
  std::shared_ptr<std::vector<CGAL::IO::Color> > fcolors;
  std::shared_ptr<std::vector<CGAL::IO::Color> > vcolors;
  std::shared_ptr<Iterative_update_emiter> emiter;
  int sid = 0;

public:
  Iterative_AW3_visualization_visitor(const bool visualize_iterations,
                                      std::shared_ptr<std::vector<Kernel::Point_3> > points,
                                      std::shared_ptr<std::vector<std::vector<std::size_t> > > faces,
                                      std::shared_ptr<std::vector<CGAL::IO::Color> > fcolors,
                                      std::shared_ptr<std::vector<CGAL::IO::Color> > vcolors,
                                      std::shared_ptr<Iterative_update_emiter> emiter)
    : visualize_iterations(visualize_iterations),
      points(points), faces(faces), fcolors(fcolors), vcolors(vcolors),
      emiter(emiter)
  { }

  template <typename AlphaWrapper, typename Point>
  void before_Steiner_point_insertion(const AlphaWrapper& wrapper,
                                      const Point& /* p */)
  {
    if(!visualize_iterations)
      return;

    if(!points || !faces || !fcolors || !vcolors)
      return;

    // If the next top of the queue has vertices on the bbox, don't draw (try to avoid producing
    // spikes in the visualization)
//    const auto& gate = wrapper.queue().top();
//    if(wrapper.triangulation().number_of_vertices() > 500 && gate.is_permissive_facet())
//      return;

    // Skip some...
    if(wrapper.triangulation().number_of_vertices() % 50 != 0)
      return;

    // Extract the wrap as a triangle soup
    points->clear();
    faces->clear();
    fcolors->clear();
    vcolors->clear();

    using Dt = typename std::decay<decltype(wrapper.triangulation())>::type;
    using Vertex_handle = typename Dt::Vertex_handle;
    using Facet = typename Dt::Facet;
    using Cell_handle = typename Dt::Cell_handle;

    std::unordered_map<Vertex_handle, std::size_t> vertex_to_id;
    std::size_t nv = 0;

#if 0
    // This is used to compute colors depending on what is old and what is new.
    // It is not currently used (a uniform gray color is used), but leaving it as it might be useful.
    std::size_t min_time_stamp = -1, max_time_stamp = 0;
    for(auto cit=wrapper.triangulation().finite_cells_begin(), cend=wrapper.triangulation().finite_cells_end(); cit!=cend; ++cit)
    {
      if(cit->time_stamp() > max_time_stamp)
        max_time_stamp = cit->time_stamp();
      if(cit->time_stamp() < min_time_stamp)
        min_time_stamp = cit->time_stamp();
    }
#endif

    for(auto fit=wrapper.triangulation().finite_facets_begin(), fend=wrapper.triangulation().finite_facets_end(); fit!=fend; ++fit)
    {
      Facet f = *fit;
      if(!f.first->is_outside())
        f = wrapper.triangulation().mirror_facet(f);

      const Cell_handle c = f.first;
      const int s = f.second;
      const Cell_handle nh = c->neighbor(s);
      if(c->is_outside() == nh->is_outside())
        continue;

      std::array<std::size_t, 3> ids;
      for(int pos=0; pos<3; ++pos)
      {
        Vertex_handle vh = c->vertex(Dt::vertex_triple_index(s, pos));
        auto insertion_res = vertex_to_id.emplace(vh, nv);
        if(insertion_res.second) // successful insertion, never-seen-before vertex
        {
          points->push_back(wrapper.triangulation().point(vh));
          vcolors->push_back(CGAL::IO::Color(0, 0, 0));
          ++nv;
        }

        ids[pos] = insertion_res.first->second;
      }

      faces->emplace_back(std::vector<std::size_t>{ids[0], ids[1], ids[2]});

#if 0
      double color_val = double(c->time_stamp() - min_time_stamp) / double(max_time_stamp - min_time_stamp);
       color_val = int(256. * color_val);
      fcolors.push_back(CGAL::IO::Color(color_val, 10, 150)); // young is red, old is blue
      // fcolors.push_back(CGAL::IO::Color(256 - color_val, 256 - color_val, 256 - color_val)); // young is light, old is dark
#endif
      fcolors->push_back(CGAL::IO::Color(100, 100, 100)); // uniform darkish gray
    }

    emiter->emit_new_iteration(sid++);
  }

  template <typename AlphaWrapper>
  void on_alpha_wrapping_end(const AlphaWrapper&)
  {
    if(!visualize_iterations)
      return;

    emiter->emit_last_iteration(sid);
  }
};

template <typename BaseVisitor>
struct AW3_interrupter_visitor
  : BaseVisitor
{
  // shared pointer because visitors are copied
  std::shared_ptr<bool> should_stop = std::make_shared<bool>(false);

  AW3_interrupter_visitor(const BaseVisitor base)
    : BaseVisitor(base)
  { }

  template <typename Wrapper>
  constexpr bool go_further(const Wrapper&)
  {
    return !(*should_stop);
  }
};

struct Wrapper_thread
  : public QThread
{
  Q_OBJECT

  using Visitor = AW3_interrupter_visitor<Iterative_AW3_visualization_visitor>;

public:
  Wrapper wrapper;
  const Oracle oracle;
  const double alpha, offset;
  const bool enforce_manifoldness;
  Visitor visitor;

  SMesh wrap;

  QTimer* timer;

public:
  Wrapper_thread(const Oracle& oracle,
                 const double alpha,
                 const double offset,
                 const bool enforce_manifoldness,
                 Visitor visitor)
    : wrapper(oracle),
      alpha(alpha), offset(offset),
      enforce_manifoldness(enforce_manifoldness),
      visitor(visitor),
      timer(new QTimer(this))
  {
    connect(timer, SIGNAL(timeout()),
            this,  SLOT(emit_status()));

    timer->start(1000);
  }

  ~Wrapper_thread()
  {
    delete timer;
  }

  void run() override
  {
    QElapsedTimer elapsed_timer;
    elapsed_timer.start();

    wrapper(alpha, offset, wrap,
            CGAL::parameters::do_enforce_manifoldness(enforce_manifoldness)
                             .visitor(visitor));

    if(wrapper.queue().empty())
      Q_EMIT done(this);
    else
      Q_EMIT interrupted(this);

    std::cout << "Wrapping took " << elapsed_timer.elapsed() / 1000. << "s" << std::endl;
  }

public Q_SLOTS:
  void emit_status()
  {
    Q_EMIT status_report(QString("%1 vertices").arg(wrapper.triangulation().number_of_vertices()));
  }

Q_SIGNALS:
  void interrupted(Wrapper_thread*);
  void done(Wrapper_thread*);
  void status_report(QString);
};

class Polyhedron_demo_alpha_wrap_3_plugin
  : public QObject,
    public CGAL::Three::Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

  using Triangles = std::vector<Kernel::Triangle_3>;
  using Segments = std::vector<Kernel::Segment_3>;
  using Points = std::vector<Kernel::Point_3>;

private:
  CGAL::Bbox_3 m_wrap_bbox;
  double m_wrap_bbox_diag_length;

  QAction* actionAlpha_wrap_3_ = nullptr;
  Ui::alpha_wrap_3_dialog ui;

  // GUI for the interruption
  QMessageBox* m_message_box = nullptr;

  // storage of intermediate wraps for iterative visualization
  std::shared_ptr<std::vector<Kernel::Point_3> > m_iter_points;
  std::shared_ptr<std::vector<std::vector<std::size_t> > > m_iter_faces;
  std::shared_ptr<std::vector<CGAL::IO::Color> > m_iter_fcolors;
  std::shared_ptr<std::vector<CGAL::IO::Color> > m_iter_vcolors;
  Scene_polygon_soup_item* m_iterative_wrap_item = nullptr;
  bool m_do_snapshot = false;

public:
  void init(QMainWindow* mainWindow,
            CGAL::Three::Scene_interface* scene_interface,
            Messages_interface*)
  {
    this->scene = scene_interface;
    this->mw = mainWindow;

    actionAlpha_wrap_3_ = new QAction("3D Alpha Wrapping", this->mw);
    if(actionAlpha_wrap_3_)
    {
      connect(actionAlpha_wrap_3_, SIGNAL(triggered()),
              this,                SLOT(on_actionAlpha_wrap_3_triggered()));
    }
  }

  bool applicable(QAction*) const
  {
    if(scene->selectionIndices().empty())
      return false;

    for(int index : scene->selectionIndices())
    {
      if(!qobject_cast<Scene_polygon_soup_item*>(scene->item(index)) &&
         !qobject_cast<Scene_surface_mesh_item*>(scene->item(index)) &&
         !qobject_cast<Scene_polyhedron_selection_item*>(scene->item(index)) &&
         !qobject_cast<Scene_polylines_item*>(scene->item(index)) &&
         !qobject_cast<Scene_points_with_normal_item*>(scene->item(index)))
      {
        return false;
      }
    }

    return true;
  }

  QList<QAction*> actions() const
  {
    return QList<QAction*>() << actionAlpha_wrap_3_;
  }

private:
  void print_message(QString message) const
  {
    CGAL::Three::Three::information(message);
  }

public Q_SLOTS:
  // This is UI stuff
  void on_alphaValue_changed(double)
  {
    QSignalBlocker block(ui.relativeAlphaValue);
    ui.relativeAlphaValue->setValue(m_wrap_bbox_diag_length / ui.alphaValue->value());
  }

  void on_relativeAlphaValue_changed(double)
  {
    QSignalBlocker block(ui.alphaValue);
    ui.alphaValue->setValue(m_wrap_bbox_diag_length / ui.relativeAlphaValue->value());
  }

  void on_offsetValue_changed(double)
  {
    QSignalBlocker block(ui.relativeOffsetValue);
    ui.relativeOffsetValue->setValue(m_wrap_bbox_diag_length / ui.offsetValue->value());
  }

  void on_relativeOffsetValue_changed(double)
  {
    QSignalBlocker block(ui.offsetValue);
    ui.offsetValue->setValue(m_wrap_bbox_diag_length / ui.relativeOffsetValue->value());
  }

  void update_iteration_snapshot_checkbox()
  {
    ui.snapshotIterations->setCheckable(ui.visualizeIterations->isChecked());
  }

  // This is for the visualization
  void update_iterative_wrap_item(int sid)
  {
    if(m_iterative_wrap_item == nullptr)
      return;

    // Update the wrap item's visualization
    m_iterative_wrap_item->load(*m_iter_points, *m_iter_faces, *m_iter_fcolors, *m_iter_vcolors);
    m_iterative_wrap_item->setName(QString("Iterative wrap #%1").arg(sid));
    m_iterative_wrap_item->setAlpha(255 / 2);
    m_iterative_wrap_item->setRenderingMode(FlatPlusEdges);

    m_iterative_wrap_item->invalidateOpenGLBuffers();
    m_iterative_wrap_item->redraw();
    m_iterative_wrap_item->itemChanged();

    // Refresh the view
    CGAL::Three::Viewer_interface* viewer = CGAL::Three::Three::activeViewer();
    viewer->update();

    CGAL::Three::Three::getWaitCondition()->wakeAll();

    if(m_do_snapshot)
    {
      std::stringstream oss;
      oss << "Wrap_iteration-" << sid << ".png" << std::ends;
      QString filename = QString::fromStdString(oss.str().c_str());

      viewer->saveSnapshot(filename, 1920, 1080, true /*expand*/, 2.0 /*oversampling*/);
    }
  }

  void finish_iterative_wrap_item(int sid)
  {
    if(m_iterative_wrap_item == nullptr)
      return;

    if(m_do_snapshot)
    {
      m_iterative_wrap_item->setName(QString("Iterative wrap #%1").arg(sid));

      m_iterative_wrap_item->setAlpha(255);
      m_iterative_wrap_item->invalidateOpenGLBuffers();
      m_iterative_wrap_item->redraw();
      m_iterative_wrap_item->itemChanged();

      // Refresh the view
      CGAL::Three::Viewer_interface* viewer = CGAL::Three::Three::activeViewer();
      viewer->update();

      std::stringstream oss;
      oss << "Wrap_iteration-" << sid << ".png" << std::ends;
      QString filename = QString::fromStdString(oss.str().c_str());

      viewer->saveSnapshot(filename);
    }

    CGAL_assertion(m_iterative_wrap_item);
    scene->erase(scene->item_id(m_iterative_wrap_item));
  }

  void reset_iterative_wrap_item()
  {
    m_iterative_wrap_item = nullptr;
  }

  // This is for the message box and thread interruption
  void wrapping_done(Wrapper_thread* wrapper_thread)
  {
    Scene_surface_mesh_item* wrap_item = new Scene_surface_mesh_item(std::move(wrapper_thread->wrap));
    wrap_item->setName(tr("Wrap with alpha %2 offset %3").arg(wrapper_thread->alpha)
                                                         .arg(wrapper_thread->offset));
    wrap_item->setColor(Qt::gray);
    const int wrap_item_id = scene->addItem(wrap_item);
    scene->setSelectedItem(wrap_item_id);

    wrapper_thread->terminate();
    wrapper_thread->wait();
    delete wrapper_thread;

    if(m_message_box)
    {
      m_message_box->done(0);
      m_message_box = nullptr;
    }
  }

  // In case we wish to do something more one day
  void wrapping_interrupted(Wrapper_thread* wrapper_thread)
  {
    wrapping_done(wrapper_thread);
  }

  void status_report(const QString& msg)
  {
    if(m_message_box == nullptr)
      return;

    m_message_box->setInformativeText(msg);
  }

  // Main call
  void on_actionAlpha_wrap_3_triggered()
  {
    QDialog dialog(mw);

    ui.setupUi(&dialog);

    connect(ui.alphaValue, SIGNAL(valueChanged(double)), this, SLOT(on_alphaValue_changed(double)));
    connect(ui.relativeAlphaValue, SIGNAL(valueChanged(double)), this, SLOT(on_relativeAlphaValue_changed(double)));
    connect(ui.offsetValue, SIGNAL(valueChanged(double)), this, SLOT(on_offsetValue_changed(double)));
    connect(ui.relativeOffsetValue, SIGNAL(valueChanged(double)), this, SLOT(on_relativeOffsetValue_changed(double)));

    connect(ui.visualizeIterations, SIGNAL(clicked(bool)), this, SLOT(update_iteration_snapshot_checkbox()));

    connect(ui.buttonBox, SIGNAL(accepted()), &dialog, SLOT(accept()));
    connect(ui.buttonBox, SIGNAL(rejected()), &dialog, SLOT(reject()));

    dialog.setWindowFlags(Qt::Dialog|Qt::CustomizeWindowHint|Qt::WindowCloseButtonHint);

    Triangles triangles;
    Segments segments;
    Points points;

    for(int index : this->scene->selectionIndices())
    {
      // ---
      Scene_surface_mesh_item* sm_item = qobject_cast<Scene_surface_mesh_item*>(this->scene->item(index));
      if(sm_item != nullptr)
      {
        SMesh* pMesh = sm_item->polyhedron();
        auto vpm = get(CGAL::vertex_point, *pMesh);

        triangles.reserve(triangles.size() + num_faces(*pMesh));
        for(boost::graph_traits<SMesh>::face_descriptor f : faces(*pMesh))
        {
          boost::graph_traits<SMesh>::halfedge_descriptor h = halfedge(f, *pMesh);
          if(!is_triangle(h, *pMesh))
          {
            print_message("Warning: a non-triangular face in input has been ignored");
            continue;
          }

          triangles.emplace_back(get(vpm, target(h, *pMesh)),
                                 get(vpm, target(next(h, *pMesh), *pMesh)),
                                 get(vpm, source(h, *pMesh)));
        }

        continue;
      }

      // ---
      Scene_polygon_soup_item* soup_item = qobject_cast<Scene_polygon_soup_item*>(scene->item(index));
      if(soup_item != nullptr)
      {
        triangles.reserve(triangles.size() + soup_item->polygons().size());

        for(const auto& p : soup_item->polygons())
        {
          if(p.size() != 3)
          {
            print_message("Warning: a non-triangular face in input has been ignored");
            continue;
          }

          triangles.emplace_back(soup_item->points()[p[0]],
                                 soup_item->points()[p[1]],
                                 soup_item->points()[p[2]]);
        }

        continue;
      }

      // ---
      Scene_polyhedron_selection_item* selection_item =
          qobject_cast<Scene_polyhedron_selection_item*>(scene->item(index));
      if(selection_item != nullptr)
      {
        SMesh* pMesh = selection_item->polyhedron();
        auto vpm = get(CGAL::vertex_point, *pMesh);

        triangles.reserve(triangles.size() + selection_item->selected_facets.size());
        for(const auto& f : selection_item->selected_facets)
        {
          boost::graph_traits<SMesh>::halfedge_descriptor h = halfedge(f, *pMesh);
          if(!is_triangle(h, *pMesh))
          {
            print_message("Warning: a non-triangular face in input has been ignored");
            continue;
          }

          triangles.emplace_back(get(vpm, target(h, *pMesh)),
                                 get(vpm, target(next(h, *pMesh), *pMesh)),
                                 get(vpm, source(h, *pMesh)));
        }

        segments.reserve(segments.size() + selection_item->selected_edges.size());
        for(const auto& e : selection_item->selected_edges)
        {
          segments.emplace_back(get(vpm, target(halfedge(e, *pMesh), *pMesh)),
                                get(vpm, target(opposite(halfedge(e, *pMesh), *pMesh), *pMesh)));
        }

        points.reserve(points.size() + selection_item->selected_vertices.size());
        for(const auto& v : selection_item->selected_vertices)
        {
          points.push_back(get(vpm, v));
        }

        continue;
      }

      // ---
      Scene_polylines_item* lines_item = qobject_cast<Scene_polylines_item*>(scene->item(index));
      if(lines_item != nullptr)
      {
        for(auto& polyline : lines_item->polylines)
        {
          for(auto it=std::begin(polyline), end=std::end(polyline); it!=end; ++it)
          {
            auto nit = std::next(it);
            if(nit != end)
              segments.emplace_back(*it, *nit);
          }
        }

        continue;
      }

      // ---
      Scene_points_with_normal_item* pts_item = qobject_cast<Scene_points_with_normal_item*>(scene->item(index));
      if(pts_item != nullptr)
      {
        points.insert(std::cend(points),
                      std::cbegin(pts_item->point_set()->points()),
                      std::cend(pts_item->point_set()->points()));

        continue;
      }
    }

    std::cout << triangles.size() << " triangles" << std::endl;
    std::cout << segments.size() << " edges" << std::endl;
    std::cout << points.size() << " points" << std::endl;

    if(!triangles.empty())
      m_wrap_bbox = triangles.front().bbox();
    else if(!segments.empty())
      m_wrap_bbox = segments.front().bbox();
    else if(!points.empty())
      m_wrap_bbox = points.front().bbox();

    for(const Kernel::Triangle_3& tr : triangles)
      m_wrap_bbox += tr.bbox();
    for(const Kernel::Segment_3& sg : segments)
      m_wrap_bbox += sg.bbox();
    for(const Kernel::Point_3& pt : points)
      m_wrap_bbox += pt.bbox();

    std::cout << "Bbox:\n" << m_wrap_bbox << std::endl;

    // The relative value uses the bbox of the full scene and not that of selected items to wrap
    // This is intentional, both because it's tedious to make it otherwise, and because it seems
    // to be simpler to compare between "all wrapped" / "some wrapped"
    m_wrap_bbox_diag_length = std::sqrt(CGAL::square(m_wrap_bbox.xmax() - m_wrap_bbox.xmin()) +
                                        CGAL::square(m_wrap_bbox.ymax() - m_wrap_bbox.ymin()) +
                                        CGAL::square(m_wrap_bbox.zmax() - m_wrap_bbox.zmin()));

    ui.relativeAlphaValue->setValue(20.);
    ui.relativeOffsetValue->setValue(600.);
    ui.alphaValue->setValue(m_wrap_bbox_diag_length / ui.relativeAlphaValue->value());
    ui.offsetValue->setValue(m_wrap_bbox_diag_length / ui.relativeOffsetValue->value());

    // EXECUTION
    int i = dialog.exec();
    if(i == QDialog::Rejected)
      return;

    QApplication::setOverrideCursor(Qt::WaitCursor);

    const bool wrap_triangles = ui.wrapTriangles->isChecked();
    const bool wrap_segments = ui.wrapSegments->isChecked();
    const bool wrap_points = ui.wrapPoints->isChecked();

    double alpha = ui.alphaValue->value();
    double offset = ui.offsetValue->value();

    const bool enforce_manifoldness = ui.runManifoldness->isChecked();
    const bool visualize_iterations = ui.visualizeIterations->isChecked();
    m_do_snapshot = ui.snapshotIterations->isChecked();

    const bool use_message_box = ui.enableMessageBox->isChecked();

    std::cout << "Wrapping faces? " << std::boolalpha << wrap_triangles << std::endl;
    std::cout << "Wrapping edges? " << std::boolalpha << wrap_segments << std::endl;

    if(!wrap_triangles)
    {
      if(wrap_segments) // add faces' edges
      {
        segments.reserve(segments.size() + 3 * triangles.size());
        for(const auto& tr : triangles)
        {
          segments.emplace_back(tr[0], tr[1]);
          segments.emplace_back(tr[1], tr[2]);
          segments.emplace_back(tr[2], tr[0]);
        }
      }
      else // triangles & segments are not wrapped, but points are -> add faces' vertices
      {
        points.reserve(points.size() + 2 * segments.size() + 3 * triangles.size());
        for(const auto& s : segments)
        {
          points.push_back(s[0]);
          points.push_back(s[1]);
        }

        for(const auto& tr : triangles)
        {
          points.push_back(tr[0]);
          points.push_back(tr[1]);
          points.push_back(tr[2]);
        }
      }
    }

    TS_Oracle ts_oracle;
    SS_Oracle ss_oracle(ts_oracle);
    Oracle oracle(ss_oracle);

    if(wrap_triangles)
      oracle.add_triangle_soup(triangles);
    if(wrap_segments)
      oracle.add_segment_soup(segments);
    if(wrap_points)
      oracle.add_point_set(points);

    if(!oracle.do_call())
    {
      print_message("Warning: empty input - nothing to wrap");
      QApplication::restoreOverrideCursor();
      return;
    }

    if(alpha <= 0. || offset <= 0.)
    {
      print_message("Warning: alpha/offset must be strictly positive");
      QApplication::restoreOverrideCursor();
      return;
    }

    for(int index : this->scene->selectionIndices())
    {
      Scene_surface_mesh_item* sm_item = qobject_cast<Scene_surface_mesh_item*>(this->scene->item(index));
      if(sm_item != nullptr)
        sm_item->setRenderingMode(Flat);
      Scene_polygon_soup_item* soup_item = qobject_cast<Scene_polygon_soup_item*>(scene->item(index));
      if(soup_item != nullptr)
        soup_item->setRenderingMode(Flat);
    }

    if(visualize_iterations)
    {
      m_iterative_wrap_item = new Scene_polygon_soup_item();
      m_iterative_wrap_item->setName(QString("Iterative wrap"));
      const int iterative_wrap_item_id = scene->addItem(m_iterative_wrap_item);
      scene->setSelectedItem(iterative_wrap_item_id);

      // Deal with independent (e.g. manual from the main window) destruction of the iterative item
      connect(m_iterative_wrap_item, SIGNAL(aboutToBeDestroyed()),
              this, SLOT(reset_iterative_wrap_item()));
    }

    // Visitors
    m_iter_points = std::make_shared<std::vector<Kernel::Point_3> >();
    m_iter_faces = std::make_shared<std::vector<std::vector<std::size_t> > >();
    m_iter_fcolors = std::make_shared<std::vector<CGAL::IO::Color> >();
    m_iter_vcolors = std::make_shared<std::vector<CGAL::IO::Color> >();
    std::shared_ptr<Iterative_update_emiter> emiter = std::make_shared<Iterative_update_emiter>();
    Iterative_AW3_visualization_visitor visu_visitor(visualize_iterations, m_iter_points, m_iter_faces, m_iter_fcolors, m_iter_vcolors, emiter);
    AW3_interrupter_visitor<Iterative_AW3_visualization_visitor> visitor(visu_visitor);

    connect(emiter.get(), SIGNAL(new_iteration(int)),
            this,         SLOT(update_iterative_wrap_item(int)));
    connect(emiter.get(), SIGNAL(last_iteration(int)),
            this,         SLOT(finish_iterative_wrap_item(int)));

    Wrapper_thread* wrapper_thread = new Wrapper_thread(oracle, alpha, offset, enforce_manifoldness, visitor);
    if(wrapper_thread == nullptr)
    {
      QMessageBox::critical(mw, tr(""), tr("ERROR: failed to create thread"));
      return;
    }

    // Connect main thread to wrapping thread
    QObject::connect(wrapper_thread, SIGNAL(done(Wrapper_thread*)),
                     this,           SLOT(wrapping_done(Wrapper_thread*)));
    QObject::connect(wrapper_thread, SIGNAL(interrupted(Wrapper_thread*)),
                     this,           SLOT(wrapping_interrupted(Wrapper_thread*)));
    QObject::connect(wrapper_thread, SIGNAL(status_report(QString)),
                     this,           SLOT(status_report(QString)));

    // Launch thread
    CGAL::Three::Three::getMutex()->lock();
    CGAL::Three::Three::isLocked() = true;
    CGAL::Three::Three::getMutex()->unlock();

    // Create message box with stop button
    if(use_message_box)
    {
      // Switch from 'wait' to 'busy'
      QApplication::restoreOverrideCursor();
      QApplication::setOverrideCursor(Qt::BusyCursor);

      m_message_box = new QMessageBox(QMessageBox::NoIcon,
                                     "Wrapping",
                                     "Wrapping in progress...",
                                     QMessageBox::Cancel,
                                     mw);
      m_message_box->setDefaultButton(QMessageBox::Cancel);
      QAbstractButton* cancelButton = m_message_box->button(QMessageBox::Cancel);
      cancelButton->setText(tr("Stop"));

      // Connect the message box to the thread
      connect(cancelButton, &QAbstractButton::clicked,
              this, [wrapper_thread]() { *(wrapper_thread->visitor.should_stop) = true; });

      m_message_box->open();
    }

    // Actual start
    QApplication::setOverrideCursor(Qt::WaitCursor);

    wrapper_thread->start();

    CGAL::Three::Three::getMutex()->lock();
    CGAL::Three::Three::isLocked() = false;
    CGAL::Three::Three::getMutex()->unlock();

    QApplication::restoreOverrideCursor();
  }
};

#include "Alpha_wrap_3_plugin.moc"
