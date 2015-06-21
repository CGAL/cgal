#include "polyhedral_surface.h"
#include "get_polyhedral_surface.h"

#include <QMainWindow>
#include <QString>
#include <QStatusBar>
#include <QApplication>
#include <QAction>
#include <QDialog>
#include "ui_optionsdialog.h"
#include <QDoubleSpinBox>
#include <iostream>
#include <fstream>
#include <boost/format.hpp>
#include "mainwindow.h"

typedef CGAL_polyhedral_surface::Polyhedron Polyhedron;

Polyhedral_surface::Polyhedral_surface(QObject* parent,
                                       double sharp_edges_angle_lower_bound,
                                       double sharp_edges_angle_upper_bound)
  : Surface(parent), 
    m_inverse_normals(false),
    surface_ptr(0),
    parent(parent),
    display_octree(false),
    display_edges_octree(false),
    display_surface(true),
    display_all_edges(true),
    display_control_edges(false),
    sharp_edges_angle_lower_bound(sharp_edges_angle_lower_bound),
    sharp_edges_angle_upper_bound(sharp_edges_angle_upper_bound),
    is_octree_initialized(false),
    selected_edge(-1),
    selected_facet(-1),
    is_dirty(true),
    list_id(0)
{
  connection_map["actionDisplay_octree"] = 
    std::make_pair(SIGNAL(toggled(bool)),
                   SLOT(toggle_display_octree(bool)));

  connection_map["actionDisplay_edges_octree"] = 
    std::make_pair(SIGNAL(toggled(bool)),
                   SLOT(toggle_display_edges_octree(bool)));

  connection_map["actionDisplay_surface"] = 
    std::make_pair(SIGNAL(toggled(bool)),
                   SLOT(toggle_display_surface(bool)));

  connection_map["actionDisplay_all_edges"] = 
    std::make_pair(SIGNAL(toggled(bool)),
                   SLOT(toggle_display_all_edges(bool)));

  connection_map["actionDisplay_control_edges"] = 
    std::make_pair(SIGNAL(toggled(bool)),
                   SLOT(toggle_display_control_edges(bool)));

  connection_map["actionInverse_normals"] = 
    std::make_pair(SIGNAL(toggled(bool)),
                   SLOT(set_inverse_normals(bool)));

  connection_map["actionSubdivision"] = 
    std::make_pair(SIGNAL(triggered()),
                   SLOT(make_one_subdivision_step()));
  connection_map["action_Options"] =
    std::make_pair(SIGNAL(triggered()),
                   SLOT(on_action_Options_triggered()));
}

Polyhedral_surface::~Polyhedral_surface()
{
  clear();
  delete surface_ptr;
}

void Polyhedral_surface::on_action_Options_triggered()
{
  QDialog *options_dialog = new QDialog(qobject_cast<QWidget*>(parent));
  Ui::OptionDialog ui;
  ui.setupUi(options_dialog);

  QDoubleSpinBox* sb_upper = options_dialog->findChild<QDoubleSpinBox*>("angle_upper_bound");
  QDoubleSpinBox* sb_lower = options_dialog->findChild<QDoubleSpinBox*>("angle_lower_bound");

  if(!sb_lower || !sb_upper) 
    return;

  sb_lower->setValue(sharp_edges_angle_lower_bound);
  sb_upper->setValue(sharp_edges_angle_upper_bound);
  if(options_dialog->exec() == QDialog::Accepted)
  {
    sharp_edges_angle_upper_bound = sb_upper->value();
    sharp_edges_angle_lower_bound = sb_lower->value();
    set_sharp_edges_angle_bounds(sharp_edges_angle_lower_bound,
                                 sharp_edges_angle_upper_bound);
  }
}

void Polyhedral_surface::clear() {
  for(Connection_map::const_iterator
        it = connection_map.begin(),
        end = connection_map.end();
      it != end;
      ++it)
  {
    QAction* action = parent->findChild<QAction*>(it->first);
    action->setVisible(false);
  }
}

void Polyhedral_surface::connect_actions()
{
  for(Connection_map::const_iterator
        it = connection_map.begin(),
        end = connection_map.end();
      it != end;
      ++it)
  {
    QAction* action = parent->findChild<QAction*>(it->first);
    action->setVisible(true);
    if(action)
      connect(action, it->second.first,
              this, it->second.second);
  }
  MainWindow* mw = qobject_cast<MainWindow *>(parent);
  if(mw) {
//     mw->fix_menus_visibility();
    mw->show_only("polyhedral");
  }


  connect(this, SIGNAL(changed()), this, SLOT(display_nb_elements_in_status_bar()));
}

void Polyhedral_surface::display_nb_elements_in_status_bar() const 
{
  QMainWindow* mw = qobject_cast<QMainWindow *>(parent);
  if(surface_ptr && mw)
  {
    mw->statusBar()->showMessage(QString("%1 vertices. %2 edges. %3 facets.")
                                 .arg(surface_ptr->size_of_vertices())
                                 .arg(surface_ptr->size_of_halfedges()/2)
                                 .arg(surface_ptr->size_of_facets()));
  }
}

void Polyhedral_surface::set_dirty()
{
  is_dirty = true;
  Q_EMIT changed();
}

void Polyhedral_surface::busy() const 
{
  QMainWindow* mw = qobject_cast<QMainWindow *>(parent);
  if(mw)
  {
    mw->statusBar()->showMessage(QString("Constructing octree..."));
  }
  QApplication::setOverrideCursor(QCursor(Qt::WaitCursor));
}

void Polyhedral_surface::not_busy() const 
{
  QApplication::restoreOverrideCursor();
  QMainWindow* mw = qobject_cast<QMainWindow *>(parent);
  if(mw)
  {
    mw->statusBar()->clearMessage();
  }
}

void Polyhedral_surface::set_sharp_edges_angle_bounds(const double lower_bound,
                                                      const double upper_bound)
{
  sharp_edges_angle_lower_bound = lower_bound;
  sharp_edges_angle_upper_bound = upper_bound;
  if(surface_ptr) {
    surface_ptr->set_sharp_edges_angle_bounds(lower_bound, upper_bound);
    surface_ptr->set_sharp_vertices_angle_bounds(lower_bound, upper_bound);
    update_data_structures();
    Q_EMIT set_dirty();
  }
}

void Polyhedral_surface::update_data_structures() 
{
  surface_ptr->compute_sharp_edges_incidence_graph();
  if(display_octree || display_edges_octree) {
    construct_octree();
    is_octree_initialized = true;
  }
  else
    is_octree_initialized = false;
}

void Polyhedral_surface::construct_octree() 
{
  busy();
  surface_ptr->construct_octree();
  not_busy();
}

void Polyhedral_surface::toggle_display_octree(bool b)
{
  if(surface_ptr && b && !is_octree_initialized) {
    is_octree_initialized = true;
    construct_octree();
  }
  display_octree = b;
  Q_EMIT set_dirty();
}

void Polyhedral_surface::toggle_display_edges_octree(bool b)
{
  if(surface_ptr && b && !is_octree_initialized) {
    is_octree_initialized = true;
    construct_octree();
  }
  display_edges_octree = b;
  Q_EMIT set_dirty();
}

void Polyhedral_surface::toggle_display_surface(bool b)
{
  display_surface = b;
  Q_EMIT set_dirty();
}

void Polyhedral_surface::toggle_display_all_edges(bool b)
{
  display_all_edges = b;
  Q_EMIT set_dirty();
}

void Polyhedral_surface::toggle_display_control_edges(bool b)
{
  display_control_edges = b;
  Q_EMIT set_dirty();
}

void Polyhedral_surface::make_one_subdivision_step()
{
  if(surface_ptr)
  {
    Polyhedron output;
    CSubdivider_loop<Polyhedron , Poly_kernel> pw_loop_subdiviser;

    pw_loop_subdiviser.subdivide(*surface_ptr, output);
    static_cast<Polyhedron&>(*surface_ptr) = output;
    surface_ptr->compute_normals();
    surface_ptr->compute_type();
    update_data_structures();
    Q_EMIT set_dirty();
  }
}

bool Polyhedral_surface::open(const QString& filename)
{
  clear();

  std::cerr << "Opening file \"" << qPrintable(filename) << "\"...";
  std::ifstream in(filename.toUtf8());
  if(!in) return false;

  QApplication::setOverrideCursor(QCursor(Qt::WaitCursor));

  if(surface_ptr)
    delete surface_ptr;
  surface_ptr = new CGAL_polyhedral_surface(in, 
                                            sharp_edges_angle_lower_bound, 
                                            sharp_edges_angle_upper_bound,
                                            false /*do not construct 
                                                    octree*/);
  if(!in) {
    QApplication::restoreOverrideCursor();
    return false;
  }
  is_octree_initialized = false;
  selected_facet = selected_edge = -1;
  update_data_structures();

  connect_actions();
  std::cerr << " Done.\n";
  QApplication::restoreOverrideCursor();
  float xmin, ymin, zmin, xmax, ymax, zmax;
  get_bbox(xmin, ymin, zmin, xmax, ymax, zmax);
  const float xcenter = (xmin + xmax) / 2;
  const float ycenter = (ymin + ymax) / 2;
  const float zcenter = (zmin + zmax) / 2;
  const float xdelta = (-xmin + xmax);
  const float ydelta = (-ymin + ymax);
  const float zdelta = (-zmin + zmax);
//   const float radius = std::max(std::max(xdelta, ydelta), zdelta) * std::sqrt(3.)/ 2.;
  std::cerr << boost::format("Bounding box: xmin=%1%, ymin=%2%, zmin=%3%\n"
                             "              xmax=%4%, ymax=%5%, zmax=%6%\n"
                             "              center=(%7%, %8%, %9%)\n")
    % xmin % ymin % zmin % xmax % ymax % zmax
    % xcenter % ycenter % zcenter
            << boost::format("              span=(%1%,%2%,%3%)\n")
    % xdelta % ydelta % zdelta;
  viewer->camera()->setSceneBoundingBox(qglviewer::Vec(xmin, ymin, zmin),
                                        qglviewer::Vec(xmax, ymax, zmax));
  viewer->setBackgroundColor(Qt::white);
  viewer->showEntireScene();
  
  QAction* actionInverse_normals = qFindChild<QAction*>(this, "actionInverse_normals");
  if(actionInverse_normals) actionInverse_normals->setChecked(false);
  Q_EMIT set_dirty();
  return true;
}

void Polyhedral_surface::close() 
{
  delete surface_ptr;
  surface_ptr = 0;
}

void Polyhedral_surface::draw() {
  draw(false);
}

void Polyhedral_surface::drawWithNames() {
  draw(true);
}

void Polyhedral_surface::postSelection(const QPoint&)
{
  if(!surface_ptr) return;

  selected_facet = selected_edge = -1;
    
  const int nb_vertices = surface_ptr->incidence_graph.vertices.size();
  const int nb_edges = surface_ptr->incidence_graph.edges.size();
  if(viewer->selectedName() >= nb_edges + nb_vertices)
    selected_facet = viewer->selectedName() - nb_edges - nb_vertices;
  else if(viewer->selectedName() >= nb_vertices)
    selected_edge = viewer->selectedName() - nb_vertices;

  std::cerr << boost::format("post-selection.\n"
                             "selectedName()=%1%\n"
                             "selected edge=%2%\n"
                             "selected facet=%3%\n")
    % viewer->selectedName() % selected_edge % selected_facet;

  Q_EMIT set_dirty();
}
  
void Polyhedral_surface::draw(bool with_names)
{
  if(!list_id)
  {
    std::cerr << "Generating OpenGL display list ID: ";
    std::cerr << (list_id = ::glGenLists(1)) << "\n";
  }
  if(!with_names)
  {
    if(is_dirty)
      ::glNewList(list_id, GL_COMPILE_AND_EXECUTE);
    else if(::glIsList(list_id))
    {
      ::glCallList(list_id);
      return;
    }
    else
      std::cerr << "Call list (" << list_id << ")failed.\n";
  }

  if(surface_ptr)
  {
    if(display_surface)
    {
      // enable polygon offset
      ::glEnable(GL_POLYGON_OFFSET_FILL);
      ::glPolygonOffset(1.0f,1.0f);
      ::glEnable(GL_LIGHTING);
      ::glColor3f(0.2f, 0.2f, 1.f);
      if(with_names)
        surface_ptr->gl_draw_direct_triangles_with_name(false,
                                                        true,
                                                        inverse_normals());
      else
        surface_ptr->gl_draw_almost_all_triangles(selected_facet,
                                                  false,
                                                  true,
                                                  inverse_normals());
      if(!with_names && selected_facet >= 0)
      {
        ::glColor3f(1., 1.f, 0.f);
        surface_ptr->incidence_graph.gl_draw_facet(selected_facet,
                                                   false,
                                                   true,
                                                   inverse_normals());
      }
      ::glDisable(GL_LIGHTING);
      ::glLineWidth(1.0f);
      if(display_all_edges)
      {
        // superimpose ordinary edges
        ::glColor3d(0.,0.,.8);
        surface_ptr->superimpose_edges(false,display_control_edges);
      }
      // superimpose control edges
      if(display_control_edges)
      {
        ::glDisable(GL_LIGHTING);
        ::glColor3d(.0, .0, .0);
        ::glLineWidth(1.0f);
        surface_ptr->superimpose_edges(true,false);
      }

      // draw sharp edges
      ::glColor3ub(128,128,128);
      if(with_names)
        surface_ptr->gl_draw_sharp_edges_with_names(3.0f,255,0,0);
      else
        surface_ptr->gl_draw_sharp_edges(3.0f,255,0,0);
      if(!with_names && selected_edge >= 0)
      {
        ::glLineWidth(3.0f);
        ::glColor3d(0., 1., 0.);
        surface_ptr->incidence_graph.gl_draw_edge(selected_edge);
      }
    } // end if display_surface
    if(!with_names && (display_octree||display_edges_octree) )
    {
      ::glColor3ub(0,0,0);
      ::glLineWidth(1.0f);
      ::glDisable(GL_LIGHTING);
      ::glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
      ::glEnable(GL_LINE_STIPPLE);
      if(display_octree)
      {
        ::glColor3ub(0,0,0);
        surface_ptr->gl_draw_facet_octree();
      }
      if(display_edges_octree)
      {
        ::glColor3d(1.,0.,0.);
        surface_ptr->gl_draw_edges_octree();
      }
      ::glDisable(GL_LINE_STIPPLE);
      ::glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
      ::glEnable(GL_LIGHTING);
    } // end if display_octree
      
    ::glDisable(GL_POLYGON_OFFSET_FILL);
      
    if(!with_names && is_dirty)
    {
      ::glEndList();
      is_dirty = false;
    }
  } // end if(surface_ptr)
}
void Polyhedral_surface::get_bbox(float& xmin, float& ymin, float& zmin,
                                  float& xmax, float& ymax, float& zmax)
{
  if(surface_ptr) {
    xmin=surface_ptr->bbox().xmin();
    ymin=surface_ptr->bbox().ymin();
    ymin=surface_ptr->bbox().zmin();
    xmax=surface_ptr->bbox().xmax();
    ymax=surface_ptr->bbox().ymax();
    zmax=surface_ptr->bbox().zmax();
  }
  else 
  {
    xmin = ymin = zmin = 0.f;
    xmax = ymax = zmax = 1.f;
  }
}

void Polyhedral_surface::set_inverse_normals(const bool b) {
  m_inverse_normals = b;
  set_dirty();
}

bool Polyhedral_surface::inverse_normals() const {
  return m_inverse_normals;
}


Surface* get_polyhedral_surface(QObject* parent,
				double sharp_edges_angle_lower_bound,
				double sharp_edges_angle_upper_bound = 180.)
{
  return new Polyhedral_surface(parent,
				sharp_edges_angle_lower_bound,
				sharp_edges_angle_upper_bound);
}

#include "polyhedral_surface.moc"
