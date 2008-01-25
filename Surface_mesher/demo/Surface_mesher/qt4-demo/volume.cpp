#include <CGAL/basic.h>
#include "volume.h"
#include "viewer.h"
#include <CGAL/Bbox_3.h>

#include <QApplication>
#include <QAction>
#include <QMainWindow>
#include <QStatusBar>
#include <QDoubleSpinBox>

#include <CGAL/Timer.h>

Volume::Volume(QObject* parent) : 
  Surface(parent),
  m_sm_angle(30),
  m_sm_radius(0),
  m_sm_distance(0),
  m_isovalue(0),
  m_relative_precision(0.0001),
  m_view_surface(false),
  m_view_mc(false),
  parent(qobject_cast<QMainWindow*>(parent)),
  m_inverse_normals(false),
  two_sides(false)
{
  QAction* marching_cube_action = parent->findChild<QAction*>("actionMarching_cubes");
  QAction* surface_mesher_action = parent->findChild<QAction*>("actionSurface_mesher");
  spinBox_isovalue = parent->findChild<QDoubleSpinBox*>("spinBox_isovalue");
  spinBox_radius_bound = parent->findChild<QDoubleSpinBox*>("spinBox_radius_bound");
  spinBox_distance_bound = parent->findChild<QDoubleSpinBox*>("spinBox_distance_bound");

  if(spinBox_distance_bound && spinBox_radius_bound && spinBox_isovalue)
  {
    connect(spinBox_isovalue, SIGNAL(valueChanged(double)),
            this, SLOT(set_isovalue(double)));
    connect(spinBox_radius_bound, SIGNAL(valueChanged(double)),
            this, SLOT(set_radius_bound(double)));
    connect(spinBox_distance_bound, SIGNAL(valueChanged(double)),
            this, SLOT(set_distance_bound(double)));
  }
  else
    CGAL_error_msg("Cannot find spinboxes");

  connect(marching_cube_action, SIGNAL(triggered()),
          this, SLOT(display_marchin_cube()));
  connect(surface_mesher_action, SIGNAL(triggered()),
          this, SLOT(display_surface_mesher_result()));
  QAction* inverse_normals = parent->findChild<QAction*>("actionInverse_normals");
  if(inverse_normals) {
    inverse_normals->setVisible(true);
    connect(inverse_normals, SIGNAL(toggled(bool)),
            this, SLOT(set_inverse_normals(bool)));
  }
  else
    CGAL_error_msg("Cannot find action actionInverse_normals!");
  QAction* two_sides_action = parent->findChild<QAction*>("actionDisplay_front_and_back");
  if(two_sides_action) {
    two_sides_action->setVisible(true);
    connect(two_sides_action, SIGNAL(toggled(bool)),
            this, SLOT(set_two_sides(bool)));
  }
  else
    CGAL_error_msg("Cannot find action actionDisplay_front_and_back!");

  Viewer* viewer = parent->findChild<Viewer*>("viewer");
  if(viewer)
    connect(this, SIGNAL(new_bounding_box(double, double, double, double, double, double)),
            viewer, SLOT(interpolateToFitBoundingBox(double, double, double, double, double, double)));
  else
    CGAL_error_msg("Cannot find the viewer!");
}

void Volume::set_inverse_normals(const bool b) {
  m_inverse_normals = b;
  emit changed();
}

void Volume::set_two_sides(const bool b) {
  two_sides = b;
  emit changed();
}

void Volume::open(const QString& filename)
{
  if(!m_image.read(filename.toStdString().c_str(), 0.f))
    status_message(QString("Opening of file %1 failed!").arg(filename));
  else {
    status_message(QString("File %1 successfully opened.").arg(filename));
    viewer->camera()->setSceneBoundingBox(qglviewer::Vec(0, 0, 0),
                                          qglviewer::Vec(m_image.xmax(),
                                                         m_image.ymax(),
                                                         m_image.zmax()));

    viewer->showEntireScene();
    spinBox_isovalue->setRange(m_image.min_value, m_image.max_value);
    emit changed();
  }
}

void Volume::status_message(QString string)
{
  parent->statusBar()->showMessage(string, 20000);
}

void Volume::busy() const 
{
  QApplication::setOverrideCursor(QCursor(Qt::WaitCursor));
}

void Volume::not_busy() const 
{
  QApplication::restoreOverrideCursor();
}

void Volume::display_marchin_cube()
{
  MarchingCubes mc ;
  unsigned int nx = m_image.xdim();
  unsigned int ny = m_image.ydim();
  unsigned int nz = m_image.zdim();
  if(nx * ny * nz == 0)
  {
    status_message("No volume loaded.");
    return;
  }

  CGAL::Timer timer;
  busy();
  status_message("Marching cubes...");

  timer.start();
  mc.set_resolution(nx,ny,nz) ;
  mc.init_all() ;

  // set data
  for(unsigned int i=0;i<nx;i++)
    for(unsigned int j=0;j<ny;j++)
      for(unsigned int k=0;k<nz;k++)
      {
        const float& value = m_image.value(i,j,k);
        mc.set_data(value,i,j,k);
      }
  mc.run(m_image.isovalue());

  // compute scaling ratio
   const double xr = m_image.xmax() / nx;
   const double yr = m_image.ymax() / ny;
   const double zr = m_image.zmax() / nz;

  std::vector<double> facets;
  mc.get_facets(facets,xr,yr,zr);
  m_surface_mc.clear();

  timer.stop();
  not_busy();

  CGAL::Bbox_3 bbox(0,0,0,0,0,0);

  const unsigned int nbt = facets.size() / 9;
  for(unsigned int i=0;i<nbt;i++)
  {
    const Point a(facets[9*i],   facets[9*i+1], facets[9*i+2]);
    const Point b(facets[9*i+3], facets[9*i+4], facets[9*i+5]);
    const Point c(facets[9*i+6], facets[9*i+7], facets[9*i+8]);
    const Triangle_3 t(a,b,c);
    const Vector u = t[1] - t[0];
    const Vector v = t[2] - t[0];
    Vector n = CGAL::cross_product(u,v);
    n = n / std::sqrt(n*n);
    m_surface_mc.push_back(Facet(t,n));
    bbox = bbox + t.bbox();
  }

  status_message(QString("Marching cubes...done (%2 facets in %1 s)").arg(timer.time()).arg(m_surface_mc.size()));

  mc.clean_temps() ;
  mc.clean_all() ;

  m_view_mc = true;
  m_view_surface = false;
  emit changed();
  emit new_bounding_box(bbox.xmin(),
                        bbox.ymin(),
                        bbox.zmin(),
                        bbox.xmax(),
                        bbox.ymax(),
                        bbox.zmax());
}

void Volume::display_surface_mesher_result()
{
  unsigned int nx = m_image.xdim();
  unsigned int ny = m_image.ydim();
  unsigned int nz = m_image.zdim();
  if(nx * ny * nz == 0)
  {
    status_message("No volume loaded.");
    return;
  }

  m_surface.clear();
  CGAL::Timer timer;
  busy();

  status_message("Surface meshing...");

  timer.start();

  Tr tr;            // 3D-Delaunay triangulation
  C2t3 c2t3 (tr);   // 2D-complex in 3D-Delaunay triangulation

  // Carefully choosen bounding sphere: the center must be inside the
  // surface defined by 'image' and the radius must be high enough so that
  // the sphere actually bounds the whole image.

//   bool probe = m_image.probe_sink(10000);
//   if(!probe)
//   {
//     status_message("Sink probing failed");
//     not_busy();
//     return;
//   }

  Sphere bounding_sphere(m_image.center(),m_image.radius()*m_image.radius());

  // definition of the surface
  Surface_3 surface(m_image, bounding_sphere, m_relative_precision);

  std::vector<Point> seeds;
  search_for_connected_components(std::back_inserter(seeds));

  Oracle oracle;

  for(std::vector<Point>::const_iterator it = seeds.begin(), end = seeds.end();
      it != end; ++it)
  {
    CGAL::Random_points_on_sphere_3<Point> random_points_on_sphere_3(2*m_image.radius());
    Oracle::Intersect_3 intersect = oracle.intersect_3_object();
    for(int i = 0; i < 20; ++i)
    {
      const Point test = *it + (*random_points_on_sphere_3++ - CGAL::ORIGIN);
      CGAL::Object o = intersect(surface, Segment_3(*it, test));
      if (const Point* intersection = CGAL::object_cast<Point>(&o))
        tr.insert(*intersection);
      else 
      {
        std::cerr << 
          boost::format("Error. Segment (%1%, %2%) does not intersect the surface! values=(%3%, %4%)\n")
          % *it % test
          % surface(*it) % surface(test);
      }
    }
  }

  std::cerr << boost::format("Number of initial points: %1%\n") % tr.number_of_vertices();

  // defining meshing criteria
  CGAL::Surface_mesh_default_criteria_3<Tr> criteria(m_sm_angle,  // angular bound
                                                     m_sm_radius,  // radius bound
                                                     m_sm_distance); // distance bound
  std::cerr << "Surface_mesher... angle=" << m_sm_angle << ", radius= " << m_sm_radius
            << ", distance=" << m_sm_distance << "\n";
  // meshing surface
  make_surface_mesh(c2t3, surface, criteria, CGAL::Manifold_tag(), 0);
  timer.stop();
  not_busy();

  CGAL::Bbox_3 bbox(0,0,0,0,0,0);

  // get output surface
  for(C2t3::Facet_iterator
        fit = c2t3.facets_begin(), end = c2t3.facets_end();
      fit != end; ++fit)
  {
    const Tr::Cell_handle& cell = fit->first;
    const int index = fit->second;
    const Triangle_3 t = 
      Triangle_3(cell->vertex(tr.vertex_triple_index(index, 0))->point(),
                 cell->vertex(tr.vertex_triple_index(index, 1))->point(),
                 cell->vertex(tr.vertex_triple_index(index, 2))->point());
    bbox = bbox + t.bbox();
    const Vector u = t[1] - t[0];
    const Vector v = t[2] - t[0];
    Vector n = CGAL::cross_product(u,v);
    n = n / std::sqrt(n*n);
    m_surface.push_back(Facet(t,n));
  }

  const unsigned int nbt = m_surface.size();
  status_message(QString("Surface meshing...done (%1 triangles, %2 s)").arg(nbt).arg(timer.time()));

  // toggle visualization
  m_view_mc = false;
  m_view_surface = true;

  emit changed();
  emit new_bounding_box(bbox.xmin(),
                        bbox.ymin(),
                        bbox.zmin(),
                        bbox.xmax(),
                        bbox.ymax(),
                        bbox.zmax());
}

void Volume::draw()
{
  float	ambient[]  =   { 0.25f,
                         0.20725f,
                         0.20725f,
                         0.922f };
  float	diffuse[]  =   { 1.0f,
                         0.829f,
                         0.829f,
                         0.922f };

  float	specular[]  = {  0.296648f,
                         0.296648f,
                         0.296648f,
                         0.522f };

  float	emission[]  = {  0.3f,
                         0.3f,
                         0.3f,
                         1.0f };
  float shininess[] = {  11.264f };

  // apply
  ::glMaterialfv( GL_FRONT_AND_BACK, GL_AMBIENT,   ambient);
  ::glMaterialfv( GL_FRONT_AND_BACK, GL_DIFFUSE,   diffuse);
  ::glMaterialfv( GL_FRONT_AND_BACK, GL_SPECULAR,  specular);
  ::glMaterialfv( GL_FRONT_AND_BACK, GL_SHININESS, shininess);
  ::glMaterialfv( GL_FRONT_AND_BACK, GL_EMISSION,  emission);

  ::glEnable(GL_LINE_SMOOTH);

  if(two_sides)
    ::glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
  else
    ::glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE); // default

  // draw surface mesh
  if(m_view_surface)
  {
    ::glEnable(GL_LIGHTING);
    ::glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
    ::glColor3f(0.2f, 0.2f, 1.f);
    ::glEnable(GL_POLYGON_OFFSET_FILL);
    ::glPolygonOffset(3.0f,-3.0f);
    gl_draw_surface(m_inverse_normals);

    ::glDisable(GL_LIGHTING);
    ::glLineWidth(1.);
    ::glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
    ::glColor3ub(0,0,0);
    ::glDisable(GL_POLYGON_OFFSET_FILL);
    gl_draw_surface(m_inverse_normals);
  }

  // draw MC surface mesh
  if(m_view_mc)
  {
    ::glEnable(GL_LIGHTING);
    ::glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
    ::glColor3f(0.2f, 0.2f, 1.f);
    ::glEnable(GL_POLYGON_OFFSET_FILL);
    ::glPolygonOffset(3.0f,-3.0f);
    gl_draw_surface_mc(m_inverse_normals);

    ::glDisable(GL_LIGHTING);
    ::glLineWidth(1.);
    ::glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
    ::glColor3ub(0,0,0);
    ::glDisable(GL_POLYGON_OFFSET_FILL);
    gl_draw_surface_mc(m_inverse_normals);
  }

  ::glDisable(GL_LIGHTING);
  m_image.gl_draw_bbox(3.0f,0,0,0);
}

void Volume::set_isovalue(double d)
{
  m_isovalue = d;
  m_image.isovalue() = d;
  emit changed();
}

void Volume::set_radius_bound(double d)
{ 
  m_sm_radius = FT(d);
}

void Volume::set_distance_bound(double d)
{ 
  m_sm_distance = FT(d);
}

void Volume::gl_draw_surface(bool inverse_normals)
{
  gl_draw_surface(m_surface.begin(),
                  m_surface.end(),
                  inverse_normals);
}

void Volume::gl_draw_surface_mc(bool inverse_normals)
{
  gl_draw_surface(m_surface_mc.begin(),
                  m_surface_mc.end(),
                  inverse_normals);
}
