#include <CGAL/basic.h>

#include  <algorithm> // std::sort
#include <boost/shared_ptr.hpp>

#include <CGAL/Bbox_3.h>
#include <CGAL/Timer.h>

#include "volume.h"
#include "viewer.h"
#include "mainwindow.h"
#include "isovalues_list.h"

#include <QApplication>
#include <QAction>
#include <QStatusBar>
#include <QDoubleSpinBox>
#include <QMessageBox>
#include <QTreeWidgetItem>
#include <QTime>
#include <QColor>
#include <QColorDialog>

#include <GL/glu.h>

#include <CGAL/Surface_mesher/Standard_criteria.h>
#include <CGAL/Surface_mesher/Vertices_on_the_same_psc_element_criterion.h>


struct Threshold : public std::unary_function<FT, unsigned char> {
  double isovalue;

  Threshold(double isovalue) : isovalue(isovalue) {}

  result_type operator()(FT value)
  {
    if(value >=  isovalue)
      return 1;
    else
      return 0;
  }
};

class Classify_from_isovalue_list :
  public std::unary_function<FT, unsigned char> 
{
  typedef std::pair<FT, result_type> Isovalue;
  typedef std::vector<Isovalue> Isovalues;
  boost::shared_ptr<Isovalues> isovalues;

  struct Sort_isovalues : std::binary_function<Isovalue, Isovalue, bool> 
  {
    bool operator()(const Isovalue& isoval1, const Isovalue& isoval2)
    {
      return isoval1.first < isoval2.first;
    }
  };
public:
  Classify_from_isovalue_list(Isovalues_list * list)
  {
    isovalues = boost::shared_ptr<Isovalues>(new Isovalues(list->numberOfIsoValues()));
    for(int i = 0, nbs = list->numberOfIsoValues(); i < nbs; ++i )
      (*isovalues)[i] = std::make_pair(list->isovalue(i), i);
    std::sort(isovalues->begin(), isovalues->end(), Sort_isovalues());
  }

  result_type operator()(FT value)
  {
    result_type result = 0;
//     std::cerr << "isovalues: ";
    for(int i = 1, end = isovalues->size(); i <= end; ++i)
    {
//       std::cerr << (*isovalues)[i-1] << ", ";
      if(value >= (*isovalues)[i-1].first &&
         i >= result)
      {
        result = i;
      }
    }
//     if(result>1)
//       std::cerr << "result = "  << (int)result << "/" << list->numberOfIsoValues() << std::endl;
//     else
//       std::cerr << std::endl;
    if(result>0)
      return (*isovalues)[result-1].second + 1;
    else
      return 0;
  }
};

class Generate_surface_identifiers :
  std::binary_function<Classify_from_isovalue_list::result_type,
                       Classify_from_isovalue_list::result_type,
                       const QTreeWidgetItem*>
{
  Isovalues_list* list;
public:
  Generate_surface_identifiers(Isovalues_list* list) : list(list) {};

  result_type operator()(const Classify_from_isovalue_list::result_type& a,
                         const Classify_from_isovalue_list::result_type& b)
  {
    return list->item((std::min)(a, b));
  }
};

// class Classify_from_isovalue_list :
//   public std::unary_function<FT, const QTreeWidgetItem*> 
// {
//   typedef std::pair<FT, result_type> Isovalue;
//   typedef std::vector<Isovalue> Isovalues;
//   boost::shared_ptr<Isovalues> isovalues;

//   struct Sort_isovalues : std::binary_function<Isovalue, Isovalue, bool> 
//   {
//     bool operator()(const Isovalue& isoval1, const Isovalue& isoval2)
//     {
//       return isoval1.first < isoval2.first;
//     }
//   };
// public:
//   Classify_from_isovalue_list(Isovalues_list * list)
//   {
//     isovalues = boost::shared_ptr<Isovalues>(new Isovalues(list->numberOfIsoValues()));
//     for(int i = 0, nbs = list->numberOfIsoValues(); i < nbs; ++i )
//       (*isovalues)[i] = std::make_pair(list->isovalue(i), list->item(i));
//     std::sort(isovalues->begin(), isovalues->end(), Sort_isovalues());
//   }

//   result_type operator()(FT value)
//   {
//     int result = 0;
// //     std::cerr << "isovalues: ";
//     for(int i = 1, end = isovalues->size(); i <= end; ++i)
//     {
// //       std::cerr << (*isovalues)[i-1] << ", ";
//       if(value >= (*isovalues)[i-1].first &&
//          i >= result)
//       {
//         result = i;
//       }
//     }
//     if(result>1)
//       std::cerr << boost::format("result = %1%/%2%\n") % result % isovalues->size();
//     if(result>0)
//       return (*isovalues)[result-1].second;
//     else
//       return 0;
//   }
// };
Volume::Volume(MainWindow* mw) : 
  Surface(mw),
  m_sm_angle(30),
  m_sm_radius(0),
  m_sm_distance(0),
  m_relative_precision(0.0001),
  m_view_surface(false),
  m_view_mc(false),
  m_triangulation_color(QColor(Qt::green)),
  mw(mw),
  m_inverse_normals(false),
  two_sides(false),
  list_draw_marching_cube(0),
  lists_draw_surface(),
  lists_draw_surface_is_valid(false),
  lists_draw_surface_mc(),
  list_draw_marching_cube_is_valid(false)
{
  spinBox_radius_bound = mw->findChild<QDoubleSpinBox*>("spinBox_radius_bound");
  spinBox_distance_bound = mw->findChild<QDoubleSpinBox*>("spinBox_distance_bound");
  Q_ASSERT_X(spinBox_radius_bound && spinBox_distance_bound,
             "Volume::Volume()", "Cannot find spinboxes!");

  isovalues_list = mw->isovalues;

  connect(spinBox_radius_bound, SIGNAL(valueChanged(double)),
          this, SLOT(set_radius_bound(double)));
  connect(spinBox_distance_bound, SIGNAL(valueChanged(double)),
          this, SLOT(set_distance_bound(double)));

  connect(mw->actionMarching_cubes, SIGNAL(triggered()),
          this, SLOT(display_marchin_cube()));
  connect(mw->actionSurface_mesher, SIGNAL(triggered()),
          this, SLOT(display_surface_mesher_result()));

  mw->actionInverse_normals->setVisible(true);
  connect(mw->actionInverse_normals, SIGNAL(toggled(bool)),
          this, SLOT(set_inverse_normals(bool)));
  m_inverse_normals = mw->actionInverse_normals->isChecked();

  mw->actionDisplay_front_and_back->setVisible(true);
  connect(mw->actionDisplay_front_and_back, SIGNAL(toggled(bool)),
          this, SLOT(set_two_sides(bool)));
  two_sides = mw->actionDisplay_front_and_back->isChecked();

  mw->actionDraw_triangles_edges->setVisible(true);
  connect(mw->actionDraw_triangles_edges, SIGNAL(toggled(bool)),
          this, SLOT(set_draw_triangles_edges(bool)));
  draw_triangles_edges = mw->actionDraw_triangles_edges->isChecked();

  mw->actionUse_Gouraud_shading->setVisible(true);
  connect(mw->actionUse_Gouraud_shading, SIGNAL(toggled(bool)),
          this, SLOT(set_use_gouraud(bool)));
  use_gouraud = mw->actionUse_Gouraud_shading->isChecked();

  mw->actionShow_triangulation->setVisible(true);
  connect(mw->actionShow_triangulation, SIGNAL(toggled(bool)),
          this, SLOT(set_draw_triangulation(bool)));
  m_draw_triangulation = mw->actionShow_triangulation->isChecked();

  mw->actionTriangulation_edges_color->setVisible(true);
  connect(mw->actionTriangulation_edges_color, SIGNAL(triggered()),
          this, SLOT(set_triangulation_edges_color()));

  connect(this, SIGNAL(new_bounding_box(double, double, double, double, double, double)),
          mw->viewer, SLOT(interpolateToFitBoundingBox(double, double, double, double, double, double)));

  connect(isovalues_list, SIGNAL(isovalues_changed()),
          this, SLOT(changed_parameters()));
  connect(isovalues_list, SIGNAL(changed()),
          mw->viewer, SLOT(updateGL()));
}

void Volume::set_inverse_normals(const bool b) {
  m_inverse_normals = b;

  list_draw_marching_cube = 0; // Invalidate the display list for the
                               // marching cube. See gl_draw_marchingcube()
                               // for an explanation.

  emit changed();
}

void Volume::set_two_sides(const bool b) {
  two_sides = b;
  emit changed();
}

void Volume::set_draw_triangles_edges(const bool b) {
  draw_triangles_edges = b;
  emit changed();
}

void Volume::set_draw_triangulation(const bool b) {
  m_draw_triangulation = b;
  emit changed();
}

void Volume::set_triangulation_edges_color() {
  const QColor color = QColorDialog::getColor(m_triangulation_color, mw);
  if (color.isValid()) {
    m_triangulation_color = color;
    emit changed();
  }
}

void Volume::set_use_gouraud(const bool b) {
  use_gouraud = b;
  emit changed();
}

void Volume::open(const QString& filename)
{
  fileinfo.setFile(filename);
  if(!fileinfo.isReadable())
  {
    QMessageBox::warning(mw, mw->windowTitle(),
                         QString(tr("Cannot read file <tt>%1</tt>!")).arg(filename));
    status_message(QString("Opening of file %1 failed!").arg(filename));
  }
  else
  {
    if(!m_image.read(filename.toStdString().c_str()))
    {
      QMessageBox::warning(mw, mw->windowTitle(),
                           QString(tr("Error with file <tt>%1</tt>:\nunknown file format!")).arg(filename));
      status_message(QString("Opening of file %1 failed!").arg(filename));
    }
    else
    {
      status_message(QString("File %1 successfully opened.").arg(filename));
      mw->viewer->camera()->setSceneBoundingBox(qglviewer::Vec(0, 0, 0),
                                                qglviewer::Vec(m_image.xmax(),
                                                               m_image.ymax(),
                                                               m_image.zmax()));

      mw->viewer->showEntireScene();
      isovalues_list->load_values(fileinfo.absoluteFilePath());
      changed_parameters();
      emit changed();
    }
  }
}

void Volume::status_message(QString string)
{
  std::cerr << qPrintable(string) << std::endl;
  mw->statusBar()->showMessage(string, 20000);
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
  if(m_surface_mc.empty())
  {
    QTime total_time;
    total_time.start();

    isovalues_list->save_values(fileinfo.absoluteFilePath());

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
    m_surface_mc.clear();

    if(mc.ntrigs()!=0)
      mc.clean_all();
    mc.set_resolution(nx,ny,nz);
    mc.init_all();
    mc.set_ext_data(static_cast<unsigned char*>(m_image.image()->data));

    const double xr = m_image.xmax() / nx;
    const double yr = m_image.ymax() / ny;
    const double zr = m_image.zmax() / nz;

    nbs_of_mc_triangles.resize(isovalues_list->numberOfIsoValues());

    for(int isovalue_id = 0; 
        isovalue_id < isovalues_list->numberOfIsoValues();
        ++isovalue_id)
    {
      status_message(QString("Marching cubes, isovalue #%1...").arg(isovalue_id));

      // set data
//       for(unsigned int i=0;i<nx;i++)
//         for(unsigned int j=0;j<ny;j++)
//           for(unsigned int k=0;k<nz;k++)
//           {
//             const float& value = m_image.value(i,j,k);
//             mc.set_data(value,i,j,k);
//           }
      // compute scaling ratio
      if(isovalue_id > 0)
        mc.init_temps();
      mc.run(isovalues_list->isovalue(isovalue_id), xr, yr, zr);
      mc.clean_temps();

      std::vector<double> facets;
      mc.get_facets(facets);

      timer.stop();
      const unsigned int begin = isovalue_id == 0 ? 0 : nbs_of_mc_triangles[isovalue_id-1];
      const unsigned int nbt = facets.size() / 9;
      for(unsigned int i=begin;i<nbt;i++)
      {
        const Point a(facets[9*i],   facets[9*i+1], facets[9*i+2]);
        const Point b(facets[9*i+3], facets[9*i+4], facets[9*i+5]);
        const Point c(facets[9*i+6], facets[9*i+7], facets[9*i+8]);
        const Triangle_3 t(a,b,c);
        const Vector u = t[1] - t[0];
        const Vector v = t[2] - t[0];
        Vector n = CGAL::cross_product(u,v);
        n = n / std::sqrt(n*n);
        m_surface_mc.push_back(Facet(t,n,isovalues_list->item(isovalue_id)));
      }
      nbs_of_mc_triangles[isovalue_id]=m_surface_mc.size();
      timer.start();
    }
    timer.stop();
    not_busy();

    status_message(QString("Marching cubes...done. %2 facets in %1s (CPU time), total time is %3s.")
                   .arg(timer.time())
                   .arg(m_surface_mc.size())
                   .arg(total_time.elapsed()/1000.));

    // invalidate the display list
    lists_draw_surface_mc_is_valid = false;
    list_draw_marching_cube_is_valid = false;
  }
  CGAL::Bbox_3 bbox(0,0,0,0,0,0);
  for(std::vector<Facet>::const_iterator
        it = m_surface_mc.begin(), end = m_surface_mc.end();
      it != end; ++it)
  {
    bbox = bbox + it->first.bbox();
  }

  m_view_mc = true;
  m_view_surface = false;
  emit changed();
  if(!m_surface_mc.empty())
  {
    emit new_bounding_box(bbox.xmin(),
                          bbox.ymin(),
                          bbox.zmin(),
                          bbox.xmax(),
                          bbox.ymax(),
                          bbox.zmax());
  }
}

void Volume::display_surface_mesher_result()
{
  if(m_surface.empty() || // Either the surface is not computed.
     m_view_surface) // Or it is computed and displayed, and one want
                     // to recompute it.
  {
    QTime total_time;
    total_time.start();

    isovalues_list->save_values(fileinfo.absoluteFilePath());

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

    del.clear();
    C2t3 c2t3(del);   // 2D-complex in 3D-Delaunay triangulation
    Sphere bounding_sphere(m_image.center(),m_image.radius()*m_image.radius());

    // definition of the surface
    Surface_3 surface(m_image, bounding_sphere, m_relative_precision);
//     Threshold threshold(m_image.isovalue());
    Classify_from_isovalue_list classify(isovalues_list);
    Generate_surface_identifiers generate_ids(isovalues_list);

    std::vector<Point> seeds;
    search_for_connected_components(std::back_inserter(seeds), classify);

    // surface mesh traits class
    typedef CGAL::Surface_mesher::Implicit_surface_oracle_3<Kernel,
      Surface_3, 
      Classify_from_isovalue_list,
      Generate_surface_identifiers> Oracle;
    Oracle oracle(classify, generate_ids);

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
          del.insert(*intersection);
        else 
        {
          std::cerr << 
            boost::format("Error. Segment (%1%, %2%) does not intersect the surface! values=(%3%, %4%)\n")
            % *it % test
            % surface(*it) % surface(test);
        }
      }
    }

    std::cerr << boost::format("Number of initial points: %1%\n") % del.number_of_vertices();

    // defining meshing criteria
    typedef CGAL::Surface_mesher::Refine_criterion<Tr> Criterion;
    CGAL::Surface_mesher::Curvature_size_criterion<Tr>
      curvature_size_criterion (m_sm_distance);
    CGAL::Surface_mesher::Uniform_size_criterion<Tr>
      uniform_size_criterion (m_sm_radius);
    CGAL::Surface_mesher::Aspect_ratio_criterion<Tr>
      aspect_ratio_criterion (m_sm_angle);
    CGAL::Surface_mesher::Vertices_on_the_same_psc_element_criterion<Tr, Surface_3>
      vertices_on_the_same_psc_element_criterion(surface);
    
    std::vector<Criterion*> criterion_vector;
    criterion_vector.push_back(&aspect_ratio_criterion);
    criterion_vector.push_back(&uniform_size_criterion);
    criterion_vector.push_back(&curvature_size_criterion);
    criterion_vector.push_back(&vertices_on_the_same_psc_element_criterion);

    CGAL::Surface_mesher::Standard_criteria<Criterion> criteria(criterion_vector);
    std::cerr << "Surface_mesher... angle=" << m_sm_angle << ", radius= " << m_sm_radius
              << ", distance=" << m_sm_distance << "\n";

    // meshing surface
    make_surface_mesh(c2t3, surface, oracle, criteria, CGAL::Manifold_tag(), 0);
    timer.stop();
    not_busy();

    // get output surface
    for(C2t3::Facet_iterator
          fit = c2t3.facets_begin(), end = c2t3.facets_end();
        fit != end; ++fit)
    {
      const Tr::Cell_handle& cell = fit->first;
      const int index = fit->second;
      const Triangle_3 t = 
        Triangle_3(cell->vertex(del.vertex_triple_index(index, 0))->point(),
                   cell->vertex(del.vertex_triple_index(index, 1))->point(),
                   cell->vertex(del.vertex_triple_index(index, 2))->point());
      const Vector u = t[1] - t[0];
      const Vector v = t[2] - t[0];
      Vector n = CGAL::cross_product(u,v);
      n = n / std::sqrt(n*n);
      m_surface.push_back(Facet(t,n,cell->vertex(del.vertex_triple_index(index, 0))->point().element_index()));
    }

    const unsigned int nbt = m_surface.size();
    status_message(QString("Surface meshing...done. %1 facets in %2s (CPU time), total time is %3s.)")
                   .arg(nbt)
                   .arg(timer.time())
                   .arg(total_time.elapsed()/1000.));

    // invalidate the display list
    lists_draw_surface_is_valid = false;
  }

  CGAL::Bbox_3 bbox(0,0,0,0,0,0);
  for(std::vector<Facet>::const_iterator
        it = m_surface.begin(), end = m_surface.end();
      it != end; ++it)
  {
    bbox = bbox + it->first.bbox();
  }

  // toggle visualization
  m_view_mc = false;
  m_view_surface = true;
  emit changed();
  if(!m_surface.empty())
  {
    emit new_bounding_box(bbox.xmin(),
                          bbox.ymin(),
                          bbox.zmin(),
                          bbox.xmax(),
                          bbox.ymax(),
                          bbox.zmax());
  }
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
    gl_draw_surface();

    if(draw_triangles_edges)
    {
      ::glDisable(GL_LIGHTING);
      ::glLineWidth(1.);
      ::glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
      ::glColor3ub(0,0,0);
      ::glDisable(GL_POLYGON_OFFSET_FILL);
      gl_draw_surface();
    }
  }

  // draw MC surface mesh
  if(m_view_mc)
  {
    ::glEnable(GL_LIGHTING);
    ::glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
    ::glColor3f(0.2f, 0.2f, 1.f);
    ::glEnable(GL_POLYGON_OFFSET_FILL);
    ::glPolygonOffset(3.0f,-3.0f);
    gl_draw_surface_mc();

    if(draw_triangles_edges)
    {
      ::glDisable(GL_LIGHTING);
      ::glLineWidth(1.);
      ::glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
      ::glColor3ub(0,0,0);
      ::glDisable(GL_POLYGON_OFFSET_FILL);
      gl_draw_surface_mc();
    }
  }

  ::glDisable(GL_LIGHTING);
  m_image.gl_draw_bbox(3.0f,0,0,0);


  if(!m_view_mc && m_draw_triangulation)
  {
    // draw the triangualtion
    mw->viewer->qglColor(m_triangulation_color);
    ::glLineWidth(1.0);
    ::glBegin(GL_LINES);
    for(Tr::Finite_edges_iterator 
          eit = del.finite_edges_begin(), 
          end = del.finite_edges_end();
        eit != end; ++eit) 
    {
      const Point p1 = eit->first->vertex(eit->second)->point();
      const Point p2 = eit->first->vertex(eit->third)->point();
      ::glVertex3d(p1.x(),p1.y(),p1.z());
      ::glVertex3d(p2.x(),p2.y(),p2.z());
    }
    ::glEnd();
  }
}

void Volume::set_radius_bound(double d)
{ 
  m_sm_radius = FT(d);
  changed_parameters();
}

void Volume::set_distance_bound(double d)
{ 
  m_sm_distance = FT(d);
  changed_parameters();
}

void Volume::gl_draw_surface_mc()
{
  if(use_gouraud)
  {
    gl_draw_marchingcube();
    return;
  }
  
  if(lists_draw_surface_mc_is_valid)
  {
    for(int i = 0, nbs = isovalues_list->numberOfIsoValues(); i < nbs; ++i )
    {
      if(isovalues_list->enabled(i))
      {
        mw->viewer->qglColor(isovalues_list->color(i));
        ::glCallList(lists_draw_surface_mc[i]);
      }
    }
  }
  else
  {
    lists_draw_surface_mc.resize(isovalues_list->numberOfIsoValues(), 0);
    for(int i = 0, nbs = isovalues_list->numberOfIsoValues(); i < nbs; ++i )
    {
      if(!lists_draw_surface_mc[i])
      {
        lists_draw_surface_mc[i] = ::glGenLists(1);
      }

      std::cerr << boost::format("(Re-)Generating list #%1% for marching cube surface #%2%"
                                 " in gl_draw_surface(), ()\n")
        % lists_draw_surface_mc[i]
        % i;

      mw->viewer->qglColor(isovalues_list->color(i));

      if(lists_draw_surface_mc[i])             // If
        ::glNewList(lists_draw_surface_mc[i],  // lists_draw_surface[i]==0 then something
                    isovalues_list->enabled(i) // got wrong in the list generation.
                    ? GL_COMPILE_AND_EXECUTE
                    : GL_COMPILE);


      gl_draw_surface(m_surface_mc.begin(),
                      m_surface_mc.end(),
                      isovalues_list->item(i));
        
      if(lists_draw_surface_mc[i]) // If lists_draw_surface[i]==0 then
      {                            // something got wrong in the list
        ::glEndList();             // generation.
      }
    }
    lists_draw_surface_mc_is_valid = (::glGetError() == GL_NO_ERROR);
  }
}

void Volume::gl_draw_surface()
{
  if(lists_draw_surface_is_valid)
  {
    for(int i = 0, nbs = isovalues_list->numberOfIsoValues(); i < nbs; ++i )
    {
      if(isovalues_list->enabled(i))
      {
        mw->viewer->qglColor(isovalues_list->color(i));
        ::glCallList(lists_draw_surface[i]);
      }
    }
  }
  else
  {
    lists_draw_surface.resize(isovalues_list->numberOfIsoValues(), 0);
    for(int i = 0, nbs = isovalues_list->numberOfIsoValues(); i < nbs; ++i )
    {
      if(!lists_draw_surface[i])
      {
        lists_draw_surface[i] = ::glGenLists(1);
      }

      std::cerr << boost::format("(Re-)Generating list #%1% for surface #%2%"
                                 " in gl_draw_surface(), ()\n")
        % lists_draw_surface[i]
        % i;
        
      mw->viewer->qglColor(isovalues_list->color(i));

      if(lists_draw_surface[i])                 // If
        ::glNewList(lists_draw_surface[i],      // lists_draw_surface[i]==0
                    isovalues_list->enabled(i)  // then something got wrong
                    ? GL_COMPILE_AND_EXECUTE    // in the list generation.
                    : GL_COMPILE);

      gl_draw_surface(m_surface.begin(),
                      m_surface.end(),
                      isovalues_list->item(i));
        
      if(lists_draw_surface[i]) // If lists_draw_surface[i]==0 then
      {                         // something got wrong in the list
        ::glEndList();          // generation.
      }
    }
    lists_draw_surface_is_valid = (::glGetError() == GL_NO_ERROR);
  }
}

template <typename Iterator>
void Volume::gl_draw_surface(Iterator begin, Iterator end, const QTreeWidgetItem* i)
{
  ::glBegin(GL_TRIANGLES);
  unsigned int counter = 0;
  for(Iterator it = begin; it != end; ++it)
  {
    const Facet& f = *it;

    if(f.third != i) continue;

    const Vector& n = f.second;

    if(m_inverse_normals)
      ::glNormal3d(-n.x(),-n.y(),-n.z());
    else
      ::glNormal3d(n.x(),n.y(),n.z());

    const Triangle_3& t = f.first;
    const Point& a = t[0];
    const Point& b = t[1];
    const Point& c = t[2];

    ::glVertex3d(a.x(),a.y(),a.z());
    ::glVertex3d(b.x(),b.y(),b.z());
    ::glVertex3d(c.x(),c.y(),c.z());
    ++counter;
  }
  ::glEnd();
  std::cerr << boost::format("number of facets: %1%\n")
    % counter;
}

void Volume::changed_parameters()
{
  m_surface.clear();
  m_surface_mc.clear();
}

void Volume::gl_draw_one_marching_cube_vertex(int i)
{
  if(!m_inverse_normals)
    glArrayElement(i);
  else
  {
    const Vertex* const vertex = mc.vert(i);
    ::glNormal3d(-vertex->nx, -vertex->ny, -vertex->nz);
    ::glVertex3d(vertex->x, vertex->y, vertex->z);
  }
}

void Volume::gl_draw_marchingcube()
{
  if(list_draw_marching_cube_is_valid)
    ::glCallList(list_draw_marching_cube);
  else
  {
    if(!list_draw_marching_cube)
      list_draw_marching_cube = ::glGenLists(1);
    std::cerr << boost::format("(Re-)Generating list #%1% for"
                               " gl_draw_marchingcube()\n")
      % list_draw_marching_cube;

    if(list_draw_marching_cube)          // If list_draw_marching_cube==0 then
    ::glNewList(list_draw_marching_cube, // something got wrong in the list
                GL_COMPILE_AND_EXECUTE); // generation.

    ::glVertexPointer(3, GL_DOUBLE, sizeof(Vertex), mc.vertices());
    ::glNormalPointer(GL_DOUBLE, sizeof(Vertex), &(mc.vertices()->nx));
    ::glEnableClientState(GL_VERTEX_ARRAY);

    // because of that conditionnal, the display list has to be
    // reconstructed each time m_inverse_normals is toggled.
    if(!m_inverse_normals)
      ::glEnableClientState(GL_NORMAL_ARRAY);
    const int size = mc.ntrigs();

    for(int i = 0, nbs = isovalues_list->numberOfIsoValues(); i < nbs; ++i)
    {
      const int begin = i == 0 ? 0 : nbs_of_mc_triangles[i-1];
      const int end = nbs_of_mc_triangles[i];
      mw->viewer->qglColor(isovalues_list->color(i));
      ::glBegin(GL_TRIANGLES);
      for(int i = begin; i < end; ++i)
      {
        const MC_Triangle* const trig = mc.trig(i);
        gl_draw_one_marching_cube_vertex(trig->v1);
        gl_draw_one_marching_cube_vertex(trig->v2);
        gl_draw_one_marching_cube_vertex(trig->v3);
      }
      ::glEnd();
    }
    if(list_draw_marching_cube > 0) // If list_draw_marching_cube==0 then
    {                               // something got wrong in the list
      ::glEndList();                // generation.
      list_draw_marching_cube_is_valid = (::glGetError() == GL_NO_ERROR);
    }
    if(!list_draw_marching_cube_is_valid)
      std::cerr << boost::format("OpenGL error: %1%\n") 
        % ::gluErrorString(::glGetError());
  }
}

#include "volume.moc"
