#include <CGAL/Three/Viewer_interface.h>
#include <CGAL/Three/Triangle_container.h>
#include <CGAL/Three/Edge_container.h>
#include <CGAL/Three/Point_container.h>
#include <CGAL/Three/Three.h>
#include <CGAL/Buffer_for_vao.h>

#include "Color_map.h"
#include "Scene_lcc_item.h"

//todo : create a struct for facets containing useful infos for drawing and their volume, and fill it during computeElements().
using namespace CGAL::Three;
typedef Triangle_container Tri;
typedef Edge_container Ec;
typedef Point_container Pc;
typedef Viewer_interface Vi;

typedef Scene_lcc_item::LCC::Dart_const_handle Dart_const_handle;
typedef Scene_lcc_item::LCC::Dart_handle Dart_handle;
typedef Scene_lcc_item::LCC::Point Point;

struct Facet{
  Facet():normal(Scene_lcc_item::LCC::Vector(0,0,0)){}
  Dart_const_handle f_handle;
  std::vector<Point> points;
  Scene_lcc_item::LCC::Vector normal;
  std::size_t volume_id;
  std::size_t size() { return points.size(); }
};


struct lcc_priv{

  Scene_lcc_item::LCC lcc;
  std::vector<float> faces;
  std::vector<float> lines;
  std::vector<float> vertices;
  std::vector<float> colors;

  std::vector<Facet> facets;
  std::size_t nb_volumes;
  bool is_mono_color;

  std::size_t nb_lines, nb_vertices, nb_faces;

  lcc_priv(const Scene_lcc_item::LCC& lcc)
    :lcc(lcc), is_mono_color(true){}

  bool compute_face(Dart_const_handle dh, Facet& f)
  {
    f.f_handle = dh;
    const CGAL::qglviewer::Vec offset = static_cast<CGAL::Three::Viewer_interface*>(CGAL::QGLViewer::QGLViewerPool().first())->offset();
    // We fill only closed faces.
    Dart_const_handle cur=dh;
    Dart_const_handle min=dh;
    do
    {
      if (!lcc.is_next_exist(cur)) return false; // open face=>not filled
      if (cur<min) min=cur;
      cur=lcc.next(cur);
    }
    while(cur!=dh);

    cur=dh;
    do
    {
      f.points.push_back(lcc.point(cur));
      cur=lcc.next(cur);
    }
    while(cur!=dh);
    for (std::size_t i = 0; i < f.size() ; ++ i){
      const Point& pa = f.points[i];
      const Point& pb = f.points[(i+1)%f.size()];
      double x = f.normal.x() + (pa.y()-pb.y())*(pa.z()+pb.z());
      double y = f.normal.y() + (pa.z()-pb.z())*(pa.x()+pb.x());
      double z = f.normal.z() + (pa.x()-pb.x())*(pa.y()+pb.y());
      f.normal = Scene_lcc_item::LCC::Vector(x,y,z);
    }

    if (f.size()==3)
    {
      for(auto pt = f.points.begin();
          pt != f.points.end();
          ++pt)
      {
        faces.push_back(pt->x() + offset.x);
        faces.push_back(pt->y() + offset.y);
        faces.push_back(pt->z() + offset.z);
      }
    }
    else if (CGAL::Buffer_for_vao<float, std::size_t>::is_facet_convex(f.points,f.normal))
    {
      if (f.size()==4)
      {
        Point p = f.points[0];
        faces.push_back(p.x() + offset.x);
        faces.push_back(p.y() + offset.y);
        faces.push_back(p.z() + offset.z);
        p = f.points[1];
        faces.push_back(p.x() + offset.x);
        faces.push_back(p.y() + offset.y);
        faces.push_back(p.z() + offset.z);
        p = f.points[2];
        faces.push_back(p.x() + offset.x);
        faces.push_back(p.y() + offset.y);
        faces.push_back(p.z() + offset.z);

        p = f.points[0];
        faces.push_back(p.x() + offset.x);
        faces.push_back(p.y() + offset.y);
        faces.push_back(p.z() + offset.z);
        p = f.points[2];
        faces.push_back(p.x() + offset.x);
        faces.push_back(p.y() + offset.y);
        faces.push_back(p.z() + offset.z);
        p = f.points[3];
        faces.push_back(p.x() + offset.x);
        faces.push_back(p.y() + offset.y);
        faces.push_back(p.z() + offset.z);
      }
      else
      {
        for(std::size_t i=1; i<f.size()-1; ++i)
        {
          Point& p0 = f.points[0];
          Point& p1 = f.points[i];
          Point& p2 = f.points[i+1];

          // (1) add points
          faces.push_back(p0.x() + offset.x);
          faces.push_back(p0.y() + offset.y);
          faces.push_back(p0.z() + offset.z);
          faces.push_back(p1.x() + offset.x);
          faces.push_back(p1.y() + offset.y);
          faces.push_back(p1.z() + offset.z);
          faces.push_back(p2.x() + offset.x);
          faces.push_back(p2.y() + offset.y);
          faces.push_back(p2.z() + offset.z);
        }
      } // Convex face with > 4 vertices
    }
    else
    {
      struct Vertex_info
      {
        Scene_lcc_item::LCC::Vector v;
        std::size_t index;
      };

      struct Face_info
      {
        bool exist_edge[3];
        bool is_external;
        bool is_process;
      };

      typedef CGAL::Triangulation_2_projection_traits_3<CGAL::Exact_predicates_inexact_constructions_kernel> P_traits;
      typedef CGAL::Triangulation_vertex_base_with_info_2<Vertex_info, P_traits> Vb;
      typedef CGAL::Triangulation_face_base_with_info_2<Face_info, P_traits>     Fb1;
      typedef CGAL::Constrained_triangulation_face_base_2<P_traits, Fb1>         Fb;
      typedef CGAL::Triangulation_data_structure_2<Vb,Fb>                        TDS;
      typedef CGAL::Exact_predicates_tag                                         Itag;
      typedef CGAL::Constrained_Delaunay_triangulation_2<P_traits, TDS, Itag>    CDT;

      P_traits cdt_traits(f.normal);
      CDT cdt(cdt_traits);
        // (1) We insert all the edges as contraint in the CDT.
        typename CDT::Vertex_handle previous=NULL, first=NULL;
        for (unsigned int i=0; i<f.size(); ++i)
        {
          typename CDT::Vertex_handle vh = cdt.insert(f.points[i]);
          if(first==NULL)
          { first=vh; }
          vh->info().v=f.normal;

          if(previous!=NULL && previous!=vh)
          { cdt.insert_constraint(previous, vh); }
          previous=vh;
        }

        if (previous!=NULL && previous!=first)
        { cdt.insert_constraint(previous, first); }

        // (2) We mark all external triangles
        // (2.1) We initialize is_external and is_process values
        for(typename CDT::All_faces_iterator fit = cdt.all_faces_begin(),
              fitend = cdt.all_faces_end(); fit!=fitend; ++fit)
        {
          fit->info().is_external = true;
          fit->info().is_process = false;
        }
        // (2.2) We check if the facet is external or internal
        std::queue<typename CDT::Face_handle> face_queue;
        typename CDT::Face_handle face_internal = NULL;
        if (cdt.infinite_vertex()->face()!=NULL)
        { face_queue.push(cdt.infinite_vertex()->face()); }
        while(!face_queue.empty())
        {
          typename CDT::Face_handle fh = face_queue.front();
          face_queue.pop();
          if(!fh->info().is_process)
          {
            fh->info().is_process = true;
            for(int i=0; i<3; ++i)
            {
              if(!cdt.is_constrained(std::make_pair(fh, i)))
              {
                if (fh->neighbor(i)!=NULL)
                { face_queue.push(fh->neighbor(i)); }
              }
              else if (face_internal==NULL)
              {
                face_internal = fh->neighbor(i);
              }
            }
          }
        }

        if ( face_internal!=NULL )
        { face_queue.push(face_internal); }

        while(!face_queue.empty())
        {
          typename CDT::Face_handle fh = face_queue.front();
          face_queue.pop();
          if(!fh->info().is_process)
          {
            fh->info().is_process = true;
            fh->info().is_external = false;
            for(unsigned int i=0; i<3; ++i)
            {
              if(!cdt.is_constrained(std::make_pair(fh, i)))
              {
                if (fh->neighbor(i)!=NULL)
                { face_queue.push(fh->neighbor(i)); }
              }
            }
          }
        }

        // (3) Now we iterates on the internal faces to add the vertices
        //     and the normals to the appropriate vectors
        for(typename CDT::Finite_faces_iterator ffit=cdt.finite_faces_begin(),
              ffitend = cdt.finite_faces_end(); ffit!=ffitend; ++ffit)
        {
          if(!ffit->info().is_external)
          {
            for(unsigned int i=0; i<3; ++i)
            {
              Point p = ffit->vertex(i)->point();
              faces.push_back(p.x() + offset.x);
              faces.push_back(p.y() + offset.y);
              faces.push_back(p.z() + offset.z);
            }
          }
        }
    }
    return true;
  }
};

Scene_lcc_item::Scene_lcc_item(const LCC& lcc)
  :d(new lcc_priv(lcc))
{
  d->nb_faces = 0;
  d->nb_lines = 0;
  d->nb_vertices = 0;
  d->nb_volumes = 0;
  setTriangleContainer(0,
                       new Tri(Three::mainViewer()->isOpenGL_4_3() ? Vi::PROGRAM_FLAT
                                                                   : Vi::PROGRAM_OLD_FLAT, false));

  setEdgeContainer(0,
                   new Ec(Three::mainViewer()->isOpenGL_4_3() ? Vi::PROGRAM_SOLID_WIREFRAME
                                                              : Vi::PROGRAM_NO_SELECTION
                                                                , false));
  setPointContainer(0,
                    new Pc(Vi::PROGRAM_NO_SELECTION, false));
}

Scene_lcc_item::~Scene_lcc_item()
{
  delete d;
}

Scene_lcc_item* Scene_lcc_item::clone() const
{
  Scene_lcc_item* item = new Scene_lcc_item(d->lcc);
  return item;
}

bool Scene_lcc_item::supportsRenderingMode(RenderingMode m) const
{
  return m==FlatPlusEdges;
}

QString Scene_lcc_item::toolTip() const
{
  return QString();
}

void Scene_lcc_item::compute_bbox() const
{
  Scene_item::Bbox bb;
  for (LCC::Dart_range::const_iterator it=d->lcc.darts().begin(),
       itend=d->lcc.darts().end(); it!=itend; ++it )
  {
    bb+=d->lcc.point(it).bbox();
  }
  this->setBbox(bb);
}


void Scene_lcc_item::draw(CGAL::Three::Viewer_interface* viewer) const
{
  if(!isInit(viewer))
    initGL(viewer);
  if ( getBuffersFilled() &&
       ! getBuffersInit(viewer))
  {
    initializeBuffers(viewer);
    setBuffersInit(viewer, true);
  }
  if(!getBuffersFilled())
  {
    computeElements();
    initializeBuffers(viewer);
  }
  if(d->is_mono_color)
    getTriangleContainer(0)->setColor(this->color());
  getTriangleContainer(0)->draw(viewer, d->is_mono_color);
}

void Scene_lcc_item::drawEdges(CGAL::Three::Viewer_interface* viewer) const
{
  if(!isInit(viewer))
    initGL(viewer);
  if ( getBuffersFilled() &&
       ! getBuffersInit(viewer))
  {
    initializeBuffers(viewer);
    setBuffersInit(viewer, true);
  }
  if(!getBuffersFilled())
  {
    computeElements();
    initializeBuffers(viewer);
  }
  GLfloat offset_factor;
  GLfloat offset_units;
  viewer->glGetFloatv(GL_POLYGON_OFFSET_FACTOR, &offset_factor);
  viewer->glGetFloatv(GL_POLYGON_OFFSET_UNITS, &offset_units);
  viewer->glEnable(GL_POLYGON_OFFSET_LINE);
  viewer->glPolygonOffset(0.3f, 0.3f);
  if(viewer->isOpenGL_4_3())
  {
    QVector2D vp(viewer->width(), viewer->height());

    getEdgeContainer(0)->setViewport(vp);
    getEdgeContainer(0)->setWidth(2);
  }
  getEdgeContainer(0)->setColor(QColor(Qt::black));
  getEdgeContainer(0)->draw(viewer, true);
  drawPoints(viewer);
  viewer->glDisable(GL_POLYGON_OFFSET_LINE);
  viewer->glPolygonOffset(offset_factor, offset_units);
}

void Scene_lcc_item::drawPoints(CGAL::Three::Viewer_interface* viewer) const
{
  if(!visible())
    return;
  if(!isInit(viewer))
    initGL(viewer);
  if ( getBuffersFilled() &&
       ! getBuffersInit(viewer))
  {
    initializeBuffers(viewer);
    setBuffersInit(viewer, true);
  }
  GLfloat point_size;
  viewer->glGetFloatv(GL_POINT_SIZE, &point_size);
  viewer->setGlPointSize(GLfloat(5));
  if(!getBuffersFilled())
  {
    computeElements();
    initializeBuffers(viewer);
  }

  getPointContainer(0)->setColor(QColor(Qt::black));
  getPointContainer(0)->draw(viewer, true);
  viewer->setGlPointSize(point_size);
}

void Scene_lcc_item::computeElements() const
{
  CGAL::Three::Three::CursorScopeGuard guard{QCursor(Qt::WaitCursor)};
  d->facets.clear();
  const CGAL::qglviewer::Vec offset = CGAL::Three::Three::mainViewer()->offset();
  typename LCC::size_type markvolumes  = d->lcc.get_new_mark();
  typename LCC::size_type markfaces    = d->lcc.get_new_mark();
  typename LCC::size_type markedges    = d->lcc.get_new_mark();
  typename LCC::size_type markvertices = d->lcc.get_new_mark();
  std::size_t volume_id = 0;
  for (typename LCC::Dart_range::const_iterator it=d->lcc.darts().begin(),
       itend=d->lcc.darts().end(); it!=itend; ++it )
  {
    if (!d->lcc.is_marked(it, markvolumes))
    {
      for (typename LCC::template Dart_of_cell_basic_range<3>::
           const_iterator itv=d->lcc.template darts_of_cell_basic<3>(it, markvolumes).begin(),
           itvend=d->lcc.template darts_of_cell_basic<3>(it, markvolumes).end();
           itv!=itvend; ++itv)
      {
        d->lcc.mark(itv, markvolumes); // To be sure that all darts of the basic iterator will be marked
        if (!d->lcc.is_marked(itv, markfaces))
        {
          Facet f;
          if(d->compute_face(itv, f))
            d->facets.push_back(f);
          d->facets.back().volume_id = volume_id;
          for (typename LCC::template Dart_of_cell_basic_range<2>::
               const_iterator itf=d->lcc.template darts_of_cell_basic<2>(itv, markfaces).begin(),
               itfend=d->lcc.template darts_of_cell_basic<2>(itv, markfaces).end();
               itf!=itfend; ++itf)
          {
            d->lcc.mark(itf, markfaces); // To be sure that all darts of the basic iterator will be marked
            if ( !d->lcc.is_marked(itf, markedges))
            {
              Point p1 = d->lcc.point(itf);
              LCC::Dart_const_handle d2 = d->lcc.other_extremity(itf);
              Point p2 = d->lcc.point(d2);
              d->lines.push_back(p1.x() + offset.x);
              d->lines.push_back(p1.y() + offset.y);
              d->lines.push_back(p1.z() + offset.z);

              d->lines.push_back(p2.x() + offset.x);
              d->lines.push_back(p2.y() + offset.y);
              d->lines.push_back(p2.z() + offset.z);

              for (typename LCC::template Dart_of_cell_basic_range<1>::
                   const_iterator ite=d->lcc.template darts_of_cell_basic<1>(itf, markedges).begin(),
                   iteend=d->lcc.template darts_of_cell_basic<1>(itf, markedges).end();
                   ite!=iteend; ++ite)
              {
                d->lcc.mark(ite, markedges); // To be sure that all darts of the basic iterator will be marked
                if ( !d->lcc.is_marked(ite, markvertices))
                {
                  Point p1 = d->lcc.point(ite);
                  d->vertices.push_back(p1.x() + offset.x);
                  d->vertices.push_back(p1.y() + offset.y);
                  d->vertices.push_back(p1.z() + offset.z);
                  CGAL::mark_cell<LCC, 0>(d->lcc, ite, markvertices);
                }
              }
            }
          }
        }
      }
      ++volume_id;
    }
  }
  d->nb_volumes = volume_id;
  for (typename LCC::Dart_range::const_iterator it=d->lcc.darts().begin(),
       itend=d->lcc.darts().end(); it!=itend; ++it )
  {
    d->lcc.unmark(it, markvertices);
    d->lcc.unmark(it, markedges);
    d->lcc.unmark(it, markfaces);
    d->lcc.unmark(it, markvolumes);

  }

  d->lcc.free_mark(markvolumes);
  d->lcc.free_mark(markfaces);
  d->lcc.free_mark(markedges);
  d->lcc.free_mark(markvertices);

  getTriangleContainer(0)->allocate(
        Tri::Flat_vertices, d->faces.data(),
        static_cast<int>(d->faces.size()*sizeof(float)));
  if(!d->is_mono_color)
  {
    getTriangleContainer(0)->allocate(Tri::FColors, d->colors.data(),
                                            static_cast<int>(d->colors.size()*sizeof(float)));
  }
  else
    getTriangleContainer(0)->allocate(Tri::FColors, 0, 0);

  getEdgeContainer(0)->allocate(
        Ec::Vertices, d->lines.data(),
        static_cast<int>(d->lines.size()*sizeof(float)));

  getPointContainer(0)->allocate(
        Pc::Vertices, d->vertices.data(),
        static_cast<int>(d->vertices.size()*sizeof(float)));

  setBuffersFilled(true);
  d->nb_faces = d->faces.size();
  d->nb_lines = d->lines.size();
  d->nb_vertices= d->vertices.size();
}
void Scene_lcc_item::initializeBuffers(CGAL::Three::Viewer_interface *viewer) const
{
  getTriangleContainer(0)->initializeBuffers(viewer);
  getTriangleContainer(0)->setFlatDataSize(d->nb_faces);
  getEdgeContainer(0)->initializeBuffers(viewer);
  getEdgeContainer(0)->setFlatDataSize(d->nb_lines);
  getPointContainer(0)->initializeBuffers(viewer);
  getPointContainer(0)->setFlatDataSize(d->nb_vertices);

  d->faces.clear();
  d->faces.shrink_to_fit();
  d->lines.clear();
  d->lines.shrink_to_fit();
  d->vertices.clear();
  d->vertices.shrink_to_fit();
}

void Scene_lcc_item::invalidateOpenGLBuffers()
{
  setBuffersFilled(false);
  getTriangleContainer(0)->reset_vbos(ALL);
  getEdgeContainer(0)->reset_vbos(ALL);
  getPointContainer(0)->reset_vbos(ALL);
  compute_bbox();
}

bool Scene_lcc_item::isEmpty() const
{
  return false;
}

void Scene_lcc_item::randomFaceColors()
{
  d->is_mono_color = false;
  d->colors.resize(d->nb_faces);
  for(std::size_t i=0; i< d->colors.size()-3; i+=3)
  {
    QColor col = generate_random_color();
    d->colors[i] = col.redF();
    d->colors[i+1] = col.greenF();
    d->colors[i+2] = col.blueF();
  }

  invalidateOpenGLBuffers();
  redraw();
}

void Scene_lcc_item::randomVolumeColors()
{
  d->is_mono_color = false;
  d->colors.resize(d->nb_faces);
  std::vector<QColor> colors(d->nb_volumes);
  for(std::size_t i = 0; i<d->nb_volumes; ++i)
  {
    colors[i] = generate_random_color();
  }
  std::size_t color_id = 0;
  for(auto f : d->facets)//filled in the same order as GL faces
  {
    QColor col = colors[f.volume_id];
    //3 points per face.
    for(std::size_t j = 0; j < 3; ++j)
    {
      d->colors[color_id+j*3] = col.redF();
      d->colors[color_id+j*3+1] = col.greenF();
      d->colors[color_id+j*3+2] = col.blueF();
    }
    color_id += 9;
  }
  invalidateOpenGLBuffers();
  redraw();
}

void Scene_lcc_item::resetColors()
{
  d->is_mono_color = true;
  invalidateOpenGLBuffers();
  redraw();
}
QMenu* Scene_lcc_item::contextMenu()
{
  const char* prop_name = "Menu modified by Scene_lcc_item.";

  QMenu* menu = Scene_item::contextMenu();

  // Use dynamic properties:
  // https://doc.qt.io/qt-5/qobject.html#property
  bool menuChanged = menu->property(prop_name).toBool();

  if(!menuChanged) {
    menu->addSeparator();
    QAction* action = menu->addAction(tr("Set Random Colors for Faces."));
    action->setObjectName("actionRandomFaceColors");
    connect(action, &QAction::triggered,
            this, &Scene_lcc_item::randomFaceColors);
    action = menu->addAction(tr("Set Random Colors for Volumes."));
        action->setObjectName("actionRandomVolumeColors");
        connect(action, &QAction::triggered,
                this, &Scene_lcc_item::randomVolumeColors);
    menu->setProperty(prop_name, true);
  }
  QAction* action = menu->findChild<QAction*>("actionResetColors");
  if(!action)
  {
    action = menu->addAction(tr("Reset Colors."));
    action->setObjectName("actionResetColors");
    connect(action, &QAction::triggered,
            this, &Scene_lcc_item::resetColors);
  }
  action->setVisible(!d->is_mono_color);

  return menu;
}
