#include "Scene_nef_polyhedron_item.h"
#include <CGAL/Three/Viewer_interface.h>
#include <CGAL/Three/Three.h>
#include <CGAL/Three/Triangle_container.h>
#include <CGAL/Three/Edge_container.h>
#include <CGAL/Three/Point_container.h>
#include "Scene_surface_mesh_item.h"
#include "Nef_type.h"
#include "SMesh_type.h"
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#include <CGAL/Inverse_index.h>

#include <QObject>
#include <QApplication>

#include <CGAL/minkowski_sum_3.h>
#include <CGAL/convex_decomposition_3.h>

#include <CGAL/Triangulation_2_projection_traits_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/boost/graph/convert_nef_polyhedron_to_polygon_mesh.h>
#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>

using namespace CGAL::Three;
typedef Viewer_interface Vi;
typedef Triangle_container Tc;
typedef Edge_container Ec;
typedef Point_container Pc;

typedef Nef_polyhedron::Traits Traits;
typedef Nef_polyhedron::Halffacet Facet;
typedef CGAL::Triangulation_2_projection_traits_3<Traits>   P_traits;
typedef Nef_polyhedron::Halfedge_const_handle Halfedge_handle;
struct Face_info {
    Nef_polyhedron::Halfedge_const_handle e[3];
    int nesting_level;
};
typedef CGAL::Triangulation_vertex_base_with_info_2<Halfedge_handle,
P_traits>        Vb;
typedef CGAL::Triangulation_face_base_with_info_2<Face_info,
P_traits>          Fb1;
typedef CGAL::Constrained_triangulation_face_base_2<P_traits, Fb1>   Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb>                  TDS;
typedef CGAL::Exact_predicates_tag                                    Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<P_traits,
TDS,
Itag>             CDTbase;
typedef CGAL::Constrained_triangulation_plus_2<CDTbase>              CDT;


struct DPoint {
    DPoint(GLdouble x, GLdouble y, GLdouble z)
    {
        coords[0] = x;
        coords[1] = y;
        coords[2] = z;
    }
    GLdouble coords[3];
};



struct Scene_nef_polyhedron_item_priv
{
  typedef CGAL::Three::Scene_item Base;
  typedef std::vector<QColor> Color_vector;

  Scene_nef_polyhedron_item_priv(Scene_nef_polyhedron_item* parent)
    :nef_poly(new Nef_polyhedron)
  {
    item = parent;
    nb_points = 0;
    nb_lines = 0;
    nb_facets = 0;
  }
  Scene_nef_polyhedron_item_priv(Nef_polyhedron* const p, Scene_nef_polyhedron_item *parent)
    :nef_poly(p)
  {
    item = parent;
    nb_points = 0;
    nb_lines = 0;
    nb_facets = 0;
  }
  Scene_nef_polyhedron_item_priv(const Nef_polyhedron& p, Scene_nef_polyhedron_item* parent)
    :nef_poly(new Nef_polyhedron(p))
  {
    item = parent;
    nb_points = 0;
    nb_lines = 0;
    nb_facets = 0;
  }

~Scene_nef_polyhedron_item_priv()
  {
      delete nef_poly;
  }

  void mark_domains(CDT& ct,
                    CDT::Face_handle start,
                    int index,
                    std::list<CDT::Edge>& border ) const;
  void compute_normals_and_vertices(void) const;
  Nef_polyhedron* nef_poly;


  mutable std::vector<float> positions_lines;
  mutable std::vector<float> positions_facets;
  mutable std::vector<float> positions_points;
  mutable std::vector<float> normals;
  mutable std::vector<float> color_lines;
  mutable std::vector<float> color_facets;
  mutable std::vector<float> color_points;
  mutable std::size_t nb_points;
  mutable std::size_t nb_lines;
  mutable std::size_t nb_facets;
  Scene_nef_polyhedron_item *item;

};

void Scene_nef_polyhedron_item::common_constructor()
{
  setTriangleContainer(0, new Tc(
                         Vi::PROGRAM_WITH_LIGHT, false));
  setEdgeContainer(0, new Ec(
                     Vi::PROGRAM_NO_SELECTION, false));
  setPointContainer(0, new Pc(
                      Vi::PROGRAM_NO_SELECTION, false));
}
Scene_nef_polyhedron_item::Scene_nef_polyhedron_item()
{
  is_selected = true;
  d = new Scene_nef_polyhedron_item_priv(this);
  common_constructor();
}

Scene_nef_polyhedron_item::Scene_nef_polyhedron_item(Nef_polyhedron* const p)
{
    is_selected = true;
    d = new Scene_nef_polyhedron_item_priv(p, this);
    common_constructor();

}

Scene_nef_polyhedron_item::Scene_nef_polyhedron_item(const Nef_polyhedron& p)
{
     is_selected = true;
     d = new Scene_nef_polyhedron_item_priv(p, this);
     common_constructor();
}

Scene_nef_polyhedron_item::~Scene_nef_polyhedron_item()
{
  delete d;
}

void
Scene_nef_polyhedron_item_priv::mark_domains(CDT& ct,
                                             CDT::Face_handle start,
                                             int index,
                                             std::list<CDT::Edge>& border ) const
{
  if(start->info().nesting_level != -1){
    return;
  }
  std::list<CDT::Face_handle> queue;
  queue.push_back(start);
  while(! queue.empty()){
    CDT::Face_handle fh = queue.front();
    queue.pop_front();
    if(fh->info().nesting_level == -1){
      fh->info().nesting_level = index;
      for(int i = 0; i < 3; i++){
        CDT::Edge e(fh,i);
        CDT::Face_handle n = fh->neighbor(i);
        if(n->info().nesting_level == -1){
          if(ct.is_constrained(e)) border.push_back(e);
          else queue.push_back(n);
        }
      }
    }
  }
}


void Scene_nef_polyhedron_item_priv::compute_normals_and_vertices(void) const
{
    QApplication::setOverrideCursor(Qt::WaitCursor);
    int count = 0;
    positions_facets.resize(0);
    positions_points.resize(0);
    normals.resize(0);
    positions_lines.resize(0);
    const CGAL::qglviewer::Vec offset = Three::mainViewer()->offset();

    //The Facets
    {
        for(Nef_polyhedron::Halffacet_const_iterator
            f = nef_poly->halffacets_begin (),
            end = nef_poly->halffacets_end();
            f != end; ++f)
        {
            if(f->is_twin()) continue;
            bool incident_volume_marked = f->incident_volume()->mark();
            count++;
            Nef_polyhedron::Vector_3 v = (incident_volume_marked? -1:1) * f->plane().orthogonal_vector();
            P_traits cdt_traits(v);
            CDT cdt(cdt_traits);

            for(Nef_polyhedron::Halffacet_cycle_const_iterator
                fc = f->facet_cycles_begin(),
                end = f->facet_cycles_end();
                fc != end; ++fc)
            {
                if ( fc.is_shalfedge() )
                {

                    Nef_polyhedron::SHalfedge_const_handle h = fc;
                    Nef_polyhedron::SHalfedge_around_facet_const_circulator hc(h), he(hc);

                    CDT::Vertex_handle previous, first;

                    do {
                        Nef_polyhedron::SVertex_const_handle v = hc->source();
                        const Nef_polyhedron::Point_3& point = v->source()->point();
                        CDT::Vertex_handle vh = cdt.insert(point);
                        if(first == 0) {
                            first = vh;
                        }
                        vh->info() = hc->source();
                        if(previous != 0 && previous != vh) {
                            cdt.insert_constraint(previous, vh);
                        }
                        previous = vh;
                    } while( ++hc != he );

                    cdt.insert_constraint(previous, first);
                }
            }
            // sets mark is_external
            for(CDT::All_faces_iterator
                fit = cdt.all_faces_begin(),
                end = cdt.all_faces_end();
                fit != end; ++fit)
            {
                fit->info().nesting_level = -1;

            }

            //check if the facet is external or internal
            std::queue<CDT::Face_handle> face_queue;
            face_queue.push(cdt.infinite_vertex()->face());

            std::list<CDT::Edge> border;
            mark_domains(cdt, cdt.infinite_face(), 0, border);
            while(! border.empty()){
              CDT::Edge e = border.front();
              border.pop_front();
              CDT::Face_handle n = e.first->neighbor(e.second);
              if(n->info().nesting_level == -1){
                mark_domains(cdt, n, e.first->info().nesting_level+1, border);
              }
            }

            //iterates on the internal faces to add the vertices to the positions
            //and the normals to the appropriate vectors

            for(CDT::Finite_faces_iterator
                ffit = cdt.finite_faces_begin(),
                end = cdt.finite_faces_end();
                ffit != end; ++ffit)
            {


                if(ffit->info().nesting_level%2 != 1){ continue;}
                for(int i = 0; i<3; i++)
                {
                    positions_facets.push_back(CGAL::to_double(ffit->vertex(i)->point().x())+offset.x);
                    positions_facets.push_back(CGAL::to_double(ffit->vertex(i)->point().y())+offset.y);
                    positions_facets.push_back(CGAL::to_double(ffit->vertex(i)->point().z())+offset.z);

                }

                GLfloat normal[3];
                normal[0] = CGAL::to_double(v.x());
                normal[1] = CGAL::to_double(v.y());
                normal[2] = CGAL::to_double(v.z());
                GLfloat norm = normal[0]*normal[0]
                        + normal[1]*normal[1]
                        + normal[2]*normal[2];
                norm = CGAL::sqrt(norm);
                normal[0] /= norm;
                normal[1] /= norm;
                normal[2] /= norm;

                normals.push_back(normal[0]);
                normals.push_back(normal[1]);
                normals.push_back(normal[2]);

                normals.push_back(normal[0]);
                normals.push_back(normal[1]);
                normals.push_back(normal[2]);

                normals.push_back(normal[0]);
                normals.push_back(normal[1]);
                normals.push_back(normal[2]);

            }
        }

    } // end facets

    //The Lines
    {
       for(Nef_polyhedron::Halfedge_const_iterator
            e = nef_poly->halfedges_begin(),
            end = nef_poly->halfedges_end();
            e != end; ++e)
        {
            if (e->is_twin()) continue;
            const Nef_polyhedron::Vertex_const_handle& s = e->source();
            const Nef_polyhedron::Vertex_const_handle& t = e->twin()->source();
            const Nef_polyhedron::Point_3& a = s->point();
            const Nef_polyhedron::Point_3& b = t->point();

            positions_lines.push_back(CGAL::to_double(a.x())+offset.x);
            positions_lines.push_back(CGAL::to_double(a.y())+offset.y);
            positions_lines.push_back(CGAL::to_double(a.z())+offset.z);

            positions_lines.push_back(CGAL::to_double(b.x())+offset.x);
            positions_lines.push_back(CGAL::to_double(b.y())+offset.y);
            positions_lines.push_back(CGAL::to_double(b.z())+offset.z);


        }
    }
    //The points
    {
        for(Nef_polyhedron::Vertex_const_iterator
            v = nef_poly->vertices_begin(),
            end = nef_poly->vertices_end();
            v != end; ++v)
        {
            const Nef_polyhedron::Point_3& p = v->point();
            positions_points.push_back(CGAL::to_double(p.x())+offset.x);
            positions_points.push_back(CGAL::to_double(p.y())+offset.y);
            positions_points.push_back(CGAL::to_double(p.z())+offset.z);

                color_points.push_back(item->color().lighter(50).redF());
                color_points.push_back(item->color().lighter(50).greenF());
                color_points.push_back(item->color().lighter(50).blueF());

                color_points.push_back(item->color().lighter(50).redF());
                color_points.push_back(item->color().lighter(50).greenF());
                color_points.push_back(item->color().lighter(50).blueF());

        }

    } //end points
    QApplication::restoreOverrideCursor();
}
Scene_nef_polyhedron_item*
Scene_nef_polyhedron_item::clone() const {
    return new Scene_nef_polyhedron_item(*d->nef_poly);
}

bool
Scene_nef_polyhedron_item::load_from_off(std::istream& in)
{
    //   const std::size_t discarded = CGAL::OFF_to_nef_3(in, *nef_poly);
    //   return discarded != 0;

    Exact_polyhedron exact_poly;
    in >> exact_poly;
    *d->nef_poly = Nef_polyhedron(exact_poly);

    //   Polyhedron poly;
    //   in >> poly;
    //   *nef_poly = Nef_polyhedron(poly);
    invalidateOpenGLBuffers();
    return (bool) in;
}

QFont
Scene_nef_polyhedron_item::font() const {
    QFont font;
    font.setItalic(!font.italic());
    return font;
}

bool
Scene_nef_polyhedron_item::load(std::istream& in)
{
    in >> *d->nef_poly;
    invalidateOpenGLBuffers();
    return (bool) in;
}

bool
Scene_nef_polyhedron_item::save(std::ostream& in) const
{
    in << *d->nef_poly;
    return (bool) in;
}

QString
Scene_nef_polyhedron_item::toolTip() const
{
    if(!d->nef_poly)
        return QString();

    return QObject::tr("<p><b>%1</b> (mode: %5, color: %6)<br />"
                       "<i>Nef_polyhedron_3</i></p>"
                       "<p>Number of vertices: %2<br />"
                       "Number of edges: %3<br />"
                       "Number of facets: %4<br />"
                       "number of volumes: %7</p>")
            .arg(this->name())
            .arg(d->nef_poly->number_of_vertices())
            .arg(d->nef_poly->number_of_edges())
            .arg(d->nef_poly->number_of_facets())
            .arg(this->renderingModeName())
            .arg(this->color().name())
            .arg(d->nef_poly->number_of_volumes());
}

void Scene_nef_polyhedron_item::draw(CGAL::Three::Viewer_interface* viewer) const
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

  Tc* tc = getTriangleContainer(0);
  tc->setColor(this->color());
  tc->getVao(viewer)->bind();
  tc->getVao(viewer)->program->setUniformValue("is_two_side", 1);
  tc->getVao(viewer)->release();
  tc->draw(viewer, true);
  drawPoints(viewer);

}
void Scene_nef_polyhedron_item::drawEdges(CGAL::Three::Viewer_interface* viewer) const
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
  Ec* ec = getEdgeContainer(0);
  ec->setColor(QColor(Qt::black));
  ec->draw(viewer, true);
  if(renderingMode() == PointsPlusNormals)
  {
    drawPoints(viewer);
  }
}

void Scene_nef_polyhedron_item::drawPoints(CGAL::Three::Viewer_interface* viewer) const
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
  GLfloat point_size;
  viewer->glGetFloatv(GL_POINT_SIZE, &point_size);
  viewer->setGlPointSize(10.f);
  Pc* pc = getPointContainer(0);
  pc->setColor(this->color().lighter(50));
  pc->draw(viewer, true);
  viewer->setGlPointSize(point_size);
}

Nef_polyhedron*
Scene_nef_polyhedron_item::nef_polyhedron() {
    return d->nef_poly;
}

Nef_polyhedron*
Scene_nef_polyhedron_item::nef_polyhedron()const {
    return d->nef_poly;
}

bool
Scene_nef_polyhedron_item::isEmpty() const {
    return (d->nef_poly == 0) || d->nef_poly->is_empty();
}

void
Scene_nef_polyhedron_item::compute_bbox() const {
    if(isEmpty())
    {
        setBbox(Bbox());
        return;
    }
    CGAL::Bbox_3 bbox(d->nef_poly->vertices_begin()->point().bbox());
    for(Nef_polyhedron::Vertex_const_iterator it = d->nef_poly->vertices_begin();
        it != d->nef_poly->vertices_end();
        ++it) {
        bbox = bbox + it->point().bbox();
    }
    setBbox(Bbox(bbox.xmin(),bbox.ymin(),bbox.zmin(),
                bbox.xmax(),bbox.ymax(),bbox.zmax()));
}

template <class Poly_A, class Poly_B>
void copy_to(const Poly_A& poly_a, Poly_B& poly_b)
{
  CGAL::copy_face_graph(poly_a, poly_b);
}

void from_exact(Exact_polyhedron& in,
                SMesh& out)
{
    copy_to(in, out);
    CGAL_assertion(out.is_valid());
}

void to_exact(SMesh& in,
              Exact_polyhedron& out)
{
    copy_to(in, out);
    CGAL_assertion(out.is_valid());
}

bool Scene_nef_polyhedron_item::is_simple() const
{
    return d->nef_poly->is_simple();
}
template<typename FaceGraph>
struct Halfedge_index_pmap
{

};

template<typename FaceGraph>
struct Face_index_pmap
{

};
// [static]
Scene_nef_polyhedron_item*
Scene_nef_polyhedron_item::from_polygon_mesh(Scene_surface_mesh_item* item)
{
  SMesh* sm = item->polyhedron();
  if(!sm) return 0;
  CGAL::Surface_mesh<Exact_Kernel::Point_3> exact_sm;
  CGAL::copy_face_graph(*sm, exact_sm);
  Nef_polyhedron* nef_poly = new Nef_polyhedron(exact_sm);

  return new Scene_nef_polyhedron_item(nef_poly);
}

Scene_surface_mesh_item* Scene_nef_polyhedron_item::convert_to_surface_mesh() const
{
  SMesh* poly = new SMesh();
  CGAL::convert_nef_polyhedron_to_polygon_mesh(*this->nef_polyhedron(), *poly);
  CGAL::Polygon_mesh_processing::triangulate_faces(*poly);
  return new Scene_surface_mesh_item(poly);
}

Scene_nef_polyhedron_item&
Scene_nef_polyhedron_item::
operator+=(const Scene_nef_polyhedron_item& other)
{
    (*d->nef_poly) += (*other.d->nef_poly);
    return *this;
}

Scene_nef_polyhedron_item&
Scene_nef_polyhedron_item::
operator*=(const Scene_nef_polyhedron_item& other)
{
    (*d->nef_poly) *= (*other.d->nef_poly);
    return *this;
}

Scene_nef_polyhedron_item&
Scene_nef_polyhedron_item::
operator-=(const Scene_nef_polyhedron_item& other)
{
    (*d->nef_poly) -= (*other.d->nef_poly);
    return *this;
}

Scene_nef_polyhedron_item*
Scene_nef_polyhedron_item::
sum(const Scene_nef_polyhedron_item& a,
    const Scene_nef_polyhedron_item& b)
{
    return new Scene_nef_polyhedron_item(CGAL::minkowski_sum_3(*a.d->nef_poly,
                                                               *b.d->nef_poly));
}

void
Scene_nef_polyhedron_item::
convex_decomposition(std::list< Scene_surface_mesh_item*>& convex_parts)
{
    // copy the Nef polyhedron, as markers are added
    Nef_polyhedron N(*d->nef_poly);
    CGAL::convex_decomposition_3(N);

    typedef Nef_polyhedron::Volume_const_iterator Volume_const_iterator;

    Volume_const_iterator ci = ++N.volumes_begin();
    for( ; ci != N.volumes_end(); ++ci) {
        if(ci->mark()) {
            Exact_polyhedron P;
            N.convert_inner_shell_to_polyhedron(ci->shells_begin(), P);
            SMesh* poly = new SMesh;
            from_exact(P, *poly);
            Scene_surface_mesh_item *spoly = new Scene_surface_mesh_item(poly);
            convex_parts.push_back(spoly);
            spoly->invalidateOpenGLBuffers();
        }
    }
}

void
Scene_nef_polyhedron_item::
invalidateOpenGLBuffers()
{
    compute_bbox();
    setBuffersFilled(false);
    getTriangleContainer(0)->reset_vbos(ALL);
    getEdgeContainer(0)->reset_vbos(ALL);
    getPointContainer(0)->reset_vbos(ALL);
}

void
Scene_nef_polyhedron_item::selection_changed(bool p_is_selected)
{

    if(p_is_selected != is_selected)
        is_selected = p_is_selected;

}

void Scene_nef_polyhedron_item::computeElements() const
{
  d->compute_normals_and_vertices();

  Tc* tc = getTriangleContainer(0);
  tc->allocate(
        Tc::Flat_vertices,
        d->positions_facets.data(),
        static_cast<int>(d->positions_facets.size()*sizeof(float)));

  tc->allocate(
        Tc::Flat_normals,
        d->normals.data(),
        static_cast<int>(d->normals.size()*sizeof(float)));

  d->nb_facets = d->positions_facets.size();

  getEdgeContainer(0)->allocate(
        Ec::Vertices,
        d->positions_lines.data(),
        static_cast<int>(d->positions_lines.size()*sizeof(float)));
  d->nb_lines = d->positions_lines.size();


  getPointContainer(0)->allocate(
        Pc::Vertices,
        d->positions_points.data(),
        static_cast<int>(d->positions_points.size()*sizeof(float)));
  d->nb_points = d->positions_points.size();
  setBuffersFilled(true);
}

void Scene_nef_polyhedron_item::initializeBuffers(Viewer_interface *v) const
{
  getTriangleContainer(0)->initializeBuffers(v);
  getEdgeContainer(0)->initializeBuffers(v);
  getPointContainer(0)->initializeBuffers(v);
  getTriangleContainer(0)->setFlatDataSize(d->nb_facets);
  getEdgeContainer(0)->setFlatDataSize(d->nb_lines);
  getPointContainer(0)->setFlatDataSize(d->nb_points);
  d->positions_facets.clear();
  d->positions_facets.shrink_to_fit();
  d->normals.clear();
  d->normals.shrink_to_fit();
  d->positions_lines.clear();
  d->positions_lines.shrink_to_fit();
  d->positions_points.clear();
  d->positions_points.shrink_to_fit();
}
