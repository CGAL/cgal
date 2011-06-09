#include "Scene_nef_polyhedron_item.h"
#include "Scene_polyhedron_item.h"
#include "Nef_type.h"
#include "Polyhedron_type.h"
#include <CGAL/Polyhedron_incremental_builder_3.h>
// #include <CGAL/OFF_to_nef_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#include <CGAL/Inverse_index.h>

#include <QObject>
#include "Scene_nef_rendering.h"

#include <CGAL/minkowski_sum_3.h>
#include <CGAL/convex_decomposition_3.h> 


Scene_nef_polyhedron_item::Scene_nef_polyhedron_item()
  : Scene_item_with_display_list(),
    nef_poly(new Nef_polyhedron)
{
}

Scene_nef_polyhedron_item::Scene_nef_polyhedron_item(Nef_polyhedron* const p)
  : Scene_item_with_display_list(),
    nef_poly(p)
{
}

Scene_nef_polyhedron_item::Scene_nef_polyhedron_item(const Nef_polyhedron& p)
  : Scene_item_with_display_list(),
    nef_poly(new Nef_polyhedron(p))
{
}

// Scene_nef_polyhedron_item::Scene_nef_polyhedron_item(const Scene_nef_polyhedron_item& item)
//   : Scene_item_with_display_list(item),
//     nef_poly(new Nef_polyhedron(*item.nef_poly))
// {
// }

Scene_nef_polyhedron_item::~Scene_nef_polyhedron_item()
{
  delete nef_poly;
}

Scene_nef_polyhedron_item* 
Scene_nef_polyhedron_item::clone() const {
  return new Scene_nef_polyhedron_item(*nef_poly);
}

bool
Scene_nef_polyhedron_item::load_from_off(std::istream& in)
{
//   const std::size_t discarded = CGAL::OFF_to_nef_3(in, *nef_poly);
//   return discarded != 0;

  Exact_polyhedron exact_poly;
  in >> exact_poly;
  *nef_poly = Nef_polyhedron(exact_poly);

//   Polyhedron poly;
//   in >> poly;
//   *nef_poly = Nef_polyhedron(poly);
  return in;
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
  in >> *nef_poly;
  return in;
}

bool
Scene_nef_polyhedron_item::save(std::ostream& in) const
{
  in << *nef_poly;
  return in;
}

QString 
Scene_nef_polyhedron_item::toolTip() const
{
  if(!nef_poly)
    return QString();

  return QObject::tr("<p><b>%1</b> (mode: %5, color: %6)<br />"
                     "<i>Nef_3 polyhedron</i></p>"
                    "<p>Number of vertices: %2<br />"
                     "Number of edges: %3<br />"
                     "Number of facets: %4<br />"
                     "number of volumes: %7</p>")
    .arg(this->name())
    .arg(nef_poly->number_of_vertices())
    .arg(nef_poly->number_of_edges())
    .arg(nef_poly->number_of_facets())
    .arg(this->renderingModeName())
    .arg(this->color().name())
    .arg(nef_poly->number_of_volumes());
}

void
Scene_nef_polyhedron_item::direct_draw() const {
  gl_render_nef_facets(nef_poly);

  GLboolean lighting;
  glGetBooleanv(GL_LIGHTING, &lighting);
  glDisable(GL_LIGHTING);

  GLfloat point_size;
  glGetFloatv(GL_POINT_SIZE, &point_size);
  glPointSize(10.f);

  gl_render_nef_vertices(nef_poly);

  if(lighting) {
    glEnable(GL_LIGHTING);
  }
  glPointSize(point_size);
}

void
Scene_nef_polyhedron_item::draw_edges() const {
  gl_render_nef_edges(nef_poly);
}

Nef_polyhedron* 
Scene_nef_polyhedron_item::nef_polyhedron() {
  return nef_poly;
}

bool
Scene_nef_polyhedron_item::isEmpty() const {
  return (nef_poly == 0) || nef_poly->is_empty();
}

Scene_nef_polyhedron_item::Bbox
Scene_nef_polyhedron_item::bbox() const {
  if(isEmpty())
    return Bbox();
  CGAL::Bbox_3 bbox(nef_poly->vertices_begin()->point().bbox());
  for(Nef_polyhedron::Vertex_const_iterator it = nef_poly->vertices_begin();
      it != nef_poly->vertices_end();
      ++it) {
    bbox = bbox + it->point().bbox();
  }
  return Bbox(bbox.xmin(),bbox.ymin(),bbox.zmin(),
              bbox.xmax(),bbox.ymax(),bbox.zmax());
}

// quick hacks to convert polyhedra from exact to inexact and vice-versa
template <class Polyhedron_input,
class Polyhedron_output>
struct Copy_polyhedron_to 
  : public CGAL::Modifier_base<typename Polyhedron_output::HalfedgeDS>
{
  Copy_polyhedron_to(const Polyhedron_input& in_poly) 
    : in_poly(in_poly) {}

  void operator()(typename Polyhedron_output::HalfedgeDS& out_hds)
  {
    typedef typename Polyhedron_output::HalfedgeDS Output_HDS;
    typedef typename Polyhedron_input::HalfedgeDS Input_HDS;

    CGAL::Polyhedron_incremental_builder_3<Output_HDS> builder(out_hds);

    typedef typename Polyhedron_input::Vertex_const_iterator Vertex_const_iterator;
    typedef typename Polyhedron_input::Facet_const_iterator  Facet_const_iterator;
    typedef typename Polyhedron_input::Halfedge_around_facet_const_circulator HFCC;

    builder.begin_surface(in_poly.size_of_vertices(),
      in_poly.size_of_facets(),
      in_poly.size_of_halfedges());

    for(Vertex_const_iterator
      vi = in_poly.vertices_begin(), end = in_poly.vertices_end();
      vi != end ; ++vi) 
    {
      typename Polyhedron_output::Point_3 p(::CGAL::to_double( vi->point().x()),
	::CGAL::to_double( vi->point().y()),
	::CGAL::to_double( vi->point().z()));
      builder.add_vertex(p);
    }

    typedef CGAL::Inverse_index<Vertex_const_iterator> Index;
    Index index( in_poly.vertices_begin(), in_poly.vertices_end());

    for(Facet_const_iterator 
      fi = in_poly.facets_begin(), end = in_poly.facets_end();
      fi != end; ++fi) 
    {
      HFCC hc = fi->facet_begin();
      HFCC hc_end = hc;
      //     std::size_t n = circulator_size( hc);
      //     CGAL_assertion( n >= 3);
      builder.begin_facet ();
      do {
	builder.add_vertex_to_facet(index[hc->vertex()]);
	++hc;
      } while( hc != hc_end);
      builder.end_facet();
    }
    builder.end_surface();
  } // end operator()(..)
private:
  const Polyhedron_input& in_poly;
}; // end Copy_polyhedron_to<>

template <class Poly_A, class Poly_B>
void copy_to(const Poly_A& poly_a, Poly_B& poly_b)
{
  Copy_polyhedron_to<Poly_A, Poly_B> modifier(poly_a);
  poly_b.delegate(modifier);
}

void from_exact(Exact_polyhedron& in,
		Polyhedron& out)
{
  copy_to(in, out);
  CGAL_assertion(out.is_valid());
}

void to_exact(Polyhedron& in,
	      Exact_polyhedron& out)
{
  copy_to(in, out);
  CGAL_assertion(out.is_valid());
}

bool Scene_nef_polyhedron_item::is_simple() const
{
  return nef_poly->is_simple();
}

// [static]
Scene_nef_polyhedron_item* 
Scene_nef_polyhedron_item::from_polyhedron(Scene_polyhedron_item* item)
{
  Polyhedron* poly = item->polyhedron();
  if(!poly) return 0;

  Exact_polyhedron exact_poly;
  to_exact(*poly, exact_poly);
  Nef_polyhedron* nef_poly = new Nef_polyhedron(exact_poly);
  exact_poly.clear();

  return new Scene_nef_polyhedron_item(nef_poly);
}

Scene_polyhedron_item*
Scene_nef_polyhedron_item::convert_to_polyhedron() const {
  Exact_polyhedron exact_poly;
  nef_poly->convert_to_Polyhedron(exact_poly);
  Polyhedron* poly = new Polyhedron;
  from_exact(exact_poly, *poly);
  exact_poly.clear();
  return new Scene_polyhedron_item(poly);
}

Scene_nef_polyhedron_item&
Scene_nef_polyhedron_item::
operator+=(const Scene_nef_polyhedron_item& other)
{
  (*nef_poly) += (*other.nef_poly);
  return *this;
}

Scene_nef_polyhedron_item&
Scene_nef_polyhedron_item::
operator*=(const Scene_nef_polyhedron_item& other)
{
  (*nef_poly) *= (*other.nef_poly);
  return *this;
}

Scene_nef_polyhedron_item&
Scene_nef_polyhedron_item::
operator-=(const Scene_nef_polyhedron_item& other)
{
  (*nef_poly) -= (*other.nef_poly);
  return *this;
}

Scene_nef_polyhedron_item*
Scene_nef_polyhedron_item::
sum(const Scene_nef_polyhedron_item& a, 
    const Scene_nef_polyhedron_item& b)
{
  return new Scene_nef_polyhedron_item(CGAL::minkowski_sum_3(*a.nef_poly,
                                                             *b.nef_poly));
}

void
Scene_nef_polyhedron_item::
convex_decomposition(std::list< Scene_polyhedron_item*>& convex_parts)
{
  // copy the Nef polyhedron, as markers are added
  Nef_polyhedron N(*nef_poly);
  CGAL::convex_decomposition_3(N);
  
  typedef Nef_polyhedron::Volume_const_iterator Volume_const_iterator;
  
  Volume_const_iterator ci = ++N.volumes_begin();
  for( ; ci != N.volumes_end(); ++ci) {
    if(ci->mark()) {
      Exact_polyhedron P;
      N.convert_inner_shell_to_polyhedron(ci->shells_begin(), P);
      Polyhedron* poly = new Polyhedron;
      from_exact(P, *poly);
      convex_parts.push_back(new Scene_polyhedron_item(poly));
    }
  }
}


#include "Scene_nef_polyhedron_item.moc"
