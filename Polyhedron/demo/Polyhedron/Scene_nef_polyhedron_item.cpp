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

#include <CGAL/minkowski_sum_3.h>
#include <CGAL/convex_decomposition_3.h>

#include <CGAL/Triangulation_2_filtered_projection_traits_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>

typedef Nef_polyhedron::Traits Traits;
typedef Nef_polyhedron::Halffacet Facet;
typedef CGAL::Triangulation_2_filtered_projection_traits_3<Traits>   P_traits;
typedef Nef_polyhedron::Halfedge_const_handle Halfedge_handle;
struct Face_info {
    Nef_polyhedron::Halfedge_const_handle e[3];
    bool is_external;
};
typedef CGAL::Triangulation_vertex_base_with_info_2<Halfedge_handle,
P_traits>        Vb;
typedef CGAL::Triangulation_face_base_with_info_2<Face_info,
P_traits>          Fb1;
typedef CGAL::Constrained_triangulation_face_base_2<P_traits, Fb1>   Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb>                  TDS;
typedef CGAL::No_intersection_tag                                    Itag;
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
Scene_nef_polyhedron_item::Scene_nef_polyhedron_item()
    : Scene_item(7,3),
      nef_poly(new Nef_polyhedron)
{
    is_selected = true;
    nb_points = 0;
    nb_lines = 0;
    nb_facets = 0;
}

Scene_nef_polyhedron_item::Scene_nef_polyhedron_item(Nef_polyhedron* const p)
    : Scene_item(7,3),
      nef_poly(p)
{
    is_selected = true;
    nb_points = 0;
    nb_lines = 0;
    nb_facets = 0;
}

Scene_nef_polyhedron_item::Scene_nef_polyhedron_item(const Nef_polyhedron& p)
    : Scene_item(7,3),
      nef_poly(new Nef_polyhedron(p))
{
     is_selected = true;
     nb_points = 0;
     nb_lines = 0;
     nb_facets = 0;
}

Scene_nef_polyhedron_item::~Scene_nef_polyhedron_item()
{
    delete nef_poly;
}

void Scene_nef_polyhedron_item::initialize_buffers(Viewer_interface *viewer) const
{
    //vao for the facets
    {
        program = getShaderProgram(PROGRAM_WITH_LIGHT, viewer);
        program->bind();

        vaos[0]->bind();
        buffers[0].bind();
        buffers[0].allocate(positions_facets.data(),
                            static_cast<int>(positions_facets.size()*sizeof(double)));
        program->enableAttributeArray("vertex");
        program->setAttributeBuffer("vertex",GL_DOUBLE,0,3);
        buffers[0].release();



        buffers[1].bind();
        buffers[1].allocate(normals.data(),
                            static_cast<int>(normals.size()*sizeof(double)));
        program->enableAttributeArray("normals");
        program->setAttributeBuffer("normals",GL_DOUBLE,0,3);
        buffers[1].release();

        buffers[2].bind();
        buffers[2].allocate(color_facets.data(),
                            static_cast<int>(color_facets.size()*sizeof(double)));
        program->enableAttributeArray("colors");
        program->setAttributeBuffer("colors",GL_DOUBLE,0,3);
        buffers[2].release();
        vaos[0]->release();
        nb_facets = positions_facets.size();
        positions_facets.resize(0);
        std::vector<double>(positions_facets).swap(positions_facets);

        normals.resize(0);
        std::vector<double>(normals).swap(normals);

        color_facets.resize(0);
        std::vector<double>(color_facets).swap(color_facets);
        program->release();

    }
    //vao for the edges
    {
        program = getShaderProgram(PROGRAM_WITHOUT_LIGHT, viewer);
        program->bind();

        vaos[1]->bind();
        buffers[3].bind();
        buffers[3].allocate(positions_lines.data(),
                            static_cast<int>(positions_lines.size()*sizeof(double)));
        program->enableAttributeArray("vertex");
        program->setAttributeBuffer("vertex",GL_DOUBLE,0,3);
        buffers[3].release();



        buffers[4].bind();
        buffers[4].allocate(color_lines.data(),
                            static_cast<int>(color_lines.size()*sizeof(double)));
        program->enableAttributeArray("colors");
        program->setAttributeBuffer("colors",GL_DOUBLE,0,3);
        buffers[4].release();

        nb_lines = positions_lines.size();
        positions_lines.resize(0);
        std::vector<double>(positions_lines).swap(positions_lines);

        color_lines.resize(0);
        std::vector<double>(color_lines).swap(color_lines);
        vaos[1]->release();
        program->release();
    }
    //vao for the points
    {
        program = getShaderProgram(PROGRAM_WITHOUT_LIGHT, viewer);
        program->bind();

        vaos[2]->bind();
        buffers[5].bind();
        buffers[5].allocate(positions_points.data(),
                            static_cast<int>(positions_points.size()*sizeof(double)));
        program->enableAttributeArray("vertex");
        program->setAttributeBuffer("vertex",GL_DOUBLE,0,3);
        buffers[5].release();



        buffers[6].bind();
        buffers[6].allocate(color_points.data(),
                            static_cast<int>(color_points.size()*sizeof(double)));
        program->enableAttributeArray("colors");
        program->setAttributeBuffer("colors",GL_DOUBLE,0,3);
        buffers[6].release();

        vaos[2]->release();

        nb_points = positions_points.size();
        positions_points.resize(0);
        std::vector<double>(positions_points).swap(positions_points);

        color_points.resize(0);
        std::vector<double>(color_points).swap(color_points);
        program->release();
    }
    are_buffers_filled = true;
}
void Scene_nef_polyhedron_item::compute_normals_and_vertices(void)
{
     int count = 0;
    positions_facets.resize(0);
    positions_points.resize(0);
    color_lines.resize(0);
    color_facets.resize(0);
    color_points.resize(0);
    normals.resize(0);
    positions_lines.resize(0);
    //The Facets
    {
        for(Nef_polyhedron::Halffacet_const_iterator
            f = nef_poly->halffacets_begin (),
            end = nef_poly->halffacets_end();
            f != end; ++f)
        {
            if(f->is_twin()) continue;
            count++;
            Nef_polyhedron::Vector_3 v = f->plane().orthogonal_vector();
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

                    // sets mark is_external
                    for(CDT::All_faces_iterator
                        fit = cdt.all_faces_begin(),
                        end = cdt.all_faces_end();
                        fit != end; ++fit)
                    {
                        fit->info().is_external = false;

                    }
                    //check if the facet is external or internal
                    std::queue<CDT::Face_handle> face_queue;
                    face_queue.push(cdt.infinite_vertex()->face());

                    while(! face_queue.empty() ) {
                        CDT::Face_handle fh = face_queue.front();
                        face_queue.pop();
                        if(fh->info().is_external) continue;
                        fh->info().is_external = true;
                        for(int i = 0; i <3; ++i) {
                            if(!cdt.is_constrained(std::make_pair(fh, i)))
                            {
                                face_queue.push(fh->neighbor(i));
                            }
                        }

                    }
                    //iterates on the internal faces to add the vertices to the positions
                    //and the normals to the appropriate vectors

                    for(CDT::Finite_faces_iterator
                        ffit = cdt.finite_faces_begin(),
                        end = cdt.finite_faces_end();
                        ffit != end; ++ffit)
                    {


                        if(ffit->info().is_external){ continue;}
                        for(int i = 0; i<3; i++)
                        {
                            positions_facets.push_back(CGAL::to_double(ffit->vertex(i)->point().x()));
                            positions_facets.push_back(CGAL::to_double(ffit->vertex(i)->point().y()));
                            positions_facets.push_back(CGAL::to_double(ffit->vertex(i)->point().z()));

                        }



                        Nef_polyhedron::Vector_3 v = f->plane().orthogonal_vector();
                        GLdouble normal[3];
                        normal[0] = CGAL::to_double(v.x());
                        normal[1] = CGAL::to_double(v.y());
                        normal[2] = CGAL::to_double(v.z());
                        GLdouble norm = normal[0]*normal[0]
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

                        if(is_selected)
                        {
                            color_facets.push_back(this->color().lighter(120).redF());
                            color_facets.push_back(this->color().lighter(120).greenF());
                            color_facets.push_back(this->color().lighter(120).blueF());

                            color_facets.push_back(this->color().lighter(120).redF());
                            color_facets.push_back(this->color().lighter(120).greenF());
                            color_facets.push_back(this->color().lighter(120).blueF());

                            color_facets.push_back(this->color().lighter(120).redF());
                            color_facets.push_back(this->color().lighter(120).greenF());
                            color_facets.push_back(this->color().lighter(120).blueF());
                        }
                        else
                        {
                            color_facets.push_back(this->color().redF());
                            color_facets.push_back(this->color().greenF());
                            color_facets.push_back(this->color().blueF());

                            color_facets.push_back(this->color().redF());
                            color_facets.push_back(this->color().greenF());
                            color_facets.push_back(this->color().blueF());

                            color_facets.push_back(this->color().redF());
                            color_facets.push_back(this->color().greenF());
                            color_facets.push_back(this->color().blueF());

                        }

                    }
                }
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

            positions_lines.push_back(CGAL::to_double(a.x()));
            positions_lines.push_back(CGAL::to_double(a.y()));
            positions_lines.push_back(CGAL::to_double(a.z()));

            positions_lines.push_back(CGAL::to_double(b.x()));
            positions_lines.push_back(CGAL::to_double(b.y()));
            positions_lines.push_back(CGAL::to_double(b.z()));

            if(is_selected)
            {
                color_lines.push_back(this->color().lighter(50).redF());
                color_lines.push_back(this->color().lighter(50).greenF());
                color_lines.push_back(this->color().lighter(50).blueF());

                color_lines.push_back(this->color().lighter(50).redF());
                color_lines.push_back(this->color().lighter(50).greenF());
                color_lines.push_back(this->color().lighter(50).blueF());
            }
            else
            {
                color_lines.push_back(0.0);
                color_lines.push_back(0.0);
                color_lines.push_back(0.0);

                color_lines.push_back(0.0);
                color_lines.push_back(0.0);
                color_lines.push_back(0.0);
            }
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
            positions_points.push_back(CGAL::to_double(p.x()));
            positions_points.push_back(CGAL::to_double(p.y()));
            positions_points.push_back(CGAL::to_double(p.z()));

                color_points.push_back(this->color().lighter(50).redF());
                color_points.push_back(this->color().lighter(50).greenF());
                color_points.push_back(this->color().lighter(50).blueF());

                color_points.push_back(this->color().lighter(50).redF());
                color_points.push_back(this->color().lighter(50).greenF());
                color_points.push_back(this->color().lighter(50).blueF());

        }

    } //end points
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
    invalidate_buffers();
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
    in >> *nef_poly;
    invalidate_buffers();
    return (bool) in;
}

bool
Scene_nef_polyhedron_item::save(std::ostream& in) const
{
    in << *nef_poly;
    return (bool) in;
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

void Scene_nef_polyhedron_item::draw(Viewer_interface* viewer) const
{
    if(!are_buffers_filled)
        initialize_buffers(viewer);
    vaos[0]->bind();

    // tells the GPU to use the program just created
    program=getShaderProgram(PROGRAM_WITH_LIGHT);
    attrib_buffers(viewer,PROGRAM_WITH_LIGHT);
    program->bind();
    viewer->glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(nb_facets/3));
    vaos[0]->release();
    program->release();
    GLfloat point_size;
    viewer->glGetFloatv(GL_POINT_SIZE, &point_size);
    viewer->glPointSize(10.f);

    draw_points(viewer);
    viewer->glPointSize(point_size);

}
void Scene_nef_polyhedron_item::draw_edges(Viewer_interface* viewer) const
{
    if(!are_buffers_filled)
        initialize_buffers(viewer);

    vaos[1]->bind();
    program = getShaderProgram(PROGRAM_WITHOUT_LIGHT);
    attrib_buffers(viewer ,PROGRAM_WITHOUT_LIGHT);
    program->bind();
    viewer->glDrawArrays(GL_LINES,0,static_cast<GLsizei>(nb_lines/3));
    vaos[1]->release();
    program->release();
    if(renderingMode() == PointsPlusNormals)
    {
        GLfloat point_size;
        viewer->glGetFloatv(GL_POINT_SIZE, &point_size);
        viewer->glPointSize(10.f);

        draw_points(viewer);
        viewer->glPointSize(point_size);
    }
}
void Scene_nef_polyhedron_item::draw_points(Viewer_interface* viewer) const
{
    if(!are_buffers_filled)
        initialize_buffers(viewer);
    vaos[2]->bind();
    program=getShaderProgram(PROGRAM_WITHOUT_LIGHT);
    attrib_buffers(viewer ,PROGRAM_WITHOUT_LIGHT);
    program->bind();
    viewer->glDrawArrays(GL_POINTS,0,static_cast<GLsizei>(nb_points/3));
    vaos[2]->release();
    program->release();

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
            Scene_polyhedron_item *spoly = new Scene_polyhedron_item(poly);
            convex_parts.push_back(spoly);
            spoly->invalidate_buffers();
        }
    }
}

void
Scene_nef_polyhedron_item::
invalidate_buffers()
{
    //  init();
    Base::invalidate_buffers();
    compute_normals_and_vertices();
    are_buffers_filled = false;
}
void
Scene_nef_polyhedron_item::selection_changed(bool p_is_selected)
{

    if(p_is_selected != is_selected)
    {
        is_selected = p_is_selected;
        invalidate_buffers();
    }

}
