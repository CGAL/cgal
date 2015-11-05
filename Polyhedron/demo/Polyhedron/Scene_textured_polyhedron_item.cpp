#include "Scene_textured_polyhedron_item.h"
#include "Textured_polyhedron_type.h"
#include <CGAL/IO/Polyhedron_iostream.h>

#include <QObject>
#include <CGAL/gl_render.h>

typedef EPIC_kernel::Point_3 Point;


void Scene_textured_polyhedron_item::initialize_buffers(Viewer_interface *viewer = 0) const
{
    if(GLuint(-1) == textureId) {
        viewer->glGenTextures(1, &textureId);
    }
    //vao for the facets
    {
        program = getShaderProgram(PROGRAM_WITH_TEXTURE, viewer);
        program->bind();
        vaos[0]->bind();
        buffers[0].bind();
        buffers[0].allocate(positions_facets.data(),
                            static_cast<int>(positions_facets.size()*sizeof(float)));
        program->enableAttributeArray("vertex");
        program->setAttributeBuffer("vertex",GL_FLOAT,0,4);
        buffers[0].release();

        buffers[1].bind();
        buffers[1].allocate(normals.data(),
                            static_cast<int>(normals.size()*sizeof(float)));
        program->enableAttributeArray("normal");
        program->setAttributeBuffer("normal",GL_FLOAT,0,3);
        buffers[1].release();


        buffers[2].bind();
        buffers[2].allocate(textures_map_facets.data(),
                            static_cast<int>(textures_map_facets.size()*sizeof(float)));
        program->enableAttributeArray("v_texCoord");
        program->setAttributeBuffer("v_texCoord",GL_FLOAT,0,2);
        buffers[2].release();
        vaos[0]->release();
        program->release();
    }

    //vao for the lines
    {
        program = getShaderProgram(PROGRAM_WITH_TEXTURED_EDGES, viewer);
        program->bind();
        vaos[1]->bind();
        buffers[3].bind();
        buffers[3].allocate(positions_lines.data(),
                            static_cast<int>(positions_lines.size()*sizeof(float)));
        program->enableAttributeArray("vertex");
        program->setAttributeBuffer("vertex",GL_FLOAT,0,4);
        buffers[3].release();


        buffers[4].bind();
        buffers[4].allocate(textures_map_lines.data(), 
                            static_cast<int>(textures_map_lines.size()*sizeof(float)));
        program->enableAttributeArray("v_texCoord");
        program->setAttributeBuffer("v_texCoord",GL_FLOAT,0,2);
        buffers[4].release();
        vaos[1]->release();
        program->release();
    }

    viewer->glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

    viewer->glActiveTexture(GL_TEXTURE0);
    viewer->glBindTexture(GL_TEXTURE_2D, textureId);
    viewer->glTexImage2D(GL_TEXTURE_2D,
                 0,
                 GL_RGB,
                 texture.GetWidth(),
                 texture.GetHeight(),
                 0,
                 GL_RGB,
                 GL_UNSIGNED_BYTE,
                 texture.GetData());
    viewer->glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    viewer->glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    viewer->glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    viewer->glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

    nb_facets = positions_facets.size();
    positions_facets.resize(0);
    std::vector<float>(positions_facets).swap(positions_facets);
    normals.resize(0);
    std::vector<float>(normals).swap(normals);
    textures_map_facets.resize(0);
    std::vector<float>(textures_map_facets).swap(textures_map_facets);
    nb_lines = positions_lines.size();
    positions_lines.resize(0);
    std::vector<float>(positions_lines).swap(positions_lines);
    textures_map_lines.resize(0);
    std::vector<float>(textures_map_lines).swap(textures_map_lines);

    are_buffers_filled = true;
}

void
Scene_textured_polyhedron_item::compute_normals_and_vertices(void)
{
    positions_facets.resize(0);
    positions_lines.resize(0);
    textures_map_facets.resize(0);
    textures_map_lines.resize(0);
    normals.resize(0);

    typedef ::EPIC_kernel Kernel;
    typedef CGAL::Textured_items Items;
    typedef Kernel::Point_3 Point;
    typedef Kernel::Vector_3 Vector;
    typedef CGAL::Polyhedron_3<Kernel,Items> Base;

    typedef Base::Halfedge_around_facet_circulator Halfedge_around_facet_circulator;
    typedef Base::Edge_iterator Edge_iterator;
    typedef Base::Facet Facet;
    typedef Base::Facet_iterator Facet_iterator;

    //Facets
    Facet_iterator f = poly->facets_begin();

    for(f = poly->facets_begin();
        f != poly->facets_end();
        f++)
    {

        Halfedge_around_facet_circulator he = f->facet_begin();
        Halfedge_around_facet_circulator end = he;
        CGAL_For_all(he,end)
        {

            // If Flat shading:1 normal per polygon added once per vertex
            if (cur_shading == Flat || cur_shading == FlatPlusEdges)
            {

                Vector n = CGAL::Polygon_mesh_processing::
                  compute_face_normal(f, static_cast<Base&>(*poly));
                normals.push_back(n[0]);
                normals.push_back(n[1]);
                normals.push_back(n[2]);
            }

            // If Gouraud shading: 1 normal per vertex
            else if(cur_shading == Gouraud)
            {

                const Facet::Normal_3& n = he->vertex()->normal();
                normals.push_back(n[0]);
                normals.push_back(n[1]);
                normals.push_back(n[2]);
            }

            //position
            const Point& p = he->vertex()->point();
            positions_facets.push_back(p.x());
            positions_facets.push_back(p.y());
            positions_facets.push_back(p.z());
            positions_facets.push_back(1.0);

            const double u = he->vertex()->u();
            const double v = he->vertex()->v();
            textures_map_facets.push_back(u);
            textures_map_facets.push_back(v);
        }


    }
    //Lines
    typedef Kernel::Point_3		Point;
    typedef Base::Edge_iterator	Edge_iterator;

    Edge_iterator he;

    for(he = poly->edges_begin();
        he != poly->edges_end();
        he++)
    {

        const Point& a = he->vertex()->point();
        const Point& b = he->opposite()->vertex()->point();
        positions_lines.push_back(a.x());
        positions_lines.push_back(a.y());
        positions_lines.push_back(a.z());
        positions_lines.push_back(1.0);

        const double u = he->vertex()->u();
        const double v = he->vertex()->v();
        textures_map_lines.push_back(u);
        textures_map_lines.push_back(v);

        positions_lines.push_back(b.x());
        positions_lines.push_back(b.y());
        positions_lines.push_back(b.z());
        positions_lines.push_back(1.0);

        const double ou = he->opposite()->vertex()->u();
        const double ov = he->opposite()->vertex()->v();
        textures_map_lines.push_back(ou);
        textures_map_lines.push_back(ov);

    }

}

Scene_textured_polyhedron_item::Scene_textured_polyhedron_item()
    : Scene_item(5,2),poly(new Textured_polyhedron), textureId(-1)
{
    texture.GenerateCheckerBoard(2048,2048,128,0,0,0,250,250,255);
    cur_shading=FlatPlusEdges;
    is_selected=false;
    nb_facets = 0;
    nb_lines = 0;
    invalidate_buffers();
}

Scene_textured_polyhedron_item::Scene_textured_polyhedron_item(Textured_polyhedron* const p)
    : Scene_item(5,2),poly(p),textureId(-1),smooth_shading(true)
{
    cur_shading=FlatPlusEdges;
    is_selected=false;
    texture.GenerateCheckerBoard(2048,2048,128,0,0,0,250,250,255);
    nb_facets = 0;
    nb_lines = 0;
    invalidate_buffers();
}

Scene_textured_polyhedron_item::Scene_textured_polyhedron_item(const Textured_polyhedron& p)
    : Scene_item(5,2), poly(new Textured_polyhedron(p)),textureId(-1),smooth_shading(true)
{
    texture.GenerateCheckerBoard(2048,2048,128,0,0,0,250,250,255);
    cur_shading=FlatPlusEdges;
    is_selected=false;
    nb_facets = 0;
    nb_lines = 0;
    invalidate_buffers();
}

Scene_textured_polyhedron_item::~Scene_textured_polyhedron_item()
{
    delete poly;
}

Scene_textured_polyhedron_item* 
Scene_textured_polyhedron_item::clone() const {
    return new Scene_textured_polyhedron_item(*poly);
}

// Load textured_polyhedron from .OFF file
bool
Scene_textured_polyhedron_item::load(std::istream& in)
{
    std::cout<<"LOAD"<<std::endl;
    in >> *poly;
    invalidate_buffers();
    return in && !isEmpty();
}

// Write textured_polyhedron to .OFF file
bool 
Scene_textured_polyhedron_item::save(std::ostream& out) const
{
    out << *poly;
    return (bool) out;
}

QString 
Scene_textured_polyhedron_item::toolTip() const
{
    if(!poly)
        return QString();

    return QObject::tr("<p>Textured polyhedron <b>%1</b> (mode: %5, color: %6)</p>"
                       "<p>Number of vertices: %2<br />"
                       "Number of edges: %3<br />"
                       "Number of facets: %4</p>")
            .arg(this->name())
            .arg(poly->size_of_vertices())
            .arg(poly->size_of_halfedges()/2)
            .arg(poly->size_of_facets())
            .arg(this->renderingModeName())
            .arg(this->color().name());
}

// Points/Wireframe/Flat/Gouraud OpenGL drawing in a display list
void Scene_textured_polyhedron_item::draw(Viewer_interface* viewer) const {

    if(!are_buffers_filled)
        initialize_buffers(viewer);

    vaos[0]->bind();
    viewer->glActiveTexture(GL_TEXTURE0);
    viewer->glBindTexture(GL_TEXTURE_2D, textureId);
    attrib_buffers(viewer, PROGRAM_WITH_TEXTURE);
    program=getShaderProgram(PROGRAM_WITH_TEXTURE);
    program->bind();
    viewer->glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(nb_facets/4));
    //Clean-up
    program->release();
    vaos[0]->release();
}
void Scene_textured_polyhedron_item::draw_edges(Viewer_interface* viewer) const {
    if(!are_buffers_filled)
        initialize_buffers(viewer);

    vaos[1]->bind();
    viewer->glActiveTexture(GL_TEXTURE0);
    viewer->glBindTexture(GL_TEXTURE_2D, textureId);
    attrib_buffers(viewer, PROGRAM_WITH_TEXTURED_EDGES);

    program=getShaderProgram(PROGRAM_WITH_TEXTURED_EDGES);
    program->bind();
    viewer->glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(nb_lines/4));
    //Clean-up
    program->release();
    vaos[1]->release();
}

Textured_polyhedron* 
Scene_textured_polyhedron_item::textured_polyhedron()       { return poly; }
const Textured_polyhedron* 
Scene_textured_polyhedron_item::textured_polyhedron() const { return poly; }

bool
Scene_textured_polyhedron_item::isEmpty() const {
    return (poly == 0) || poly->empty();
}

Scene_textured_polyhedron_item::Bbox
Scene_textured_polyhedron_item::bbox() const {
    const Point& p = *(poly->points_begin());
    CGAL::Bbox_3 bbox(p.x(), p.y(), p.z(), p.x(), p.y(), p.z());
    for(Textured_polyhedron::Point_iterator it = poly->points_begin();
        it != poly->points_end();
        ++it) {
        bbox = bbox + it->bbox();
    }
    return Bbox(bbox.xmin(),bbox.ymin(),bbox.zmin(),
                bbox.xmax(),bbox.ymax(),bbox.zmax());
}
void
Scene_textured_polyhedron_item::invalidate_buffers()
{
    compute_normals_and_vertices();
    are_buffers_filled = false;}
void
Scene_textured_polyhedron_item::
contextual_changed()
{
    prev_shading = cur_shading;
    cur_shading = renderingMode();
    if(prev_shading != cur_shading)
    {
        invalidate_buffers();
    }
}
void
Scene_textured_polyhedron_item::selection_changed(bool p_is_selected)
{
    if(p_is_selected != is_selected)
    {
        is_selected = p_is_selected;
        initialize_buffers();
    }
    else
        is_selected = p_is_selected;
}
