#include "Scene_textured_polyhedron_item.h"
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <QApplication>
#include "Textured_polyhedron_type.h"
#include <QObject>

typedef EPIC_kernel::Point_3 Point;

struct Scene_textured_polyhedron_item_priv
{
  Scene_textured_polyhedron_item_priv(Scene_textured_polyhedron_item* parent)
    :poly(new Textured_polyhedron), textureId(-1)
  {
    item = parent;
    texture.GenerateCheckerBoard(2048,2048,128,0,0,0,250,250,255);
    nb_facets = 0;
    nb_lines = 0;
    nb_border = 0;
    positions_border.resize(0);
  }
  Scene_textured_polyhedron_item_priv(const Textured_polyhedron& p, Scene_textured_polyhedron_item* parent)
    : poly(new Textured_polyhedron(p)),textureId(-1),smooth_shading(true)

  {
    item = parent;
    texture.GenerateCheckerBoard(2048,2048,128,0,0,0,250,250,255);
    nb_facets = 0;
    nb_lines = 0;
  }
  Scene_textured_polyhedron_item_priv(Textured_polyhedron* const p,Scene_textured_polyhedron_item* parent)
    :poly(p),textureId(-1),smooth_shading(true)
  {
    item = parent;
    texture.GenerateCheckerBoard(2048,2048,128,0,0,0,250,250,255);
    nb_facets = 0;
    nb_lines = 0;

  }
  ~Scene_textured_polyhedron_item_priv()
  {
    delete poly;
  }
  void initializeBuffers(CGAL::Three::Viewer_interface *viewer) const;
  void compute_normals_and_vertices(void) const;

  enum VAOs {
      Facets=0,
      Edges,
      Border_edges,
      NbOfVaos
  };
  enum VBOs {
      Facets_Vertices=0,
      Facets_Normals,
      Facets_Texmap,
      Edges_Vertices,
      Edges_Texmap,
      Edges_border,
      NbOfVbos
  };

  Textured_polyhedron* poly;
  Texture texture;
  mutable std::vector<float> positions_lines;
  mutable std::vector<float> positions_facets;
  mutable std::vector<float> positions_border;
  mutable std::vector<float> normals;
  mutable std::vector<float> textures_map_facets;
  mutable std::vector<float> textures_map_lines;
  mutable std::size_t nb_facets;
  mutable std::size_t nb_lines;
  mutable std::size_t nb_border;

  mutable GLuint textureId;
  mutable QOpenGLShaderProgram* program;
  bool smooth_shading;

  Scene_textured_polyhedron_item* item;
};
void Scene_textured_polyhedron_item_priv::initializeBuffers(CGAL::Three::Viewer_interface *viewer = 0) const
{
    if(GLuint(-1) == textureId) {
        viewer->glGenTextures(1, &textureId);
    }
    //vao for the facets
    {
        program = item->getShaderProgram(Scene_textured_polyhedron_item::PROGRAM_WITH_TEXTURE, viewer);
        program->bind();
        item->vaos[Facets]->bind();
        item->buffers[Facets_Vertices].bind();
        item->buffers[Facets_Vertices].allocate(positions_facets.data(),
                            static_cast<int>(positions_facets.size()*sizeof(float)));
        program->enableAttributeArray("vertex");
        program->setAttributeBuffer("vertex",GL_FLOAT,0,4);
        item->buffers[Facets_Vertices].release();

        item->buffers[Facets_Normals].bind();
        item->buffers[Facets_Normals].allocate(normals.data(),
                            static_cast<int>(normals.size()*sizeof(float)));
        program->enableAttributeArray("normal");
        program->setAttributeBuffer("normal",GL_FLOAT,0,3);
        item->buffers[Facets_Normals].release();


        item->buffers[Facets_Texmap].bind();
        item->buffers[Facets_Texmap].allocate(textures_map_facets.data(),
                            static_cast<int>(textures_map_facets.size()*sizeof(float)));
        program->enableAttributeArray("v_texCoord");
        program->setAttributeBuffer("v_texCoord",GL_FLOAT,0,2);
        item->buffers[Facets_Texmap].release();
        item->vaos[Facets]->release();
        program->release();
    }

    //vao for the lines
    {
        program = item->getShaderProgram(Scene_textured_polyhedron_item::PROGRAM_WITH_TEXTURED_EDGES, viewer);
        program->bind();
        item->vaos[Edges]->bind();
        item->buffers[Edges_Vertices].bind();
        item->buffers[Edges_Vertices].allocate(positions_lines.data(),
                            static_cast<int>(positions_lines.size()*sizeof(float)));
        program->enableAttributeArray("vertex");
        program->setAttributeBuffer("vertex",GL_FLOAT,0,4);
        item->buffers[Edges_Vertices].release();


        item->buffers[Edges_Texmap].bind();
        item->buffers[Edges_Texmap].allocate(textures_map_lines.data(),
                            static_cast<int>(textures_map_lines.size()*sizeof(float)));
        program->enableAttributeArray("v_texCoord");
        program->setAttributeBuffer("v_texCoord",GL_FLOAT,0,2);
        item->buffers[Edges_Texmap].release();
        item->vaos[Edges]->release();
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
    //viewer->glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

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

    item->are_buffers_filled = true;
}

void
Scene_textured_polyhedron_item_priv::compute_normals_and_vertices(void) const
{
    QApplication::setOverrideCursor(Qt::WaitCursor);
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

    const CGAL::qglviewer::Vec offset = static_cast<CGAL::Three::Viewer_interface*>(CGAL::QGLViewer::QGLViewerPool().first())->offset();

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
            if (item->cur_shading == Flat || item->cur_shading == FlatPlusEdges)
            {

                Vector n = CGAL::Polygon_mesh_processing::
                  compute_face_normal(f, static_cast<Base&>(*poly));
                normals.push_back(n[0]);
                normals.push_back(n[1]);
                normals.push_back(n[2]);
            }

            // If Gouraud shading: 1 normal per vertex
            else if(item->cur_shading == Gouraud)
            {

                const Facet::Normal_3& n = he->vertex()->normal();
                normals.push_back(n[0]);
                normals.push_back(n[1]);
                normals.push_back(n[2]);
            }

            //position
            const Point& p = he->vertex()->point();
            positions_facets.push_back(p.x() + offset.x);
            positions_facets.push_back(p.y() + offset.y);
            positions_facets.push_back(p.z() + offset.z);
            positions_facets.push_back(1.0);

            const double u = he->u();
            const double v = he->v();
            textures_map_facets.push_back(u);
            textures_map_facets.push_back(v);
        }


    }
    //Lines
    typedef Kernel::Point_3                Point;
    typedef Base::Edge_iterator        Edge_iterator;

    Edge_iterator he;

    for(he = poly->edges_begin();
        he != poly->edges_end();
        he++)
    {

        const Point& a = he->vertex()->point();
        const Point& b = he->opposite()->vertex()->point();
        positions_lines.push_back(a.x() + offset.x);
        positions_lines.push_back(a.y() + offset.y);
        positions_lines.push_back(a.z() + offset.z);
        positions_lines.push_back(1.0);

        const double u = he->u();
        const double v = he->v();
        textures_map_lines.push_back(u);
        textures_map_lines.push_back(v);

        positions_lines.push_back(b.x()+ offset.x);
        positions_lines.push_back(b.y()+ offset.y);
        positions_lines.push_back(b.z()+ offset.z);
        positions_lines.push_back(1.0);

        const double ou = he->opposite()->u();
        const double ov = he->opposite()->v();
        textures_map_lines.push_back(ou);
        textures_map_lines.push_back(ov);

    }
    QApplication::restoreOverrideCursor();
}

Scene_textured_polyhedron_item::Scene_textured_polyhedron_item()
    : Scene_item(Scene_textured_polyhedron_item_priv::NbOfVbos,Scene_textured_polyhedron_item_priv::NbOfVaos)
{
    cur_shading=FlatPlusEdges;
    is_selected=false;
    d = new Scene_textured_polyhedron_item_priv(this);
    invalidateOpenGLBuffers();
}

Scene_textured_polyhedron_item::Scene_textured_polyhedron_item(Textured_polyhedron* const p)
    : Scene_item(Scene_textured_polyhedron_item_priv::NbOfVbos,Scene_textured_polyhedron_item_priv::NbOfVaos)
{
    cur_shading=FlatPlusEdges;
    is_selected=false;
    d = new Scene_textured_polyhedron_item_priv(p,this);
    invalidateOpenGLBuffers();
}

Scene_textured_polyhedron_item::Scene_textured_polyhedron_item(const Textured_polyhedron& p)
    : Scene_item(Scene_textured_polyhedron_item_priv::NbOfVbos,Scene_textured_polyhedron_item_priv::NbOfVaos)
{
    cur_shading=FlatPlusEdges;
    is_selected=false;
    d = new Scene_textured_polyhedron_item_priv(p,this);
    invalidateOpenGLBuffers();
}

Scene_textured_polyhedron_item::~Scene_textured_polyhedron_item()
{
    delete d;
}

Scene_textured_polyhedron_item*
Scene_textured_polyhedron_item::clone() const {
    return new Scene_textured_polyhedron_item(*d->poly);
}

// Load textured_polyhedron from .OFF file
bool
Scene_textured_polyhedron_item::load(std::istream& in)
{
    std::cout<<"LOAD"<<std::endl;
    in >> *d->poly;
    invalidateOpenGLBuffers();
    return in && !isEmpty();
}

// Write textured_polyhedron to .OFF file
bool
Scene_textured_polyhedron_item::save(std::ostream& out) const
{
    out << *d->poly;
    return (bool) out;
}

QString
Scene_textured_polyhedron_item::toolTip() const
{
    if(!d->poly)
        return QString();

    return QObject::tr("<p>Textured Polyhedron_3 <b>%1</b> (mode: %5, color: %6)</p>"
                       "<p>Number of vertices: %2<br />"
                       "Number of edges: %3<br />"
                       "Number of facets: %4</p>")
            .arg(this->name())
            .arg(d->poly->size_of_vertices())
            .arg(d->poly->size_of_halfedges()/2)
            .arg(d->poly->size_of_facets())
            .arg(this->renderingModeName())
            .arg(this->color().name());
}

// Points/Wireframe/Flat/Gouraud OpenGL drawing in a display list
void Scene_textured_polyhedron_item::draw(CGAL::Three::Viewer_interface* viewer) const {

    if(!are_buffers_filled)
    {
        d->compute_normals_and_vertices();
        d->initializeBuffers(viewer);
    }

    vaos[Scene_textured_polyhedron_item_priv::Facets]->bind();
    viewer->glActiveTexture(GL_TEXTURE0);
    viewer->glBindTexture(GL_TEXTURE_2D, d->textureId);
    attribBuffers(viewer, PROGRAM_WITH_TEXTURE);
    d->program=getShaderProgram(PROGRAM_WITH_TEXTURE);
    d->program->bind();
    viewer->glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(d->nb_facets/4));
    //Clean-up
    d->program->release();
    vaos[Scene_textured_polyhedron_item_priv::Facets]->release();
}
void Scene_textured_polyhedron_item::drawEdges(CGAL::Three::Viewer_interface* viewer) const {
    if(!are_buffers_filled)
        d->initializeBuffers(viewer);

    vaos[Scene_textured_polyhedron_item_priv::Edges]->bind();
    viewer->glActiveTexture(GL_TEXTURE0);
    viewer->glBindTexture(GL_TEXTURE_2D, d->textureId);
    attribBuffers(viewer, PROGRAM_WITH_TEXTURED_EDGES);

    d->program=getShaderProgram(PROGRAM_WITH_TEXTURED_EDGES);
    d->program->bind();
    viewer->glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(d->nb_lines/4));

    vaos[Scene_textured_polyhedron_item_priv::Edges]->release();
    d->program->release();

    vaos[Scene_textured_polyhedron_item_priv::Border_edges]->bind();
    if(!viewer->isOpenGL_4_3())
    {
      attribBuffers(viewer, PROGRAM_NO_SELECTION);
      d->program=getShaderProgram(PROGRAM_NO_SELECTION);
      d->program->bind();
    }
    else
    {
      attribBuffers(viewer, PROGRAM_SOLID_WIREFRAME);
      d->program=getShaderProgram(PROGRAM_SOLID_WIREFRAME);
      d->program->bind();
      QVector2D vp(viewer->width(), viewer->height());
      d->program->setUniformValue("viewport", vp);
      d->program->setUniformValue("near",(GLfloat)viewer->camera()->zNear());
      d->program->setUniformValue("far",(GLfloat)viewer->camera()->zFar());
      d->program->setUniformValue("width", 4.0f);
    }
    d->program->setAttributeValue("colors", QColor(Qt::blue));
    viewer->glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(d->nb_border/3));
    //Clean-up
    d->program->release();
    vaos[Scene_textured_polyhedron_item_priv::Border_edges]->release();
}

Textured_polyhedron*
Scene_textured_polyhedron_item::textured_polyhedron()       { return d->poly; }
const Textured_polyhedron*
Scene_textured_polyhedron_item::textured_polyhedron() const { return d->poly; }

bool
Scene_textured_polyhedron_item::isEmpty() const {
    return (d->poly == 0) || d->poly->empty();
}

void
Scene_textured_polyhedron_item::compute_bbox() const {
    const Point& p = *(d->poly->points_begin());
    CGAL::Bbox_3 bbox(p.x(), p.y(), p.z(), p.x(), p.y(), p.z());
    for(Textured_polyhedron::Point_iterator it = d->poly->points_begin();
        it != d->poly->points_end();
        ++it) {
        bbox = bbox + it->bbox();
    }
    _bbox = Bbox(bbox.xmin(),bbox.ymin(),bbox.zmin(),
                bbox.xmax(),bbox.ymax(),bbox.zmax());
}
void
Scene_textured_polyhedron_item::invalidateOpenGLBuffers()
{
    are_buffers_filled = false;
    compute_bbox();
}
void
Scene_textured_polyhedron_item::selection_changed(bool p_is_selected)
{
    if(p_is_selected != is_selected)
    {
        is_selected = p_is_selected;
        initializeBuffers();
    //to be replaced by a functor in the d-pointer when the merging is done
        if(p_is_selected)
          Q_EMIT selectionChanged();
    }
    else
        is_selected = p_is_selected;
}
void Scene_textured_polyhedron_item::add_border_edges(std::vector<float> border_edges,
                                                      bool is_opengl_4_3)
{
  d->positions_border = border_edges;
  d->nb_border = border_edges.size();
  d->program=is_opengl_4_3
      ? getShaderProgram(PROGRAM_SOLID_WIREFRAME)
      : getShaderProgram(PROGRAM_NO_SELECTION);
  d->program->bind();
  vaos[Scene_textured_polyhedron_item_priv::Border_edges]->bind();
  buffers[Scene_textured_polyhedron_item_priv::Edges_border].bind();
  buffers[Scene_textured_polyhedron_item_priv::Edges_Vertices].allocate(d->positions_border.data(),
                      static_cast<int>(d->positions_border.size()*sizeof(float)));
  d->program->enableAttributeArray("vertex");
  d->program->setAttributeBuffer("vertex",GL_FLOAT,0,3);
  buffers[Scene_textured_polyhedron_item_priv::Edges_Vertices].release();
  vaos[Scene_textured_polyhedron_item_priv::Border_edges]->release();

  d->program->release();
  d->positions_border.resize(0);
  std::vector<float>(d->positions_border).swap(d->positions_border);
  itemChanged();

}
