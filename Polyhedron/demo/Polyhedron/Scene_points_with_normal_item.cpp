#include "Scene_points_with_normal_item.h"
#include "Polyhedron_type.h"
#include <CGAL/internal/Operations_on_polyhedra/compute_normal.h>

#include <CGAL/IO/read_off_points.h>
#include <CGAL/IO/write_off_points.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/IO/write_xyz_points.h>
#include <CGAL/Timer.h>
#include <CGAL/Memory_sizer.h>

#include "Viewer_interface.h"
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_3.h>

#include <QObject>
#include <QMenu>

#include <set>
#include <stack>
#include <algorithm>
#include <boost/array.hpp>


struct light_info
{
    //position
    GLfloat position[4];

    //ambient
    GLfloat ambient[4];

    //diffuse
    GLfloat diffuse[4];

    //specular
    GLfloat specular[4];
};

Scene_points_with_normal_item::Scene_points_with_normal_item()
    : Scene_item(6,3),
      m_points(new Point_set),
      m_has_normals(false)
{
    setRenderingMode(Points);
    is_selected = true;
    qFunc.initializeOpenGLFunctions();
    qFunc.glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    qFunc.glGenTextures(1, &textureId);
}

// Copy constructor
Scene_points_with_normal_item::Scene_points_with_normal_item(const Scene_points_with_normal_item& toCopy)
    : Scene_item(6,3), // do not call superclass' copy constructor
      m_points(new Point_set(*toCopy.m_points)),
      m_has_normals(toCopy.m_has_normals)
{
    if (m_has_normals)
    {
        setRenderingMode(PointsPlusNormals);
        is_selected = true;
        qFunc.initializeOpenGLFunctions();
        qFunc.glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
        qFunc.glGenTextures(1, &textureId);
    }
    else
    {
        setRenderingMode(Points);
        is_selected = true;
        qFunc.initializeOpenGLFunctions();
        qFunc.glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
        qFunc.glGenTextures(1, &textureId);
    }
}

// Converts polyhedron to point set
Scene_points_with_normal_item::Scene_points_with_normal_item(const Polyhedron& input_mesh)
    : Scene_item(6,3),
      m_points(new Point_set),
      m_has_normals(true)
{
    // Converts Polyhedron vertices to point set.
    // Computes vertices normal from connectivity.
    Polyhedron::Vertex_const_iterator v;
    for (v = input_mesh.vertices_begin(); v != input_mesh.vertices_end(); v++)
    {
        const Kernel::Point_3& p = v->point();
        Kernel::Vector_3 n = compute_vertex_normal<Polyhedron::Vertex,Kernel>(*v);
        m_points->push_back(UI_point(p,n));
    }

    setRenderingMode(PointsPlusNormals);
    is_selected = true;
    qFunc.initializeOpenGLFunctions();
    //Generates an integer which will be used as ID for each buffer
    qFunc.glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    qFunc.glGenTextures(1, &textureId);
}

Scene_points_with_normal_item::~Scene_points_with_normal_item()
{
    Q_ASSERT(m_points != NULL);
    delete m_points; m_points = NULL;
}



void Scene_points_with_normal_item::initialize_buffers(Viewer_interface *viewer) const
{
    //vao for the edges
    {
        program = getShaderProgram(PROGRAM_WITHOUT_LIGHT, viewer);
        program->bind();

        vaos[0]->bind();
        buffers[0].bind();
        buffers[0].allocate(positions_lines.data(), positions_lines.size()*sizeof(double));
        program->enableAttributeArray("vertex");
        program->setAttributeBuffer("vertex",GL_DOUBLE,0,3);
        buffers[0].release();



        buffers[1].bind();
        buffers[1].allocate(color_lines.data(), color_lines.size()*sizeof(double));
        program->enableAttributeArray("colors");
        program->setAttributeBuffer("colors",GL_DOUBLE,0,3);
        buffers[1].release();

        vaos[0]->release();
        program->release();
    }
    //vao for the points
    {
        program = getShaderProgram(PROGRAM_WITHOUT_LIGHT, viewer);
        program->bind();

        vaos[1]->bind();
        buffers[2].bind();
        buffers[2].allocate(positions_points.data(), positions_points.size()*sizeof(double));
        program->enableAttributeArray("vertex");
        program->setAttributeBuffer("vertex",GL_DOUBLE,0,3);
        buffers[2].release();



        buffers[3].bind();
        buffers[3].allocate(color_points.data(), color_points.size()*sizeof(double));
        program->enableAttributeArray("colors");
        program->setAttributeBuffer("colors",GL_DOUBLE,0,3);
        buffers[3].release();

        vaos[1]->release();
        program->release();
    }
    //vao for the selected points
    {
        program = getShaderProgram(PROGRAM_WITHOUT_LIGHT, viewer);
        program->bind();

        vaos[2]->bind();
        buffers[4].bind();
        buffers[4].allocate(positions_selected_points.data(), positions_selected_points.size()*sizeof(double));
        program->enableAttributeArray("vertex");
        program->setAttributeBuffer("vertex",GL_DOUBLE,0,3);
        buffers[4].release();



        buffers[5].bind();
        buffers[5].allocate(color_selected_points.data(), color_selected_points.size()*sizeof(double));
        program->enableAttributeArray("colors");
        program->setAttributeBuffer("colors",GL_DOUBLE,0,3);
        buffers[5].release();

        vaos[2]->release();
        program->release();
    }
    are_buffers_filled = true;


}
void Scene_points_with_normal_item::compute_normals_and_vertices(void)
{
    positions_points.resize(0);
    positions_lines.resize(0);
    positions_splats.resize(0);
    positions_selected_points.resize(0);
    color_selected_points.resize(0);
    color_points.resize(0);
    color_lines.resize(0);
    normals.resize(0);
    tex_coords.resize(0);

    //The points
    {
        // The *non-selected* points
        if (m_points->nb_selected_points()< m_points->size())
        {

            for (Point_set_3<Kernel>::const_iterator it = m_points->begin(); it != m_points->end(); it++)
            {
                const UI_point& p = *it;
                if ( ! p.is_selected() )
                {
                    positions_points.push_back(p.x());
                    positions_points.push_back(p.y());
                    positions_points.push_back(p.z());

                    if(is_selected)
                    {
                        color_points.push_back(this->color().lighter(120).redF());
                        color_points.push_back(this->color().lighter(120).greenF());
                        color_points.push_back(this->color().lighter(120).blueF());
                    }
                    else
                    {
                        color_points.push_back(this->color().redF());
                        color_points.push_back(this->color().greenF());
                        color_points.push_back(this->color().blueF());
                    }

                }
            }

        }

        // Draw *selected* points
        if (m_points->nb_selected_points() > 0)
        {
            ::glPointSize(4.f);    // selected => bigger
            for (Point_set_3<Kernel>::const_iterator it = m_points->begin(); it != m_points->end(); it++)
            {
                const UI_point& p = *it;
                if (p.is_selected())
                {
                    positions_selected_points.push_back(p.x());
                    positions_selected_points.push_back(p.y());
                    positions_selected_points.push_back(p.z());

                    color_selected_points.push_back(1.0);
                    color_selected_points.push_back(0.0);
                    color_selected_points.push_back(0.0);

                }
            }

        }
    }

    //The lines
    {
        // Stock normals
        Kernel::Sphere_3 region_of_interest = m_points->region_of_interest();
        float normal_length = (float)std::sqrt(region_of_interest.squared_radius() / 1000.);

        // Stock normals of *non-selected* points
        if (m_points->nb_selected_points() < m_points->size())
        {
            // Stock normals
            for (Point_set_3<Kernel>::const_iterator it = m_points->begin(); it != m_points->end(); it++)
            {
                const UI_point& p = *it;
                const Point_set_3<Kernel>::Vector& n = p.normal();
                if (!p.is_selected())
                {
                    Point_set_3<Kernel>::Point q = p + normal_length * n;
                    positions_lines.push_back(p.x());
                    positions_lines.push_back(p.y());
                    positions_lines.push_back(p.z());

                    positions_lines.push_back(q.x());
                    positions_lines.push_back(q.y());
                    positions_lines.push_back(q.z());

                    if(is_selected)
                    {
                        color_lines.push_back(this->color().lighter(120).redF());
                        color_lines.push_back(this->color().lighter(120).greenF());
                        color_lines.push_back(this->color().lighter(120).blueF());

                        color_lines.push_back(this->color().lighter(120).redF());
                        color_lines.push_back(this->color().lighter(120).greenF());
                        color_lines.push_back(this->color().lighter(120).blueF());
                    }
                    else
                    {
                        color_lines.push_back(this->color().redF());
                        color_lines.push_back(this->color().greenF());
                        color_lines.push_back(this->color().blueF());

                        color_lines.push_back(this->color().redF());
                        color_lines.push_back(this->color().greenF());
                        color_lines.push_back(this->color().blueF());
                    }
                }
            }
        }

        // Stock normals of *selected* points
        if (m_points->nb_selected_points() > 0)
        {
            for (Point_set_3<Kernel>::const_iterator it = m_points->begin(); it != m_points->end(); it++)
            {
                const UI_point& p = *it;
                const Point_set_3<Kernel>::Vector& n = p.normal();
                if (p.is_selected())
                {
                    Point_set_3<Kernel>::Point q = p + normal_length * n;
                    positions_lines.push_back(p.x());
                    positions_lines.push_back(p.y());
                    positions_lines.push_back(p.z());

                    positions_lines.push_back(q.x());
                    positions_lines.push_back(q.y());
                    positions_lines.push_back(q.z());

                    color_lines.push_back(1.0);
                    color_lines.push_back(0.0);
                    color_lines.push_back(0.0);

                    color_lines.push_back(1.0);
                    color_lines.push_back(0.0);
                    color_lines.push_back(0.0);
                }
            }
        }
    }

    //The splats
    {
        // TODO add support for selection
        texture[0] = this->color().redF();
        texture[1] = this->color().greenF();
        texture[2] = this->color().blueF();

        // Draw splats
        bool points_have_normals = (m_points->begin() != m_points->end() &&
                m_points->begin()->normal() != CGAL::NULL_VECTOR);
        bool points_have_radii =   (m_points->begin() != m_points->end() &&
                m_points->begin()->radius() != 0);
        if(points_have_normals && points_have_radii)
        {
            for (Point_set_3<Kernel>::const_iterator it = m_points->begin(); it != m_points->end(); it++)
            {
                const UI_point& p = *it;
                normals.push_back(p.normal().x());
                normals.push_back(p.normal().y());
                normals.push_back(p.normal().z());
                tex_coords.push_back(p.radius());
                tex_coords.push_back(0);

                positions_splats.push_back(p.x());
                positions_splats.push_back(p.y());
                positions_splats.push_back(p.z());
            }
        }
    }
}
/*
void Scene_points_with_normal_item::uniform_attrib(Viewer_interface* viewer, int mode) const
{
    GLfloat mvp_mat[16];
    light_info light;
    GLint is_both_sides = 0;
    GLfloat mv_mat[16];

    GLdouble d_mat[16];
    viewer->camera()->getModelViewProjectionMatrix(d_mat);

    for (int i=0; i<16; ++i){
        mvp_mat[i] = GLfloat(d_mat[i]);
    }
    viewer->camera()->getModelViewMatrix(d_mat);
    for (int i=0; i<16; ++i)
        mv_mat[i] = GLfloat(d_mat[i]);

    qFunc.glGetIntegerv(GL_LIGHT_MODEL_TWO_SIDE, &is_both_sides);

    //Gets lighting info :

    //position
    glGetLightfv(GL_LIGHT0, GL_POSITION, light.position);

    //ambient
    glGetLightfv(GL_LIGHT0, GL_AMBIENT, light.ambient);


    //specular
    glGetLightfv(GL_LIGHT0, GL_SPECULAR, light.specular);

    //diffuse
    glGetLightfv(GL_LIGHT0, GL_DIFFUSE, light.diffuse);

    if(mode ==0)
    {
        qFunc.glUseProgram(rendering_program_lines);
        qFunc.glUniformMatrix4fv(location[0], 1, GL_FALSE, mvp_mat);
    }
    else if(mode ==1)
    {
        qFunc.glUseProgram(rendering_program_points);
        qFunc.glUniformMatrix4fv(location[1], 1, GL_FALSE, mvp_mat);
    }

    else if(mode ==2)
    {
        qFunc.glUseProgram(rendering_program_splats);
        qFunc.glUniformMatrix4fv(location[2], 1, GL_FALSE, mvp_mat);
        qFunc.glUniformMatrix4fv(location[3], 1, GL_FALSE, mv_mat);
        qFunc.glUniform3fv(location[4], 1, light.position);
        qFunc.glUniform3fv(location[5], 1, light.diffuse);
        qFunc.glUniform3fv(location[6], 1, light.specular);
        qFunc.glUniform3fv(location[7], 1, light.ambient);
        qFunc.glUniform1i(location[8], is_both_sides);
        qFunc.glUniform1i(sampler_location, 0);
    }

}

*/
// Duplicates scene item
Scene_points_with_normal_item*
Scene_points_with_normal_item::clone() const
{
    return new Scene_points_with_normal_item(*this);
}

// Is selection empty?
bool Scene_points_with_normal_item::isSelectionEmpty() const
{
    return (m_points->nb_selected_points() == 0);
}

// Delete selection
void Scene_points_with_normal_item::deleteSelection()
{
    CGAL::Timer task_timer; task_timer.start();
    std::cerr << "Delete " << m_points->nb_selected_points() << " points...";

    // Delete selected points
    m_points->delete_selection();

    std::size_t memory = CGAL::Memory_sizer().virtual_size();
    std::cerr << "done: " << task_timer.time() << " seconds, "
              << (memory>>20) << " Mb allocated"
              << std::endl;
  Q_EMIT itemChanged();
}

// Reset selection mark
void Scene_points_with_normal_item::resetSelection()
{
    // Un-select all points
    m_points->select(m_points->begin(), m_points->end(), false);
  Q_EMIT itemChanged();
}
//Select duplicated points
void Scene_points_with_normal_item::selectDuplicates()
{
    std::set<Kernel::Point_3> unique_points;
    for (Point_set::Point_iterator ptit=m_points->begin(); ptit!=m_points->end();++ptit )
        if ( !unique_points.insert(*ptit).second )
            m_points->select(&(*ptit));
  Q_EMIT itemChanged();
}

// Loads point set from .OFF file
bool Scene_points_with_normal_item::read_off_point_set(std::istream& stream)
{
    Q_ASSERT(m_points != NULL);

    m_points->clear();
    bool ok = stream &&
            CGAL::read_off_points_and_normals(stream,
                                              std::back_inserter(*m_points),
                                              CGAL::make_normal_of_point_with_normal_pmap(Point_set::value_type())) &&
            !isEmpty();

    return ok;
}

// Write point set to .OFF file
bool Scene_points_with_normal_item::write_off_point_set(std::ostream& stream) const
{
    Q_ASSERT(m_points != NULL);

    return stream &&
            CGAL::write_off_points_and_normals(stream,
                                               m_points->begin(), m_points->end(),
                                               CGAL::make_normal_of_point_with_normal_pmap(Point_set::value_type()));
}

// Loads point set from .XYZ file
bool Scene_points_with_normal_item::read_xyz_point_set(std::istream& stream)
{
    Q_ASSERT(m_points != NULL);

    m_points->clear();
    bool ok = stream &&
            CGAL::read_xyz_points_and_normals(stream,
                                              std::back_inserter(*m_points),
                                              CGAL::make_normal_of_point_with_normal_pmap(Point_set::value_type())) &&
            !isEmpty();

    if (ok)
    {
        for (Point_set::iterator it=m_points->begin(),
             end=m_points->end();it!=end; ++it)
        {
            if (it->normal() != CGAL::NULL_VECTOR)
            {
                m_has_normals=true;
                setRenderingMode(PointsPlusNormals);
                break;
            }
        }
    }
    return ok;
}

// Write point set to .XYZ file
bool Scene_points_with_normal_item::write_xyz_point_set(std::ostream& stream) const
{
    Q_ASSERT(m_points != NULL);

    return stream &&
            CGAL::write_xyz_points_and_normals(stream,
                                               m_points->begin(), m_points->end(),
                                               CGAL::make_normal_of_point_with_normal_pmap(Point_set::value_type()));
}

QString
Scene_points_with_normal_item::toolTip() const
{
    Q_ASSERT(m_points != NULL);

    return QObject::tr("<p><b>%1</b> (color: %4)<br />"
                       "<i>Point set</i></p>"
                       "<p>Number of points: %2</p>")
            .arg(name())
            .arg(m_points->size())
            .arg(color().name());
}

bool Scene_points_with_normal_item::supportsRenderingMode(RenderingMode m) const 
{
    return m==Points ||
            ( has_normals() &&
              ( m==PointsPlusNormals || m==Splatting ) );
}

void Scene_points_with_normal_item::draw_splats(Viewer_interface* viewer) const
{
    //Needs to be re-thinked because the GlSplat Renderer is deprecated and is a big part of the scene class.

   // TODO add support for selection
/*   ::glBegin(GL_POINTS);
   for ( Point_set_3<Kernel>::const_iterator it = m_points->begin(); it != m_points->end(); it++)
   {
     const UI_point& p = *it;
     ::glNormal3dv(&p.normal().x());
     ::glMultiTexCoord1d(GL_TEXTURE2, p.radius());
     ::glVertex3dv(&p.x());

   }
   ::glEnd();


*/

}

void Scene_points_with_normal_item::draw_edges(Viewer_interface* viewer) const
{
    if(!are_buffers_filled)
        initialize_buffers(viewer);
    vaos[0]->bind();
    program=getShaderProgram(PROGRAM_WITHOUT_LIGHT);
    attrib_buffers(viewer,PROGRAM_WITHOUT_LIGHT);
    program->bind();
    qFunc.glDrawArrays(GL_LINES, 0, positions_lines.size()/3);
    vaos[0]->release();
    program->release();
}
void Scene_points_with_normal_item::draw_points(Viewer_interface* viewer) const
{
    if(!are_buffers_filled)
        initialize_buffers(viewer);

    vaos[1]->bind();
    program=getShaderProgram(PROGRAM_WITHOUT_LIGHT);
    attrib_buffers(viewer,PROGRAM_WITHOUT_LIGHT);
    program->bind();
    qFunc.glDrawArrays(GL_POINTS, 0, positions_points.size()/3);
    vaos[1]->release();
    program->release();

    GLfloat point_size;
    qFunc.glGetFloatv(GL_POINT_SIZE, &point_size);
    qFunc.glPointSize(4.f);

    vaos[2]->bind();
    program=getShaderProgram(PROGRAM_WITHOUT_LIGHT);
    attrib_buffers(viewer,PROGRAM_WITHOUT_LIGHT);
    program->bind();
    qFunc.glDrawArrays(GL_POINTS, 0, positions_selected_points.size()/3);
    vaos[2]->release();
    program->release();
    qFunc.glPointSize(point_size);
}
// Gets wrapped point set
Point_set* Scene_points_with_normal_item::point_set()
{
    Q_ASSERT(m_points != NULL);
    return m_points;
}
const Point_set* Scene_points_with_normal_item::point_set() const
{
    Q_ASSERT(m_points != NULL);
    return m_points;
}

bool
Scene_points_with_normal_item::isEmpty() const
{
    Q_ASSERT(m_points != NULL);
    return m_points->empty();
}

Scene_points_with_normal_item::Bbox
Scene_points_with_normal_item::bbox() const
{
    Q_ASSERT(m_points != NULL);

    Kernel::Iso_cuboid_3 bbox = m_points->bounding_box();
    return Bbox(bbox.xmin(),bbox.ymin(),bbox.zmin(),
                bbox.xmax(),bbox.ymax(),bbox.zmax());
}

void Scene_points_with_normal_item::computes_local_spacing(int k)
{
    typedef Kernel Geom_traits;
    typedef CGAL::Search_traits_3<Geom_traits> TreeTraits;
    typedef CGAL::Orthogonal_k_neighbor_search<TreeTraits> Neighbor_search;
    typedef Neighbor_search::Tree Tree;

    Point_set::iterator end(m_points->end());

    // build kdtree
    Tree tree(m_points->begin(), end);

    // Compute the radius of each point = (distance max to k nearest neighbors)/2.
    {
        int i=0;
        for (Point_set::iterator it=m_points->begin(); it!=end; ++it, ++i)
        {
            Neighbor_search search(tree, *it, k+1);
            double maxdist2 = (--search.end())->second; // squared distance to furthest neighbor
            it->radius() = sqrt(maxdist2)/2.;
        }
    }

    m_points->set_radii_uptodate(true);
}

QMenu* Scene_points_with_normal_item::contextMenu()
{
    const char* prop_name = "Menu modified by Scene_points_with_normal_item.";

    QMenu* menu = Scene_item::contextMenu();

    // Use dynamic properties:
    // http://doc.trolltech.com/lastest/qobject.html#property
    bool menuChanged = menu->property(prop_name).toBool();

    if(!menuChanged) {
        actionDeleteSelection = menu->addAction(tr("Delete Selection"));
        actionDeleteSelection->setObjectName("actionDeleteSelection");
        connect(actionDeleteSelection, SIGNAL(triggered()),this, SLOT(deleteSelection()));

        actionResetSelection = menu->addAction(tr("Reset Selection"));
        actionResetSelection->setObjectName("actionResetSelection");
        connect(actionResetSelection, SIGNAL(triggered()),this, SLOT(resetSelection()));

        actionSelectDuplicatedPoints = menu->addAction(tr("Select duplicated points"));
        actionSelectDuplicatedPoints->setObjectName("actionSelectDuplicatedPoints");
        connect(actionSelectDuplicatedPoints, SIGNAL(triggered()),this, SLOT(selectDuplicates()));

        menu->setProperty(prop_name, true);
    }

    if (isSelectionEmpty())
    {
        actionDeleteSelection->setDisabled(true);
        actionResetSelection->setDisabled(true);
    }
    else
    {
        actionDeleteSelection->setDisabled(false);
        actionResetSelection->setDisabled(false);
    }

    return menu;
}

void Scene_points_with_normal_item::setRenderingMode(RenderingMode m)
{
    Scene_item::setRenderingMode(m);
    if (rendering_mode==Splatting && (!m_points->are_radii_uptodate()))
    {
        computes_local_spacing(6); // default value = small
    }
}

bool Scene_points_with_normal_item::has_normals() const { return m_has_normals; }

void Scene_points_with_normal_item::set_has_normals(bool b) {
    if (b!=m_has_normals){
        m_has_normals=b;
        //reset the context menu
        delete this->defaultContextMenu;
        this->defaultContextMenu = 0;
    }
}

void Scene_points_with_normal_item::changed()
{

    compute_normals_and_vertices();
    are_buffers_filled = false;
}
void Scene_points_with_normal_item::selection_changed(bool p_is_selected)
{
    if(p_is_selected != is_selected)
    {
        is_selected = p_is_selected;
        changed();
    }
}

#include "Scene_points_with_normal_item.moc"
