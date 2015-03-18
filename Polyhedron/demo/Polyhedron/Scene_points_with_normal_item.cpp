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


Scene_points_with_normal_item::Scene_points_with_normal_item()
    : Scene_item_with_display_list(),
      m_points(new Point_set),
      positions_lines(0),
      color_lines(0),
      color_points(0),
      positions_points(0),
      m_has_normals(false)
{
    setRenderingMode(Points);
    is_selected = true;
    glGenVertexArrays(1, vao);
    //Generates an integer which will be used as ID for each buffer
    glGenBuffers(7, buffer);
    compile_shaders();
}

// Copy constructor
Scene_points_with_normal_item::Scene_points_with_normal_item(const Scene_points_with_normal_item& toCopy)
    : Scene_item_with_display_list(), // do not call superclass' copy constructor
      m_points(new Point_set(*toCopy.m_points)),
      positions_lines(0),
      color_lines(0),
      color_points(0),
      positions_points(0),
      m_has_normals(toCopy.m_has_normals)
{
    if (m_has_normals)
    {
        setRenderingMode(PointsPlusNormals);
        is_selected = true;
        glGenVertexArrays(1, vao);
        //Generates an integer which will be used as ID for each buffer
        glGenBuffers(7, buffer);
        compile_shaders();
    }
    else
    {
        setRenderingMode(Points);
        is_selected = true;
        glGenVertexArrays(1, vao);
        //Generates an integer which will be used as ID for each buffer
        glGenBuffers(7, buffer);
        compile_shaders();
    }
}

// Converts polyhedron to point set
Scene_points_with_normal_item::Scene_points_with_normal_item(const Polyhedron& input_mesh)
    : Scene_item_with_display_list(),
      m_points(new Point_set),
      positions_lines(0),
      color_lines(0),
      color_points(0),
      positions_points(0),
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
    glGenVertexArrays(1, vao);
    //Generates an integer which will be used as ID for each buffer
    glGenBuffers(7, buffer);
    compile_shaders();
}

Scene_points_with_normal_item::~Scene_points_with_normal_item()
{
    Q_ASSERT(m_points != NULL);
    glDeleteBuffers(7, buffer);
    glDeleteVertexArrays(1, vao);
    glDeleteProgram(rendering_program_lines);
    glDeleteProgram(rendering_program_points);
    delete m_points; m_points = NULL;
}



void Scene_points_with_normal_item::initialize_buffers()
{
    glBindVertexArray(vao[0]);

    glBindBuffer(GL_ARRAY_BUFFER, buffer[0]);
    glBufferData(GL_ARRAY_BUFFER,
                 (positions_lines.size())*sizeof(double),
                 positions_lines.data(),
                 GL_STATIC_DRAW);
    glVertexAttribPointer(0, //number of the buffer
                          3, //number of floats to be taken
                          GL_DOUBLE, // type of data
                          GL_FALSE, //not normalized
                          0, //compact data (not in a struct)
                          NULL //no offset (seperated in several buffers)
                          );
    glEnableVertexAttribArray(0);


    glBindBuffer(GL_ARRAY_BUFFER, buffer[1]);
    glBufferData(GL_ARRAY_BUFFER,
                 (color_lines.size())*sizeof(double),
                 color_lines.data(), GL_STATIC_DRAW);
    glVertexAttribPointer(1,
                          3,
                          GL_DOUBLE,
                          GL_FALSE,
                          0,
                          NULL
                          );
    glEnableVertexAttribArray(1);

    glBindBuffer(GL_ARRAY_BUFFER, buffer[2]);
    glBufferData(GL_ARRAY_BUFFER,
                 (positions_points.size())*sizeof(double),
                 positions_points.data(), GL_STATIC_DRAW);
    glVertexAttribPointer(2,
                          3,
                          GL_DOUBLE,
                          GL_FALSE,
                          0,
                          NULL
                          );
    glEnableVertexAttribArray(2);

    glBindBuffer(GL_ARRAY_BUFFER, buffer[3]);
    glBufferData(GL_ARRAY_BUFFER,
                 (color_points.size())*sizeof(double),
                 color_points.data(), GL_STATIC_DRAW);
    glVertexAttribPointer(3,
                          3,
                          GL_DOUBLE,
                          GL_FALSE,
                          0,
                          NULL
                          );
    glEnableVertexAttribArray(3);

    glBindVertexArray(vao[1]);

    glBindBuffer(GL_ARRAY_BUFFER, buffer[4]);
    glBufferData(GL_ARRAY_BUFFER,
                 (positions_selected_points.size())*sizeof(double),
                 positions_selected_points.data(),
                 GL_STATIC_DRAW);
    glVertexAttribPointer(2, //number of the buffer
                          3, //number of floats to be taken
                          GL_DOUBLE, // type of data
                          GL_FALSE, //not normalized
                          0, //compact data (not in a struct)
                          NULL //no offset (seperated in several buffers)
                          );
    glEnableVertexAttribArray(2);


    glBindBuffer(GL_ARRAY_BUFFER, buffer[5]);
    glBufferData(GL_ARRAY_BUFFER,
                 (color_selected_points.size())*sizeof(double),
                 color_selected_points.data(), GL_STATIC_DRAW);
    glVertexAttribPointer(3,
                          3,
                          GL_DOUBLE,
                          GL_FALSE,
                          0,
                          NULL
                          );
    glEnableVertexAttribArray(3);


    // Clean-up
    glBindVertexArray(0);
}
void Scene_points_with_normal_item::compile_shaders(void)
{
    //fill the vertex shader
    //For the edges
    static const GLchar* vertex_shader_source_lines[] =
    {
        "#version 300 es \n"
        " \n"
        "layout (location = 0) in vec3 positions; \n"
        "layout (location = 1) in vec3 color_lines; \n"

        "uniform mat4 mvp_matrix; \n"

        "out highp vec3 fColors; \n"
        "vec4 positions_lines = vec4(positions, 1.0); \n"
        " \n"

        "void main(void) \n"
        "{ \n"
        "   fColors = color_lines; \n"
        "   gl_Position = mvp_matrix * positions_lines; \n"
        "} \n"
    };
    //fill the fragment shader
    static const GLchar* fragment_shader_source[]=
    {
        "#version 300 es \n"
        " \n"
        "precision mediump float; \n"
        "in vec3 fColors; \n"

        "out vec3 color; \n"
        " \n"
        "void main(void) \n"
        "{ \n"
        " color = fColors; \n"
        "} \n"
    };

    GLuint vertex_shader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertex_shader, 1, vertex_shader_source_lines, NULL);
    glCompileShader(vertex_shader);
    GLuint fragment_shader =	glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragment_shader, 1, fragment_shader_source, NULL);
    glCompileShader(fragment_shader);

    //creates the program, attaches and links the shaders
    GLuint program= glCreateProgram();
    glAttachShader(program, vertex_shader);
    glAttachShader(program, fragment_shader);
    glLinkProgram(program);

    //Clean-up
    glDeleteShader(vertex_shader);

    rendering_program_lines = program;


    //For the points
    static const GLchar* vertex_shader_source_points[] =
    {
        "#version 300 es \n"
        " \n"
        "layout (location = 2) in vec3 positions; \n"
        "layout (location = 3) in vec3 color; \n"
        "uniform mat4 mvp_matrix; \n"

        "out highp vec3 fColors; \n"
        " \n"

        "void main(void) \n"
        "{ \n"
        "   fColors = color; \n"
        "   gl_Position = mvp_matrix * vec4(positions, 1.0); \n"
        "} \n"
    };

    vertex_shader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertex_shader, 1, vertex_shader_source_points, NULL);
    glCompileShader(vertex_shader);

    glShaderSource(fragment_shader, 1, fragment_shader_source, NULL);
    glCompileShader(fragment_shader);

    program = glCreateProgram();
    glAttachShader(program, vertex_shader);
    glAttachShader(program, fragment_shader);
    glLinkProgram(program);
    //Clean-up
    glDeleteShader(vertex_shader);
    glDeleteShader(fragment_shader);
    rendering_program_points = program;
}
void Scene_points_with_normal_item::compute_normals_and_vertices(void)
{
    positions_points.clear();
    positions_lines.clear();
    color_points.clear();
    color_lines.clear();

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
    location[0] = glGetUniformLocation(rendering_program_lines, "mvp_matrix");
    location[1] = glGetUniformLocation(rendering_program_points, "mvp_matrix");
}
void Scene_points_with_normal_item::uniform_attrib(Viewer_interface* viewer, int mode) const
{
    GLfloat mvp_mat[16];

    //fills the MVP matrix.

    GLdouble d_mat[16];
    viewer->camera()->getModelViewProjectionMatrix(d_mat);
    //Convert the GLdoubles matrix in GLfloats
    for (int i=0; i<16; ++i){
        mvp_mat[i] = GLfloat(d_mat[i]);
    }

    if(mode ==0)
    {
        glUseProgram(rendering_program_lines);
        glUniformMatrix4fv(location[0], 1, GL_FALSE, mvp_mat);
    }
    else if(mode ==1)
    {
        glUseProgram(rendering_program_points);
        glUniformMatrix4fv(location[1], 1, GL_FALSE, mvp_mat);
    }
}


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
    emit itemChanged();
}

// Reset selection mark
void Scene_points_with_normal_item::resetSelection()
{
    // Un-select all points
    m_points->select(m_points->begin(), m_points->end(), false);
    emit itemChanged();
}
//Select duplicated points
void Scene_points_with_normal_item::selectDuplicates()
{
    std::set<Kernel::Point_3> unique_points;
    for (Point_set::Point_iterator ptit=m_points->begin(); ptit!=m_points->end();++ptit )
        if ( !unique_points.insert(*ptit).second )
            m_points->select(&(*ptit));
    emit itemChanged();
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

// Points OpenGL drawing in a display list
void Scene_points_with_normal_item::direct_draw() const
{
    Q_ASSERT(m_points != NULL);

    // Draw points
    m_points->gl_draw_vertices();
}

// Normals OpenGL drawing
void Scene_points_with_normal_item::draw_normals() const
{
    Q_ASSERT(m_points != NULL);

    // Draw normals
    Kernel::Sphere_3 region_of_interest = m_points->region_of_interest();
    float normal_length = (float)std::sqrt(region_of_interest.squared_radius() / 1000.);

    m_points->gl_draw_normals(normal_length);
}

void Scene_points_with_normal_item::draw_splats() const
{
    Q_ASSERT(m_points != NULL);

    // Draw splats
    bool points_have_normals = (m_points->begin() != m_points->end() &&
            m_points->begin()->normal() != CGAL::NULL_VECTOR);
    bool points_have_radii =   (m_points->begin() != m_points->end() &&
            m_points->begin()->radius() != 0);
    if(points_have_normals && points_have_radii)
    {
        m_points->gl_draw_splats();
    }
}

void Scene_points_with_normal_item::draw_edges(Viewer_interface* viewer) const
{
    glBindVertexArray(vao[0]);
    glUseProgram(rendering_program_lines);
    uniform_attrib(viewer,0);
    glDrawArrays(GL_LINES, 0, positions_lines.size()/3);
    glUseProgram(0);
    glBindVertexArray(0);
}
void Scene_points_with_normal_item::draw_points(Viewer_interface* viewer) const
{
    glBindVertexArray(vao[0]);
    glUseProgram(rendering_program_points);
    uniform_attrib(viewer,1);
    glDrawArrays(GL_POINTS, 0, positions_points.size()/3);
    glUseProgram(0);
    glBindVertexArray(0);

    GLfloat point_size;
    glGetFloatv(GL_POINT_SIZE, &point_size);
    glPointSize(4.f);
    glBindVertexArray(vao[1]);
    glUseProgram(rendering_program_points);
    uniform_attrib(viewer,1);
    glDrawArrays(GL_POINTS, 0, positions_selected_points.size()/3);
    glUseProgram(0);
    glBindVertexArray(0);
    glPointSize(point_size);
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
    Scene_item_with_display_list::setRenderingMode(m);
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
    initialize_buffers();
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
