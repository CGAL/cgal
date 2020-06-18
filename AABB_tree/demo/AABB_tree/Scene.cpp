#include "Scene.h"

#include <iostream>
#include <fstream>

#include <QtCore/qglobal.h>
#include <QString>
#include <QTextStream>
#include <QFileInfo>
#include <QInputDialog>

#include "Refiner.h"
//#include "render_edges.h"

#include <CGAL/Timer.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/subdivision_method_3.h>

#include <QOpenGLShader>
#include <QDebug>
#include "Viewer.h"


// constants
const int slow_distance_grid_size = 100;
const int fast_distance_grid_size = 20;
#define _SIGNED 0
#define _UNSIGNED 1

Scene::Scene()
    : m_frame (new ManipulatedFrame())
    , m_view_plane(false)
    , m_grid_size(slow_distance_grid_size)
    , m_cut_plane(NONE)
{
    m_pPolyhedron = NULL;

    // view options
    m_view_points = true;
    m_view_segments = true;
    m_view_polyhedron = true;

    // distance function
    m_red_ramp.build_red();
    m_blue_ramp.build_blue();
    m_max_distance_function = (FT)0.0;
    texture = new Texture(m_grid_size,m_grid_size);
    ready_to_cut = true;
    are_buffers_initialized = false;
    gl_init = false;

}

Scene::~Scene()
{
    delete m_pPolyhedron;
    delete m_frame;

    buffers[0].destroy();
    buffers[1].destroy();
    buffers[2].destroy();
    buffers[3].destroy();
    buffers[4].destroy();
    buffers[5].destroy();
    buffers[6].destroy();
    buffers[7].destroy();
    vao[0].destroy();
    vao[1].destroy();
    vao[2].destroy();
    vao[3].destroy();
    vao[4].destroy();
    vao[5].destroy();
    vao[6].destroy();


}

void Scene::compile_shaders()
{
    if(! buffers[0].create() || !buffers[1].create() || !buffers[2].create() || !buffers[3].create() || !buffers[4].create() || !buffers[5].create() || !buffers[6].create() || !buffers[7].create())
    {
        std::cerr<<"VBO Creation FAILED"<<std::endl;
    }

    if(!vao[0].create() || !vao[1].create() || !vao[2].create() || !vao[3].create() || !vao[4].create() || !vao[5].create() || !vao[6].create())
    {
        std::cerr<<"VAO Creation FAILED"<<std::endl;
    }


    //Vertex source code
    const char vertex_source[] =
    {
        "#version 120 \n"
        "attribute highp vec4 vertex;\n"
        "uniform highp mat4 mvp_matrix;\n"
        "uniform highp mat4 f_matrix;\n"
        "void main(void)\n"
        "{\n"
        "   gl_Position = mvp_matrix * f_matrix * vertex;\n"
        "}"
    };
    //Vertex source code
    const char fragment_source[] =
    {
        "#version 120 \n"
        "uniform highp vec4 color; \n"
        "void main(void) { \n"
        "gl_FragColor = color; \n"
        "} \n"
        "\n"
    };
    QOpenGLShader *vertex_shader = new QOpenGLShader(QOpenGLShader::Vertex);
    if(!vertex_shader->compileSourceCode(vertex_source))
    {
        std::cerr<<"Compiling vertex source FAILED"<<std::endl;
    }

    QOpenGLShader *fragment_shader= new QOpenGLShader(QOpenGLShader::Fragment);
    if(!fragment_shader->compileSourceCode(fragment_source))
    {
        std::cerr<<"Compiling fragmentsource FAILED"<<std::endl;
    }

    if(!rendering_program.addShader(vertex_shader))
    {
        std::cerr<<"adding vertex shader FAILED"<<std::endl;
    }
    if(!rendering_program.addShader(fragment_shader))
    {
        std::cerr<<"adding fragment shader FAILED"<<std::endl;
    }
    if(!rendering_program.link())
    {
        std::cerr<<"linking Program FAILED"<<std::endl;
    }
    rendering_program.bind();

    //Vertex source code
    const char tex_vertex_source[] =
    {
        "#version 120 \n"
        "attribute highp vec4 vertex;\n"
        "attribute highp vec2 tex_coord; \n"
        "uniform highp mat4 mvp_matrix;\n"
        "uniform highp mat4 f_matrix;\n"
        "varying highp vec2 texc;\n"
        "void main(void)\n"
        "{\n"
        "   gl_Position = mvp_matrix * f_matrix * vertex;\n"
        "    texc = tex_coord;\n"
        "}"
    };
    //Vertex source code
    const char tex_fragment_source[] =
    {
        "#version 120 \n"
        "uniform sampler2D texture;\n"
        "varying highp vec2 texc;\n"
        "void main(void) { \n"
        "gl_FragColor = texture2D(texture, texc.st);\n"
        "} \n"
        "\n"
    };
    QOpenGLShader *tex_vertex_shader = new QOpenGLShader(QOpenGLShader::Vertex);
    if(!tex_vertex_shader->compileSourceCode(tex_vertex_source))
    {
        std::cerr<<"Compiling vertex source FAILED"<<std::endl;
    }

    QOpenGLShader *tex_fragment_shader= new QOpenGLShader(QOpenGLShader::Fragment);
    if(!tex_fragment_shader->compileSourceCode(tex_fragment_source))
    {
        std::cerr<<"Compiling fragmentsource FAILED"<<std::endl;
    }

    if(!tex_rendering_program.addShader(tex_vertex_shader))
    {
        std::cerr<<"adding vertex shader FAILED"<<std::endl;
    }
    if(!tex_rendering_program.addShader(tex_fragment_shader))
    {
        std::cerr<<"adding fragment shader FAILED"<<std::endl;
    }
    if(!tex_rendering_program.link())
    {
        std::cerr<<"linking Program FAILED"<<std::endl;
    }
    tex_rendering_program.bind();

}

void Scene::initialize_buffers()
{
    //Points
    vao[0].bind();
    buffers[0].bind();
    buffers[0].allocate(pos_points.data(),
                        static_cast<int>(pos_points.size()*sizeof(float)));
    points_vertexLocation = rendering_program.attributeLocation("vertex");
    rendering_program.bind();
    rendering_program.enableAttributeArray(points_vertexLocation);
    rendering_program.setAttributeBuffer(points_vertexLocation,GL_FLOAT,0,3);
    buffers[0].release();
    rendering_program.release();
    vao[0].release();

    //Lines
    vao[1].bind();
    buffers[1].bind();
    buffers[1].allocate(pos_lines.data(),
                        static_cast<int>(pos_lines.size()*sizeof(float)));
    lines_vertexLocation = rendering_program.attributeLocation("vertex");
    rendering_program.bind();
    rendering_program.setAttributeBuffer(lines_vertexLocation,GL_FLOAT,0,3);
    buffers[1].release();
    rendering_program.enableAttributeArray(lines_vertexLocation);
    rendering_program.release();
    vao[1].release();

    //Polyhedron's edges
    vao[2].bind();
    buffers[2].bind();
    buffers[2].allocate(pos_poly.data(),
                        static_cast<int>(pos_poly.size()*sizeof(float)));
    poly_vertexLocation = rendering_program.attributeLocation("vertex");
    rendering_program.bind();
    rendering_program.setAttributeBuffer(poly_vertexLocation,GL_FLOAT,0,3);
    rendering_program.enableAttributeArray(poly_vertexLocation);
    buffers[2].release();
    rendering_program.release();
    vao[2].release();

    //cutting segments
    vao[3].bind();
    buffers[3].bind();
    buffers[3].allocate(pos_cut_segments.data(),
                        static_cast<int>(pos_cut_segments.size()*sizeof(float)));
    poly_vertexLocation = rendering_program.attributeLocation("vertex");
    rendering_program.bind();
    rendering_program.setAttributeBuffer(poly_vertexLocation,GL_FLOAT,0,3);
    rendering_program.enableAttributeArray(poly_vertexLocation);
    buffers[3].release();
    rendering_program.release();
    vao[3].release();

    //cutting plane
    vao[4].bind();
    buffers[4].bind();
    buffers[4].allocate(pos_plane.data(), static_cast<int>(pos_plane.size()*sizeof(float)));
    poly_vertexLocation = rendering_program.attributeLocation("vertex");
    rendering_program.bind();
    rendering_program.setAttributeBuffer(poly_vertexLocation,GL_FLOAT,0,3);
    rendering_program.enableAttributeArray(poly_vertexLocation);
    buffers[4].release();
    rendering_program.release();

    vao[4].release();

    //grid
    vao[5].bind();
    buffers[5].bind();
    buffers[5].allocate(pos_grid.data(), static_cast<int>(pos_grid.size()*sizeof(float)));
    poly_vertexLocation = rendering_program.attributeLocation("vertex");
    rendering_program.bind();
    rendering_program.setAttributeBuffer(poly_vertexLocation,GL_FLOAT,0,3);
    rendering_program.enableAttributeArray(poly_vertexLocation);
    buffers[5].release();
    rendering_program.release();
    vao[5].release();

    //cutting plane
    vao[6].bind();
    buffers[6].bind();
    buffers[6].allocate(pos_plane.data(), static_cast<int>(pos_plane.size()*sizeof(float)));
    poly_vertexLocation = tex_rendering_program.attributeLocation("vertex");
    tex_rendering_program.bind();
    tex_rendering_program.setAttributeBuffer(poly_vertexLocation,GL_FLOAT,0,3);
    tex_rendering_program.enableAttributeArray(poly_vertexLocation);
    buffers[6].release();
    tex_rendering_program.release();

    buffers[7].bind();
    buffers[7].allocate(tex_map.data(), static_cast<int>(tex_map.size()*sizeof(float)));
    tex_Location = tex_rendering_program.attributeLocation("tex_coord");
    tex_rendering_program.bind();
    tex_rendering_program.setAttributeBuffer(tex_Location,GL_FLOAT,0,2);
    tex_rendering_program.enableAttributeArray(tex_Location);
    buffers[7].release();
    tex_rendering_program.release();

    gl->glBindTexture(GL_TEXTURE_2D, textureId);
    gl->glTexImage2D(GL_TEXTURE_2D,
                 0,
                 GL_RGB,
                 texture->getWidth(),
                 texture->getHeight(),
                 0,
                 GL_RGB,
                 GL_UNSIGNED_BYTE,
                 texture->getData());
    gl->glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    gl->glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    gl->glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S,GL_CLAMP_TO_EDGE );
    gl->glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T,GL_CLAMP_TO_EDGE );
    vao[6].release();

    are_buffers_initialized = true;



}

void Scene::compute_elements(int mode)
{
    pos_points.resize(0);
    pos_lines.resize(0);
    pos_poly.resize(0);
    pos_cut_segments.resize(0);
    tex_map.resize(0);
    pos_grid.resize(66);
    pos_plane.resize(18);
    float diag = .6f * float(bbox_diag());
    //The Points
    {
        std::list<Point>::iterator pit;
        for(pit = m_points.begin(); pit != m_points.end(); pit++)
        {
            const Point& p = *pit;
            pos_points.push_back(p.x());
            pos_points.push_back(p.y());
            pos_points.push_back(p.z());
        }
    }
    //The Segements
    {
        std::list<Segment>::iterator sit;
        for(sit = m_segments.begin(); sit != m_segments.end(); sit++)
        {
            const Segment& s = *sit;
            const Point& p = s.source();
            const Point& q = s.target();

            pos_lines.push_back(p.x());
            pos_lines.push_back(p.y());
            pos_lines.push_back(p.z());

            pos_lines.push_back(q.x());
            pos_lines.push_back(q.y());
            pos_lines.push_back(q.z());
        }
    }
    //The Polygon's edges
    {
        Polyhedron::Edge_iterator he;
        for(he = m_pPolyhedron->edges_begin();
            he != m_pPolyhedron->edges_end();
            he++)
        {
            const Point& a = he->vertex()->point();
            const Point& b = he->opposite()->vertex()->point();
            pos_poly.push_back(a.x());
            pos_poly.push_back(a.y());
            pos_poly.push_back(a.z());

            pos_poly.push_back(b.x());
            pos_poly.push_back(b.y());
            pos_poly.push_back(b.z());
        }
    }
    //The cutting segments
    {
        for ( std::vector<Segment>::const_iterator csit = m_cut_segments.begin(),
              end = m_cut_segments.end() ; csit != end ; ++csit )
        {
            const Point& a = csit->source();
            const Point& b = csit->target();

            pos_cut_segments.push_back(a.x());
            pos_cut_segments.push_back(a.y());
            pos_cut_segments.push_back(a.z());

            pos_cut_segments.push_back(b.x());
            pos_cut_segments.push_back(b.y());
            pos_cut_segments.push_back(b.z());
        }
    }
    //The cutting plane
    {

        pos_plane[0]= -diag; pos_plane[1]=-diag; pos_plane[2]=0.0;
        pos_plane[3]= -diag; pos_plane[4]= diag; pos_plane[5]=0.;
        pos_plane[6]=  diag; pos_plane[7]= diag; pos_plane[8]=0.;
        pos_plane[9]= -diag; pos_plane[10]= -diag; pos_plane[11]=0.;
        pos_plane[12]= diag;    pos_plane[13]= diag; pos_plane[14]= 0.;
        pos_plane[15]= diag;    pos_plane[16]= -diag; pos_plane[17]= 0.;

        //UV Mapping
        tex_map.push_back(-0.11f);
        tex_map.push_back(-0.11f);

        tex_map.push_back(-0.11f);
        tex_map.push_back(1.11f);

        tex_map.push_back(1.11f);
        tex_map.push_back(1.11f);

        tex_map.push_back(-0.11f);
        tex_map.push_back(-0.11f);

        tex_map.push_back(1.11f);
        tex_map.push_back(1.11f);

        tex_map.push_back(1.11f);
        tex_map.push_back(-0.11f);





    }
    //The grid
    {
        float z = 0;
        float x = (2 * diag)/10.0;
        float y = (2 * diag)/10.0;
        for(int u = 0; u < 11; u++)
        {

            pos_grid.push_back(-diag + x* u);
            pos_grid.push_back(-diag);
            pos_grid.push_back(z);

            pos_grid.push_back(-diag + x* u);
            pos_grid.push_back(diag);
            pos_grid.push_back(z);
        }
        for(int v=0; v<11; v++)
        {

            pos_grid.push_back(-diag);
            pos_grid.push_back(-diag + v * y);
            pos_grid.push_back(z);

            pos_grid.push_back(diag);
            pos_grid.push_back(-diag + v * y);
            pos_grid.push_back(z);

        }

    }
    //The texture
    switch(mode)
    {
    case _SIGNED:
        for( int i=0 ; i < texture->getWidth(); i++ )
        {
            for( int j=0 ; j < texture->getHeight() ; j++)
            {
                compute_texture(i,j,m_red_ramp,m_blue_ramp);
            }
        }
        break;
    case _UNSIGNED:
        for( int i=0 ; i < texture->getWidth(); i++ )
        {
            for( int j=0 ; j < texture->getHeight() ; j++)
            {
                compute_texture(i,j,m_thermal_ramp,m_thermal_ramp);
            }
        }
        break;}
    sampler_location = tex_rendering_program.attributeLocation("texture");
}

void Scene::compute_texture(int i, int j,Color_ramp pos_ramp ,Color_ramp neg_ramp)
{


    const FT& d00 = m_distance_function[i][j].second;
    // determines grey level
    unsigned int i00 = 255-(unsigned)(255.0 * (double)std::fabs(d00) / m_max_distance_function);

    if(d00 > 0.0)
        texture->setData(i,j,pos_ramp.r(i00),pos_ramp.g(i00),pos_ramp.b(i00));
    else
        texture->setData(i,j,neg_ramp.r(i00),neg_ramp.g(i00),neg_ramp.b(i00));


}

void Scene::attrib_buffers(CGAL::QGLViewer* viewer)
{
    QMatrix4x4 mvpMatrix;
    double mat[16];
    viewer->camera()->getModelViewProjectionMatrix(mat);
    for(int i=0; i < 16; i++)
    {
        mvpMatrix.data()[i] = (float)mat[i];
    }
    rendering_program.bind();
    mvpLocation = rendering_program.uniformLocation("mvp_matrix");
    fLocation = rendering_program.uniformLocation("f_matrix");
    colorLocation = rendering_program.uniformLocation("color");
    rendering_program.setUniformValue(mvpLocation, mvpMatrix);
    rendering_program.release();

    tex_rendering_program.bind();
    tex_mvpLocation = tex_rendering_program.uniformLocation("mvp_matrix");
    tex_fLocation = tex_rendering_program.uniformLocation("f_matrix");
    tex_rendering_program.setUniformValue(tex_mvpLocation, mvpMatrix);
    tex_rendering_program.release();
}

void Scene::changed()
{
    if(m_cut_plane == UNSIGNED_FACETS || m_cut_plane == UNSIGNED_EDGES)
        compute_elements(_UNSIGNED);
    else
        compute_elements(_SIGNED);
    ready_to_cut=false;
    are_buffers_initialized = false;

}

int Scene::open(QString filename)
{
    QTextStream cerr(stderr);
    cerr << QString("Opening file \"%1\"\n").arg(filename);
    QApplication::setOverrideCursor(QCursor(::Qt::WaitCursor));

    QFileInfo fileinfo(filename);
    std::ifstream in(filename.toUtf8());

    if(!in || !fileinfo.isFile() || ! fileinfo.isReadable())
    {
        std::cerr << "unable to open file" << std::endl;
        QApplication::restoreOverrideCursor();
        return -1;
    }

    if(m_pPolyhedron != NULL)
        delete m_pPolyhedron;

    // allocate new polyhedron
    m_pPolyhedron = new Polyhedron;
    in >> *m_pPolyhedron;
    if(!in)
    {
        std::cerr << "invalid OFF file" << std::endl;
        QApplication::restoreOverrideCursor();

        delete m_pPolyhedron;
        m_pPolyhedron = NULL;

        return -1;
    }

    // clear tree
    clear_internal_data();

    QApplication::restoreOverrideCursor();
    changed();
    return 0;
}

void Scene::update_bbox()
{
    std::cout << "Compute bbox...";
    m_bbox = Bbox();

    if(m_pPolyhedron == NULL)
    {
        std::cout << "failed (no polyhedron)." << std::endl;
        return;
    }

    if(m_pPolyhedron->empty())
    {
        std::cout << "failed (empty polyhedron)." << std::endl;
        return;
    }

    Polyhedron::Point_iterator it = m_pPolyhedron->points_begin();
    m_bbox = (*it).bbox();
    for(; it != m_pPolyhedron->points_end();it++)
        m_bbox = m_bbox + (*it).bbox();
    std::cout << "done (" << m_pPolyhedron->size_of_facets()
              << " facets)" << std::endl;
}

void Scene::draw(CGAL::QGLViewer* viewer)
{
    if(!gl_init)
        initGL();
    if(!are_buffers_initialized)
        initialize_buffers();
    gl->glEnable(GL_DEPTH_TEST);
    QColor color;
    QMatrix4x4 fMatrix;
    fMatrix.setToIdentity();
    if(m_view_polyhedron && pos_poly.size()>0)
    {
        vao[2].bind();
        attrib_buffers(viewer);
        rendering_program.bind();
        color.setRgbF(0.0,0.0,0.0);
        rendering_program.setUniformValue(colorLocation, color);
        rendering_program.setUniformValue(fLocation, fMatrix);
        gl->glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(pos_poly.size()/3));
        rendering_program.release();
        vao[2].release();
    }
    if(m_view_points && pos_points.size()>0)
    {
        vao[0].bind();
        attrib_buffers(viewer);
        rendering_program.bind();
        color.setRgbF(0.7,0.0,0.0);
        rendering_program.setUniformValue(colorLocation, color);
        rendering_program.setUniformValue(fLocation, fMatrix);
        gl->glDrawArrays(GL_POINTS, 0, static_cast<GLsizei>(pos_points.size()/3));
        rendering_program.release();
        vao[0].release();
    }

    if(m_view_segments && pos_lines.size()>0)
    {
        vao[1].bind();
        attrib_buffers(viewer);
        rendering_program.bind();
        color.setRgbF(0.0,0.7,0.0);
        rendering_program.setUniformValue(colorLocation, color);
        rendering_program.setUniformValue(fLocation, fMatrix);
        gl->glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(pos_lines.size()/3));
        rendering_program.release();
        vao[1].release();
    }
    if (m_view_plane && pos_plane.size()>0)
    {
        switch( m_cut_plane )
        {
        case UNSIGNED_EDGES:
        case UNSIGNED_FACETS:
        case SIGNED_FACETS:

            gl->glActiveTexture(GL_TEXTURE0);
            gl->glBindTexture(GL_TEXTURE_2D, textureId);

            for(int i=0; i< 16 ; i++)
                fMatrix.data()[i] =  m_frame->matrix()[i];
            vao[6].bind();
            attrib_buffers(viewer);
            tex_rendering_program.bind();
            tex_rendering_program.setUniformValue(tex_fLocation, fMatrix);

            gl->glDrawArrays(GL_TRIANGLES, 0,static_cast<GLsizei>(pos_plane.size()/3));
            tex_rendering_program.release();
            vao[6].release();
            break;

        case CUT_SEGMENTS:

            //cutting_segments
            fMatrix.setToIdentity();
            gl->glLineWidth(2.0f);
            vao[3].bind();
            attrib_buffers(viewer);
            rendering_program.bind();
            color.setRgbF(1.0,0.0,0.0);
            rendering_program.setUniformValue(colorLocation, color);
            rendering_program.setUniformValue(fLocation, fMatrix);
            gl->glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(pos_cut_segments.size()/3));
            gl->glLineWidth(1.0f);
            rendering_program.release();
            vao[3].release();
            //grid
            for(int i=0; i< 16 ; i++)
                fMatrix.data()[i] =  m_frame->matrix()[i];
            vao[5].bind();
            attrib_buffers(viewer);
            rendering_program.bind();
            color.setRgbF(.6f, .6f, .6f);
            rendering_program.setUniformValue(colorLocation, color);
            rendering_program.setUniformValue(fLocation, fMatrix);
            gl->glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(pos_grid.size()/3));
            rendering_program.release();
            vao[5].release();

            //cutting_plane
            // for(int i=0; i< 16 ; i++)
            //     fMatrix.data()[i] =  m_frame->matrix()[i];
            gl->glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
            gl->glEnable(GL_BLEND);
            vao[4].bind();
            attrib_buffers(viewer);
            rendering_program.bind();
            color.setRgbF(.6f, .85f, 1.f, .65f);
            rendering_program.setUniformValue(colorLocation, color);
            rendering_program.setUniformValue(fLocation, fMatrix);
            gl->glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(pos_plane.size()/3));
            gl->glDisable(GL_BLEND);
            rendering_program.release();
            vao[4].release();

            break;
        case NONE: // do nothing
            break;
        }
    }


}

FT Scene::random_in(const double a,
                    const double b)
{
    double r = rand() / (double)RAND_MAX;
    return (FT)(a + (b - a) * r);
}

Point Scene::random_point(const CGAL::Bbox_3& bbox)
{
    FT x = random_in(bbox.xmin(),bbox.xmax());
    FT y = random_in(bbox.ymin(),bbox.ymax());
    FT z = random_in(bbox.zmin(),bbox.zmax());
    return Point(x,y,z);
}

Vector Scene::random_vector()
{
    FT x = random_in(0.0,1.0);
    FT y = random_in(0.0,1.0);
    FT z = random_in(0.0,1.0);
    return Vector(x,y,z);
}

Ray Scene::random_ray(const CGAL::Bbox_3& bbox)
{
    Point p = random_point(bbox);
    Point q = random_point(bbox);
    return Ray(p,q);
}

Segment Scene::random_segment(const CGAL::Bbox_3& bbox)
{
    Point p = random_point(bbox);
    Point q = random_point(bbox);
    return Segment(p,q);
}

Line Scene::random_line(const CGAL::Bbox_3& bbox)
{
    Point p = random_point(bbox);
    Point q = random_point(bbox);
    return Line(p,q);
}

Plane Scene::random_plane(const CGAL::Bbox_3& bbox)
{
    Point p = random_point(bbox);
    Vector vec = random_vector();
    return Plane(p,vec);
}

Plane Scene::frame_plane() const
{
    const CGAL::qglviewer::Vec& pos = m_frame->position();
    const CGAL::qglviewer::Vec& n = m_frame->inverseTransformOf(CGAL::qglviewer::Vec(0.f, 0.f, 1.f));

    return Plane(n[0], n[1],  n[2], - n * pos);
}

Aff_transformation Scene::frame_transformation() const
{
    const ::GLdouble* m = m_frame->matrix();

    // OpenGL matrices are row-major matrices
    return Aff_transformation (m[0], m[4], m[8], m[12],
            m[1], m[5], m[9], m[13],
            m[2], m[6], m[10], m[14]);
}

FT Scene::bbox_diag() const
{
    double dx = m_bbox.xmax()-m_bbox.xmin();
    double dy = m_bbox.ymax()-m_bbox.ymin();
    double dz = m_bbox.zmax()-m_bbox.zmin();

    return FT(std::sqrt(dx*dx + dy*dy + dz*dz));
}

void Scene::build_facet_tree()
{
    if ( NULL == m_pPolyhedron )
    {
        std::cerr << "Build facet tree failed: load polyhedron first." << std::endl;
        return;
    }

    // Don't rebuild tree if it is already built
    if ( !m_facet_tree.empty() ) { return; }

    // build tree
    CGAL::Timer timer;
    timer.start();
    std::cout << "Construct Facet AABB tree...";
    m_facet_tree.rebuild(faces(*m_pPolyhedron).first, faces(*m_pPolyhedron).second,*m_pPolyhedron);
    std::cout << "done (" << timer.time() << " s)" << std::endl;
}

void Scene::build_edge_tree()
{
    if ( NULL == m_pPolyhedron )
    {
        std::cerr << "Build edge tree failed: load polyhedron first." << std::endl;
        return;
    }

    // Don't rebuild tree if it is already built
    if ( !m_edge_tree.empty() ) { return; }

    // build tree
    CGAL::Timer timer;
    timer.start();
    std::cout << "Construct Edge AABB tree...";
    m_edge_tree.rebuild(edges(*m_pPolyhedron).first,edges(*m_pPolyhedron).second,*m_pPolyhedron);
    std::cout << "done (" << timer.time() << " s)" << std::endl;
}

void Scene::clear_internal_data()
{
    m_facet_tree.clear();
    m_edge_tree.clear();

    clear_points();
    clear_segments();
    clear_cutting_plane();
}

void Scene::clear_cutting_plane()
{
    m_cut_segments.clear();
    m_cut_plane = NONE;

    deactivate_cutting_plane();
    changed();
}

void Scene::update_grid_size()
{
    m_grid_size = m_fast_distance ? fast_distance_grid_size
                                  : slow_distance_grid_size;
    texture = new Texture(m_grid_size,m_grid_size);
}

void Scene::generate_points_in(const unsigned int nb_points,
                               const double vmin,
                               const double vmax)
{
    if(m_pPolyhedron == NULL)
    {
        std::cout << "Load polyhedron first." << std::endl;
        return;
    }

    typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Primitive;
    typedef CGAL::AABB_traits<Kernel, Primitive> Traits;
    typedef CGAL::AABB_tree<Traits> Tree;

    std::cout << "Construct AABB tree...";
    Tree tree(faces(*m_pPolyhedron).first, faces(*m_pPolyhedron).second, *m_pPolyhedron);
    std::cout << "done." << std::endl;

    CGAL::Timer timer;
    timer.start();
    std::cout << "Generate " << nb_points << " points in interval ["
              << vmin << ";" << vmax << "]";

    unsigned int nb_trials = 0;
    Vector vec = random_vector();
    while(m_points.size() < nb_points)
    {
        Point p = random_point(tree.bbox());

        // measure distance
        FT signed_distance = std::sqrt(tree.squared_distance(p));

        // measure sign
        Ray ray(p,vec);
        int nb_intersections = (int)tree.number_of_intersected_primitives(ray);
        if(nb_intersections % 2 != 0)
            signed_distance *= -1.0;

        if(signed_distance >= vmin &&
                signed_distance <= vmax)
        {
            m_points.push_back(p);
            if(m_points.size()%(nb_points/10) == 0)
                std::cout << "."; // ASCII progress bar
        }
        nb_trials++;
    }
    double speed = (double)nb_trials / timer.time();
    std::cout << "done (" << nb_trials << " trials, "
              << timer.time() << " s, "
              << speed << " queries/s)" << std::endl;
    changed();
}


void Scene::generate_inside_points(const unsigned int nb_points)
{
    if(m_pPolyhedron == NULL)
    {
        std::cout << "Load polyhedron first." << std::endl;
        return;
    }

    typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Primitive;
    typedef CGAL::AABB_traits<Kernel, Primitive> Traits;
    typedef CGAL::AABB_tree<Traits> Tree;

    std::cout << "Construct AABB tree...";
    Tree tree(faces(*m_pPolyhedron).first, faces(*m_pPolyhedron).second,*m_pPolyhedron);
    std::cout << "done." << std::endl;

    CGAL::Timer timer;
    timer.start();
    std::cout << "Generate " << nb_points << " inside points";

    unsigned int nb_trials = 0;
    Vector vec = random_vector();
    while(m_points.size() < nb_points)
    {
        Point p = random_point(tree.bbox());
        Ray ray(p,vec);
        int nb_intersections = (int)tree.number_of_intersected_primitives(ray);
        if(nb_intersections % 2 != 0)
        {
            m_points.push_back(p);
            if(m_points.size()%(nb_points/10) == 0)
                std::cout << "."; // ASCII progress bar
        }
        nb_trials++;
    }
    double speed = (double)nb_trials / timer.time();
    std::cout << "done (" << nb_trials << " trials, "
              << timer.time() << " s, "
              << speed << " queries/s)" << std::endl;
    changed();
}

void Scene::generate_boundary_segments(const unsigned int nb_slices)
{
    if(m_pPolyhedron == NULL)
    {
        std::cout << "Load polyhedron first." << std::endl;
        return;
    }

    typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Primitive;
    typedef CGAL::AABB_traits<Kernel, Primitive> Traits;
    typedef CGAL::AABB_tree<Traits> Tree;
    typedef Tree::Object_and_primitive_id Object_and_primitive_id;

    std::cout << "Construct AABB tree...";
    Tree tree(faces(*m_pPolyhedron).first,faces(*m_pPolyhedron).second,*m_pPolyhedron);
    std::cout << "done." << std::endl;

    CGAL::Timer timer;
    timer.start();
    std::cout << "Generate boundary segments from " << nb_slices << " slices: ";

    Vector normal((FT)0.0,(FT)0.0,(FT)1.0);
    unsigned int i;

    const double dz = m_bbox.zmax() - m_bbox.zmin();
    for(i=0;i<nb_slices;i++)
    {
        FT z = m_bbox.zmin() + (FT)i / (FT)nb_slices * dz;
        Point p((FT)0.0, (FT)0.0, z);
        Plane plane(p,normal);

        std::list<Object_and_primitive_id> intersections;
        tree.all_intersections(plane,std::back_inserter(intersections));

        std::list<Object_and_primitive_id>::iterator it;
        for(it = intersections.begin();
            it != intersections.end();
            it++)
        {
            Object_and_primitive_id op = *it;
            CGAL::Object object = op.first;
            Segment segment;
            if(CGAL::assign(segment,object))
                m_segments.push_back(segment);
        }
    }
    std::cout << m_segments.size() << " segments, " << timer.time() << " s." << std::endl;
    changed();
}

void Scene::generate_boundary_points(const unsigned int nb_points)
{
    if(m_pPolyhedron == NULL)
    {
        std::cout << "Load polyhedron first." << std::endl;
        return;
    }

    typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Primitive;
    typedef CGAL::AABB_traits<Kernel, Primitive> Traits;
    typedef CGAL::AABB_tree<Traits> Tree;
    typedef Tree::Object_and_primitive_id Object_and_primitive_id;

    std::cout << "Construct AABB tree...";
    Tree tree(faces(*m_pPolyhedron).first, faces(*m_pPolyhedron).second,*m_pPolyhedron);
    std::cout << "done." << std::endl;

    CGAL::Timer timer;
    timer.start();
    std::cout << "Generate boundary points: ";

    unsigned int nb = 0;
    unsigned int nb_lines = 0;
    while(nb < nb_points)
    {
        Line line = random_line(tree.bbox());

        std::list<Object_and_primitive_id> intersections;
        tree.all_intersections(line,std::back_inserter(intersections));
        nb_lines++;

        std::list<Object_and_primitive_id>::iterator it;
        for(it = intersections.begin();
            it != intersections.end();
            it++)
        {
            Object_and_primitive_id op = *it;
            CGAL::Object object = op.first;
            Point point;
            if(CGAL::assign(point,object))
            {
                m_points.push_back(point);
                nb++;
            }
        }
    }
    std::cout << nb_lines << " line queries, " << timer.time() << " s." << std::endl;
    changed();
}

void Scene::generate_edge_points(const unsigned int nb_points)
{
    if(m_pPolyhedron == NULL)
    {
        std::cout << "Load polyhedron first." << std::endl;
        return;
    }

    typedef CGAL::AABB_halfedge_graph_segment_primitive<Polyhedron> Primitive;
    typedef CGAL::AABB_traits<Kernel, Primitive> Traits;
    typedef CGAL::AABB_tree<Traits> Tree;
    typedef Tree::Object_and_primitive_id Object_and_primitive_id;

    std::cout << "Construct AABB tree...";
    Tree tree( CGAL::edges(*m_pPolyhedron).first,
               CGAL::edges(*m_pPolyhedron).second,
               *m_pPolyhedron);
    std::cout << "done." << std::endl;

    CGAL::Timer timer;
    timer.start();
    std::cout << "Generate edge points: ";

    unsigned int nb = 0;
    unsigned int nb_planes = 0;
    while(nb < nb_points)
    {
        Plane plane = random_plane(tree.bbox());

        std::list<Object_and_primitive_id> intersections;
        tree.all_intersections(plane,std::back_inserter(intersections));
        nb_planes++;

        std::list<Object_and_primitive_id>::iterator it;
        for(it = intersections.begin();
            it != intersections.end();
            it++)
        {
            Object_and_primitive_id op = *it;
            CGAL::Object object = op.first;
            Point point;
            if(CGAL::assign(point,object))
            {
                m_points.push_back(point);
                nb++;
            }
        }
    }
    std::cout << nb_planes << " plane queries, " << timer.time() << " s." << std::endl;
    changed();
}


template <typename Tree>
void Scene::compute_distance_function(const Tree& tree)
{
    // Get transformation
    Aff_transformation t = frame_transformation();
    m_max_distance_function = FT(0);
    FT diag = bbox_diag();

    const FT dx = diag;
    const FT dy = diag;
    const FT z (0);
    const FT fd =  FT(2);
    for(int i=0 ; i<m_grid_size ; ++i)
    {
        FT x = -diag/fd + FT(i)/FT(m_grid_size) * dx;

        for(int j=0 ; j<m_grid_size ; ++j)
        {
            FT y = -diag/fd + FT(j)/FT(m_grid_size) * dy;
            Point query = t( Point(x,y,z) );
            FT dist = CGAL::sqrt( tree.squared_distance(query) );

            m_distance_function[i][j] = Point_distance(query,dist);
            m_max_distance_function = (std::max)(dist, m_max_distance_function);
        }
    }
}

template <typename Tree>
void Scene::sign_distance_function(const Tree& tree)
{
    typedef typename Tree::size_type size_type;
    Vector random_vec = random_vector();

    for(int i=0 ; i<m_grid_size ; ++i)
    {
        for(int j=0 ; j<m_grid_size ; ++j)
        {
            const Point& p = m_distance_function[i][j].first;
            const FT unsigned_distance = m_distance_function[i][j].second;

            // get sign through ray casting (random vector)
            Ray ray(p, random_vec);
            size_type nbi = tree.number_of_intersected_primitives(ray);

            FT sign ( (nbi&1) == 0 ? 1 : -1);
            m_distance_function[i][j].second = sign * unsigned_distance;
        }
    }
    changed();
}


void Scene::unsigned_distance_function()
{
    // Build tree (if build fail, exit)
    build_facet_tree();
    if ( m_facet_tree.empty() ) { return; }

    compute_distance_function(m_facet_tree);

    m_cut_plane = UNSIGNED_FACETS;
    changed();
}


void Scene::unsigned_distance_function_to_edges()
{
    // Build tree (if build fail, exit)
    build_edge_tree();
    if ( m_edge_tree.empty() ) { return; }

    compute_distance_function(m_edge_tree);

    m_cut_plane = UNSIGNED_EDGES;
    changed();
}


void Scene::signed_distance_function()
{
    // Build tree (if build fail, exit)
    build_facet_tree();
    if ( m_facet_tree.empty() ) { return; }

    compute_distance_function(m_facet_tree);
    sign_distance_function(m_facet_tree);

    m_cut_plane = SIGNED_FACETS;
    changed();
}


void Scene::cut_segment_plane()
{
    // Build tree (if build fail, exit)
    build_facet_tree();
    if ( m_facet_tree.empty() ) { return; }

    Plane plane = frame_plane();

    // Compute intersections
    typedef std::vector<Facet_tree::Object_and_primitive_id> Intersections;
    Intersections intersections;
    m_facet_tree.all_intersections(plane, std::back_inserter(intersections));

    // Fill data structure
    m_cut_segments.clear();
    for ( Intersections::iterator it = intersections.begin(),
          end = intersections.end() ; it != end ; ++it )
    {
        const Segment* inter_seg = CGAL::object_cast<Segment>(&(it->first));

        if ( NULL != inter_seg )
        {
            m_cut_segments.push_back(*inter_seg);
        }
    }

    m_cut_plane = CUT_SEGMENTS;
    changed();
}
void Scene::updateCutPlane()
{
  ready_to_cut = true;
       QTimer::singleShot(0,this,SLOT(cutting_plane()));
}

void Scene::cutting_plane(bool override)
{
    if(ready_to_cut || override)
    {
        switch( m_cut_plane )
        {
        case UNSIGNED_FACETS:
            return unsigned_distance_function();
        case SIGNED_FACETS:
            return signed_distance_function();
        case UNSIGNED_EDGES:
            return unsigned_distance_function_to_edges();
        case CUT_SEGMENTS:
            return cut_segment_plane();
        case NONE: // do nothing
            return;
        }

        // Should not be here
        std::cerr << "Unknown cut_plane type" << std::endl;
        CGAL_assertion(false);
    }
}

void Scene::toggle_view_poyhedron()
{
    m_view_polyhedron = !m_view_polyhedron;
}

void Scene::toggle_view_segments()
{
    m_view_segments = !m_view_segments;
}

void Scene::toggle_view_points()
{
    m_view_points = !m_view_points;
}

void Scene::toggle_view_plane()
{
    m_view_plane = !m_view_plane;
}

void Scene::refine_bisection(const FT max_sqlen)
{
    if(m_pPolyhedron == NULL)
    {
        std::cout << "Load polyhedron first." << std::endl;
        return;
    }
    std::cout << "Refine through recursive longest edge bisection...";
    Refiner<Kernel,Polyhedron> refiner(m_pPolyhedron);
    refiner(max_sqlen);
    std::cout << "done (" << m_pPolyhedron->size_of_facets() << " facets)" << std::endl;

    clear_internal_data();
}

void Scene::refine_loop()
{
    if(m_pPolyhedron == NULL)
    {
        std::cout << "Load polyhedron first." << std::endl;
        return;
    }
    std::cout << "Loop subdivision...";
    CGAL::Subdivision_method_3::Loop_subdivision(*m_pPolyhedron);
    std::cout << "done (" << m_pPolyhedron->size_of_facets() << " facets)" << std::endl;

    clear_internal_data();
}


void Scene::activate_cutting_plane()
{
    connect(m_frame, SIGNAL(modified()), this, SLOT(updateCutPlane()));
    m_view_plane = true;
}

void Scene::deactivate_cutting_plane()
{
    disconnect(m_frame, SIGNAL(modified()), this, SLOT(updateCutPlane()));
    m_view_plane = false;
}
void Scene::initGL()
{
    gl = new QOpenGLFunctions_2_1();
   if(!gl->initializeOpenGLFunctions())
    {
        qFatal("ERROR : OpenGL Functions not initialized. Check your OpenGL Verison (should be >=3.3)");
        exit(1);
    }

    gl->glGenTextures(1, &textureId);
    compile_shaders();
    gl_init = true;
}
