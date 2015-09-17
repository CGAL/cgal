#include "Scene_implicit_function_item.h"
#include <QColor>
#include <map>
#include <CGAL/gl.h>
#include <CGAL/Simple_cartesian.h>
#include <QGLViewer/manipulatedFrame.h>

#include "Color_ramp.h"
#include <Viewer_interface.h>

#include <CGAL/double.h>

inline
bool is_nan(double d)
{
    return !CGAL::Is_valid<double>()( d );
}

void Scene_implicit_function_item::initialize_buffers(Viewer_interface *viewer = 0) const
{
    if(GLuint(-1) == textureId) {
        viewer->glGenTextures(1, &textureId);
    }

    //vao fot the cutting plane
    {
        program = getShaderProgram(PROGRAM_WITH_TEXTURE, viewer);
        program->bind();
        vaos[0]->bind();


        buffers[0].bind();
        buffers[0].allocate(positions_tex_quad.data(),
                            static_cast<int>(positions_tex_quad.size()*sizeof(float)));
        program->enableAttributeArray("vertex");
        program->setAttributeBuffer("vertex",GL_FLOAT,0,3);
        buffers[0].release();

        buffers[1].bind();
        buffers[1].allocate(texture_map.data(),
                            static_cast<int>(texture_map.size()*sizeof(float)));
        program->enableAttributeArray("v_texCoord");
        program->setAttributeBuffer("v_texCoord",GL_FLOAT,0,2);
        buffers[1].release();
        program->setAttributeValue("normal", QVector3D(0,0,0));

        program->release();
        vaos[0]->release();
    }
    //vao fot the bbox
    {
        program = getShaderProgram(PROGRAM_WITHOUT_LIGHT, viewer);
        program->bind();
        vaos[1]->bind();


        buffers[2].bind();
        buffers[2].allocate(positions_cube.data(),
                            static_cast<int>(positions_cube.size()*sizeof(float)));
        program->enableAttributeArray("vertex");
        program->setAttributeBuffer("vertex",GL_FLOAT,0,3);
        buffers[2].release();

        program->setAttributeValue("colors", QVector3D(0,0,0));
        program->release();
        vaos[1]->release();
    }
    //vao fot the grid
    {
        program = getShaderProgram(PROGRAM_WITHOUT_LIGHT, viewer);
        program->bind();
        vaos[2]->bind();


        buffers[3].bind();
        buffers[3].allocate(positions_grid.data(),
                            static_cast<int>(positions_grid.size()*sizeof(float)));
        program->enableAttributeArray("vertex");
        program->setAttributeBuffer("vertex",GL_FLOAT,0,3);
        buffers[3].release();
        program->setAttributeValue("colors", QVector3D(0.6f, 0.6f, 0.6f));
        program->release();
        vaos[2]->release();
    }



    viewer->glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    viewer->glBindTexture(GL_TEXTURE_2D, textureId);
    viewer->glTexImage2D(GL_TEXTURE_2D,
                 0,
                 GL_RGB,
                 texture->getWidth(),
                 texture->getHeight(),
                 0,
                 GL_RGB,
                 GL_UNSIGNED_BYTE,
                 texture->getData());
       viewer->glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
       viewer->glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
       viewer->glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S,GL_CLAMP_TO_EDGE );
       viewer->glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T,GL_CLAMP_TO_EDGE );

       are_buffers_filled = true;
}

void Scene_implicit_function_item::compute_vertices_and_texmap(void)
{
    positions_tex_quad.resize(0);
    positions_cube.resize(0);
    positions_grid.resize(0);
    texture_map.resize(0);

    const Bbox& b = bbox();
    float x,y,z;
    z = 0;
    x = (b.xmax-b.xmin)/10.0;
    y = (b.ymax-b.ymin)/10.0;
    // The Quad
    {


        //A
        positions_tex_quad.push_back(b.xmin);
        positions_tex_quad.push_back(b.ymin);
        positions_tex_quad.push_back(z);


        //B
        positions_tex_quad.push_back(b.xmin);
        positions_tex_quad.push_back(b.ymax);
        positions_tex_quad.push_back(z);


        //C
        positions_tex_quad.push_back(b.xmax);
        positions_tex_quad.push_back(b.ymax);
        positions_tex_quad.push_back(z);



        //A
        positions_tex_quad.push_back(b.xmin);
        positions_tex_quad.push_back(b.ymin);
        positions_tex_quad.push_back(z);


        //C
        positions_tex_quad.push_back(b.xmax);
        positions_tex_quad.push_back(b.ymax);
        positions_tex_quad.push_back(z);


        //D
        positions_tex_quad.push_back(b.xmax);
        positions_tex_quad.push_back(b.ymin);
        positions_tex_quad.push_back(z);


        //UV Mapping x2 but I don't know why.
        texture_map.push_back(0.0);
        texture_map.push_back(0.0);

        texture_map.push_back(0.0);
        texture_map.push_back(1.0);

        texture_map.push_back(1.0);
        texture_map.push_back(1.0);

        texture_map.push_back(0.0);
        texture_map.push_back(0.0);

        texture_map.push_back(1.0);
        texture_map.push_back(1.0);

        texture_map.push_back(1.0);
        texture_map.push_back(0.0);



    }
    //The grid
    {

        for(int u = 0; u < 11; u++)
        {

            positions_grid.push_back(b.xmin + x* u);
            positions_grid.push_back(b.ymin);
            positions_grid.push_back(z);

            positions_grid.push_back(b.xmin + x* u);
            positions_grid.push_back(b.ymax);
            positions_grid.push_back(z);
        }
        for(int v=0; v<11; v++)
        {

            positions_grid.push_back(b.xmin);
            positions_grid.push_back(b.ymin + v * y);
            positions_grid.push_back(z);

            positions_grid.push_back(b.xmax);
            positions_grid.push_back(b.ymin + v * y);
            positions_grid.push_back(z);
        }

    }
    //the Box
    {

        positions_cube.push_back(b.xmin);
        positions_cube.push_back(b.ymin);
        positions_cube.push_back(b.zmin);

        positions_cube.push_back(b.xmin);
        positions_cube.push_back(b.ymin);
        positions_cube.push_back(b.zmax);


        positions_cube.push_back(b.xmin);
        positions_cube.push_back(b.ymin);
        positions_cube.push_back(b.zmin);


        positions_cube.push_back(b.xmin);
        positions_cube.push_back(b.ymax);
        positions_cube.push_back(b.zmin);


        positions_cube.push_back(b.xmin);
        positions_cube.push_back(b.ymin);
        positions_cube.push_back(b.zmin);


        positions_cube.push_back(b.xmax);
        positions_cube.push_back(b.ymin);
        positions_cube.push_back(b.zmin);


        positions_cube.push_back(b.xmax);
        positions_cube.push_back(b.ymin);
        positions_cube.push_back(b.zmin);


        positions_cube.push_back(b.xmax);
        positions_cube.push_back(b.ymax);
        positions_cube.push_back(b.zmin);


        positions_cube.push_back(b.xmax);
        positions_cube.push_back(b.ymin);
        positions_cube.push_back(b.zmin);


        positions_cube.push_back(b.xmax);
        positions_cube.push_back(b.ymin);
        positions_cube.push_back(b.zmax);


        positions_cube.push_back(b.xmin);
        positions_cube.push_back(b.ymax);
        positions_cube.push_back(b.zmin);


        positions_cube.push_back(b.xmin);
        positions_cube.push_back(b.ymax);
        positions_cube.push_back(b.zmax);


        positions_cube.push_back(b.xmin);
        positions_cube.push_back(b.ymax);
        positions_cube.push_back(b.zmin);


        positions_cube.push_back(b.xmax);
        positions_cube.push_back(b.ymax);
        positions_cube.push_back(b.zmin);


        positions_cube.push_back(b.xmax);
        positions_cube.push_back(b.ymax);
        positions_cube.push_back(b.zmin);


        positions_cube.push_back(b.xmax);
        positions_cube.push_back(b.ymax);
        positions_cube.push_back(b.zmax);


        positions_cube.push_back(b.xmin);
        positions_cube.push_back(b.ymin);
        positions_cube.push_back(b.zmax);


        positions_cube.push_back(b.xmin);
        positions_cube.push_back(b.ymax);
        positions_cube.push_back(b.zmax);


        positions_cube.push_back(b.xmin);
        positions_cube.push_back(b.ymin);
        positions_cube.push_back(b.zmax);


        positions_cube.push_back(b.xmax);
        positions_cube.push_back(b.ymin);
        positions_cube.push_back(b.zmax);


        positions_cube.push_back(b.xmax);
        positions_cube.push_back(b.ymax);
        positions_cube.push_back(b.zmax);


        positions_cube.push_back(b.xmin);
        positions_cube.push_back(b.ymax);
        positions_cube.push_back(b.zmax);


        positions_cube.push_back(b.xmax);
        positions_cube.push_back(b.ymax);
        positions_cube.push_back(b.zmax);


        positions_cube.push_back(b.xmax);
        positions_cube.push_back(b.ymin);
        positions_cube.push_back(b.zmax);

    }

    //The texture
    for( int i=0 ; i < texture->getWidth() ; i++ )
    {
        for( int j=0 ; j < texture->getHeight() ; j++)
        {
            compute_texture(i,j);
        }
    }
}

Scene_implicit_function_item::
Scene_implicit_function_item(Implicit_function_interface* f)
    :Scene_item(4,3)
    , function_(f)
    , frame_(new ManipulatedFrame())
    , need_update_(true)
    , grid_size_(SCENE_IMPLICIT_GRID_SIZE)
    , max_value_(0.)
    , min_value_(0.)
    , blue_color_ramp_()
    , red_color_ramp_()
    , textureId(-1)
{
    texture = new Texture(grid_size_-1,grid_size_-1);
    blue_color_ramp_.build_blue();
    red_color_ramp_.build_red();
    //
    //Generates an integer which will be used as ID for each buffer

    compute_min_max();
    compute_function_grid();
    double offset_x = (bbox().xmin + bbox().xmax) / 2;
    double offset_y = (bbox().ymin + bbox().ymax) / 2;
    double offset_z = (bbox().zmin + bbox().zmax) / 2;
    frame_->setPosition(offset_x, offset_y, offset_z);
    frame_->setOrientation(1., 0, 0, 0);
    connect(frame_, SIGNAL(modified()), this, SLOT(plane_was_moved()));

    invalidate_buffers();
}


Scene_implicit_function_item::~Scene_implicit_function_item()
{

    delete frame_;

}


Scene_implicit_function_item::Bbox
Scene_implicit_function_item::bbox() const
{
    return function_->bbox();
}

void
Scene_implicit_function_item::draw(Viewer_interface* viewer) const
{
    if(!are_buffers_filled)
        initialize_buffers(viewer);

    if(frame_->isManipulated()) {
        if(need_update_) {
            compute_function_grid();
            need_update_ = false;
        }
    }
    vaos[0]->bind();
    viewer->glActiveTexture(GL_TEXTURE0);
    viewer->glBindTexture(GL_TEXTURE_2D, textureId);
    attrib_buffers(viewer, PROGRAM_WITH_TEXTURE);
    QMatrix4x4 f_mat;
    GLdouble d_mat[16];
    frame_->getMatrix(d_mat);
    //Convert the GLdoubles matrices in GLfloats
    for (int i=0; i<16; ++i){
        f_mat.data()[i] = GLfloat(d_mat[i]);
    }
    program = getShaderProgram(PROGRAM_WITH_TEXTURE);
    program->bind();
    program->setUniformValue("f_matrix", f_mat);
    program->setUniformValue("light_amb", QVector4D(1.0,1.0,1.0,1.0));
    program->setUniformValue("light_diff", QVector4D(0,0,0,1));
    program->setAttributeValue("color_facets", QVector3D(1.0,1.0,1.0));
    viewer->glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(positions_tex_quad.size()/3));
    vaos[0]->release();
    program->release();
}

void
Scene_implicit_function_item::draw_edges(Viewer_interface* viewer) const
{
    if(!are_buffers_filled)
        initialize_buffers(viewer);
    //  draw_aux(viewer, true);
    vaos[1]->bind();
    attrib_buffers(viewer, PROGRAM_WITHOUT_LIGHT);
    program = getShaderProgram(PROGRAM_WITHOUT_LIGHT);
    program->bind();
    viewer->glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(positions_cube.size()/3));
    vaos[1]->release();
    vaos[2]->bind();
    QMatrix4x4 f_mat;
    GLdouble d_mat[16];
    frame_->getMatrix(d_mat);
    //Convert the GLdoubles matrices in GLfloats
    for (int i=0; i<16; ++i){
        f_mat.data()[i] = GLfloat(d_mat[i]);
    }
    program->setUniformValue("f_matrix", f_mat);
    viewer->glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(positions_grid.size()/3));
    vaos[2]->release();
    program->release();
}

QString
Scene_implicit_function_item::toolTip() const
{
    return tr("<p>Function <b>%1</b>")
            .arg(this->name());
}

bool
Scene_implicit_function_item::supportsRenderingMode(RenderingMode m) const
{
    switch ( m )
    {
    case Splatting:
    case Gouraud:
        return false;

    case Points:
    case Wireframe:
    case Flat:
    case FlatPlusEdges:
        return true;

    default:
        return false;
    }

    return false;
}

void Scene_implicit_function_item::compute_texture(int i, int j)
{
    double v = (implicit_grid_[i][j]).second;

    if(is_nan(v)) {
        texture->setData(i,j,51,51,51);
    } else
        // determines grey level
        if ( v > 0 )
        {
            v = v/max_value_;
            GLdouble r = red_color_ramp_.r(v), g = red_color_ramp_.g(v), b = red_color_ramp_.b(v);
            texture->setData(i,j,255*r,255*g,255*b);
        }
        else
        {
            v = v/min_value_;
            GLdouble r = blue_color_ramp_.r(v), g = blue_color_ramp_.g(v), b = blue_color_ramp_.b(v);
            texture->setData(i,j,255*r,255*g,255*b);
        }
}


void
Scene_implicit_function_item::
compute_function_grid() const
{
    typedef CGAL::Simple_cartesian<double>  K;
    typedef K::Aff_transformation_3         Aff_transformation;
    typedef K::Point_3                      Point_3;

    // Get transformation
    const GLdouble* m = frame_->matrix();

    // OpenGL matrices are row-major matrices
    Aff_transformation t (m[0], m[4], m[8], m[12],
            m[1], m[5], m[9], m[13],
            m[2], m[6], m[10], m[14]);

    double diag = bbox().diagonal_length() * .6;

    const double dx = diag;
    const double dy = diag;
    const double z (0);

    int nb_quad = grid_size_ - 1;

    for(int i=0 ; i<grid_size_ ; ++i)
    {
        double x = -diag/2. + double(i)/double(nb_quad) * dx;

        for(int j=0 ; j<grid_size_ ; ++j)
        {
            double y = -diag/2. + double(j)/double(nb_quad) * dy;

            Point_3 query = t( Point_3(x, y, z) );
            double v = function_->operator()(query.x(), query.y(), query.z());

            implicit_grid_[i][j] = Point_value(Point(query.x(),query.y(),query.z()),v);
        }
    }

    // Update
    const_cast<Scene_implicit_function_item*>(this)->invalidate_buffers();

}

void
Scene_implicit_function_item::
compute_min_max()
{
    if(function_->get_min_max(min_value_, max_value_))
        return;

    double probes_nb = double(grid_size_) / 2;

    // Probe bounding box
    const Bbox& b = bbox();

    for ( int i = 0 ; i <= probes_nb ; ++i )
    {
        double x = b.xmin + double(i) * (b.xmax - b.xmin) / probes_nb;

        for ( int j = 0 ; j <= probes_nb ; ++j )
        {
            double y = b.ymin + double(j) * (b.ymax - b.ymin) / probes_nb;

            for ( int k = 0 ; k <= probes_nb ; ++k )
            {
                double z = b.zmin + double(k) * (b.zmax - b.zmin) / probes_nb;

                double v = (*function_)(x,y,z);
                if(is_nan(v)) continue;
                max_value_ = (std::max)(v,max_value_);
                min_value_ = (std::min)(v,min_value_);
            }
        }
    }
}

void
Scene_implicit_function_item::invalidate_buffers()
{
    Scene_item::invalidate_buffers();
    compute_vertices_and_texmap();
    are_buffers_filled = false;
}

void Scene_implicit_function_item::contextual_changed()
{
    if(!frame_->isManipulated()) {
        if(need_update_) {
            compute_function_grid();
            compute_vertices_and_texmap();
            need_update_ = false;
        }
    }
}


