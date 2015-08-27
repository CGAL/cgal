#include "config.h"

#include "Volume_plane_intersection.h"
#include "Volume_plane_interface.h"

#include <CGAL/gl.h>

void Volume_plane_intersection::compile_shaders()
{

    for(int i=0; i< vboSize; i++)
        buffers[i].create();
    for(int i=0; i< vaoSize; i++)
        vao[i].create();

    //Vertex source code
    const char vertex_source[] =
    {
        "#version 120 \n"
        "attribute highp vec4 vertex;\n"
        "uniform highp mat4 mvp_matrix;\n"
        "uniform highp mat4 f_matrix; \n"
        "void main(void)\n"
        "{\n"
        "   gl_Position = mvp_matrix* f_matrix * vertex; \n"
        "}"
    };
    //Fragment source code
    const char fragment_source[] =
    {
        "#version 120 \n"
        "void main(void) { \n"
        "gl_FragColor = vec4(1.0,0.0,0.0,1.0); \n"
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
}

void Volume_plane_intersection::compute_elements()
{
   a_vertex.resize(0);
   b_vertex.resize(0);
   c_vertex.resize(0);

   a_vertex.push_back(0.0); a_vertex.push_back(0.0); a_vertex.push_back(0.0);
   a_vertex.push_back(x);   a_vertex.push_back(0.0); a_vertex.push_back(0.0);

   b_vertex.push_back(0.0); b_vertex.push_back(0.0); b_vertex.push_back(0.0);
   b_vertex.push_back(0.0); b_vertex.push_back(y);   b_vertex.push_back(0.0);

   c_vertex.push_back(0.0); c_vertex.push_back(0.0); c_vertex.push_back(0.0);
   c_vertex.push_back(0.0); c_vertex.push_back(0.0); c_vertex.push_back(z);
}

void Volume_plane_intersection::init_buffers()
{
    rendering_program.bind();

        vao[0].bind();
        buffers[0].bind();
        buffers[0].allocate(a_vertex.data(), a_vertex.size()*sizeof(float));
        vertexLocation[0] = rendering_program.attributeLocation("vertex");
        rendering_program.enableAttributeArray(vertexLocation[0]);
        rendering_program.setAttributeBuffer(vertexLocation[0],GL_FLOAT,0,3);
        buffers[0].release();
        vao[0].release();

        vao[1].bind();
        buffers[1].bind();
        buffers[1].allocate(b_vertex.data(), b_vertex.size()*sizeof(float));
        vertexLocation[0] = rendering_program.attributeLocation("vertex");
        rendering_program.enableAttributeArray(vertexLocation[0]);
        rendering_program.setAttributeBuffer(vertexLocation[0],GL_FLOAT,0,3);
        buffers[1].release();
        vao[1].release();

        vao[2].bind();
        buffers[2].bind();
        buffers[2].allocate(c_vertex.data(), c_vertex.size()*sizeof(float));
        vertexLocation[0] = rendering_program.attributeLocation("vertex");
        rendering_program.enableAttributeArray(vertexLocation[0]);
        rendering_program.setAttributeBuffer(vertexLocation[0],GL_FLOAT,0,3);
        buffers[2].release();
        vao[2].release();



    rendering_program.release();

}

void Volume_plane_intersection::attrib_buffers(Viewer* viewer) const
{
    QMatrix4x4 mvpMatrix;
    double mat[16];
    viewer->camera()->getModelViewProjectionMatrix(mat);
    for(int i=0; i < 16; i++)
    {
        mvpMatrix.data()[i] = (float)mat[i];
    }

    rendering_program.bind();
    mvpLocation[0] = rendering_program.uniformLocation("mvp_matrix");
    rendering_program.setUniformValue(mvpLocation[0], mvpMatrix);
    rendering_program.release();
}

void Volume_plane_intersection::draw(Viewer* viewer) const {
  viewer->glLineWidth(4.0f);
  attrib_buffers(viewer);
  if(b && c) {

    vao[0].bind();
    rendering_program.bind();
    GLdouble mat[16];
    b->manipulatedFrame()->getMatrix(mat);
    QMatrix4x4 b_mat, c_mat;
    for(int i=0; i<16; i++)
    {
       b_mat.data()[i] = (float)mat[i];
    }
    c->manipulatedFrame()->getMatrix(mat);
    for(int i=0; i<16; i++)
    {
       c_mat.data()[i] = (float)mat[i];
    }
    rendering_program.setUniformValue("f_matrix", b_mat*c_mat);
    viewer->glDrawArrays(GL_LINES, 0, 2);
    rendering_program.release();
    vao[0].release();
  }

  if(a && c) {
      vao[1].bind();
      rendering_program.bind();
      GLdouble mat[16];
      a->manipulatedFrame()->getMatrix(mat);
      QMatrix4x4 a_mat, c_mat;
      for(int i=0; i<16; i++)
      {
         a_mat.data()[i] = (float)mat[i];
      }
      c->manipulatedFrame()->getMatrix(mat);
      for(int i=0; i<16; i++)
      {
         c_mat.data()[i] = (float)mat[i];
      }
      rendering_program.setUniformValue("f_matrix", a_mat*c_mat);
      viewer->glDrawArrays(GL_LINES, 0, 2);
      rendering_program.release();
      vao[1].release();
  }

  if(a && b) {
      vao[2].bind();
      rendering_program.bind();
      GLdouble mat[16];
      a->manipulatedFrame()->getMatrix(mat);
      QMatrix4x4 a_mat, b_mat;
      for(int i=0; i<16; i++)
      {
         a_mat.data()[i] = (float)mat[i];
      }
      b->manipulatedFrame()->getMatrix(mat);
      for(int i=0; i<16; i++)
      {
         b_mat.data()[i] = (float)mat[i];
      }
      rendering_program.setUniformValue("f_matrix", a_mat*b_mat);
      viewer->glDrawArrays(GL_LINES, 0, 2);
      rendering_program.release();
      vao[2].release();
  }

  viewer->glLineWidth(1.0f);
}
