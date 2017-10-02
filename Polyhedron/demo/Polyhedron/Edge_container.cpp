#include <CGAL/Three/Edge_container.h>


typedef CGAL::Three::Viewer_interface VI;
using namespace CGAL::Three;
Edge_container::Edge_container(int program, bool indexed)
  :Primitive_container(program, indexed)
{
  VBOs.resize(NbOfVbos);

}

void Edge_container::initGL(const Scene_item& item, Viewer_interface *viewer) const
{
  viewer->makeCurrent();
  if(indexed)
  {
    switch(program_id)
    {
    case VI::PROGRAM_WITHOUT_LIGHT:
    case VI::PROGRAM_NO_SELECTION:
    {
      if(!VBOs[Vertices])
        VBOs[Vertices] =
            new CGAL::Three::Vbo("vertex");
      if(!VBOs[Indices])
        VBOs[Indices] =
            new CGAL::Three::Vbo("indices",
                                 QOpenGLBuffer::IndexBuffer);
      VAOs[viewer] = new CGAL::Three::Vao(item.getShaderProgram(program_id, viewer));
      VAOs[viewer]->addVbo(VBOs[Vertices]);
      VAOs[viewer]->addVbo(VBOs[Indices]);
    }
      break;
    default:
      Q_UNUSED(viewer);
      break;
    }
  }
  else
  {
    switch(program_id)
    {
    case VI::PROGRAM_SPHERES:
    case VI::PROGRAM_CUTPLANE_SPHERES:
    case VI::PROGRAM_NO_SELECTION:
    case VI::PROGRAM_WITHOUT_LIGHT:
    case VI::PROGRAM_C3T3_EDGES:
    {
      if(!VBOs[Vertices])
        VBOs[Vertices] =
            new CGAL::Three::Vbo("vertex");
      if(!VBOs[Colors])
        VBOs[Colors] =
            new CGAL::Three::Vbo("colors");
      VAOs[viewer] = new CGAL::Three::Vao(item.getShaderProgram(program_id, viewer));
      VAOs[viewer]->addVbo(VBOs[Vertices]);
      VAOs[viewer]->addVbo(VBOs[Colors]);
    }
      break;
    default:
      Q_UNUSED(viewer);
      break;
    }

    switch(program_id)
    {
    case VI::PROGRAM_SPHERES:
    case VI::PROGRAM_CUTPLANE_SPHERES:
    {
      if(!VBOs[Normals])
        VBOs[Normals] =
            new CGAL::Three::Vbo("normals");
      if(!VBOs[Radius])
        VBOs[Radius] =
            new CGAL::Three::Vbo("radius",
                                 QOpenGLBuffer::VertexBuffer, GL_FLOAT, 0, 1);
      if(!VBOs[Barycenters])
        VBOs[Barycenters] =
            new CGAL::Three::Vbo("barycenter");
      VAOs[viewer]->addVbo(VBOs[Normals]);
      VAOs[viewer]->addVbo(VBOs[Radius]);
      VAOs[viewer]->addVbo(VBOs[Barycenters]);
      viewer->glVertexAttribDivisor(VAOs[viewer]->program->attributeLocation("barycenter"), 1);
      viewer->glVertexAttribDivisor(VAOs[viewer]->program->attributeLocation("radius"), 1);
      viewer->glVertexAttribDivisor(VAOs[viewer]->program->attributeLocation("colors"), 1);
    }
      break;
    default:
      break;
    }
  }
  is_gl_init[viewer] = true;
}
void Edge_container::draw(const Scene_item &item, Viewer_interface *viewer,
                          bool is_color_uniform, QOpenGLFramebufferObject *) const

{
  if(!is_init[viewer])
  {
    initializeBuffers(viewer);
  }
  item.attribBuffers(viewer, program_id);

  if(indexed)
  {
    VAOs[viewer]->bind();
    VAOs[viewer]->program->setAttributeValue("is_selected", is_selected);
    if(is_color_uniform)
      VAOs[viewer]->program->setAttributeValue("colors", color);
    VBOs[Indices]->bind();
    viewer->glDrawElements(GL_LINES, static_cast<GLuint>(idx_size),
                           GL_UNSIGNED_INT, 0);
    VBOs[Indices]->release();
    VAOs[viewer]->release();
  }
  else
  {
    VAOs[viewer]->bind();
    VAOs[viewer]->program->setAttributeValue("is_selected", is_selected);
    if(is_color_uniform)
      VAOs[viewer]->program->setAttributeValue("colors", color);
    if(program_id == VI::PROGRAM_WITHOUT_LIGHT)
      VAOs[viewer]->program->setUniformValue("f_matrix", f_matrix);
    if(program_id == VI::PROGRAM_C3T3_EDGES
       || program_id == VI::PROGRAM_CUTPLANE_SPHERES)
      VAOs[viewer]->program->setUniformValue("cutplane", plane);
    if(program_id == VI::PROGRAM_SPHERES
       || program_id == VI::PROGRAM_CUTPLANE_SPHERES)
    {
      viewer->glDrawArraysInstanced(GL_LINES, 0,
                                    static_cast<GLsizei>(flat_size/3),
                                    static_cast<GLsizei>(center_size/3));
    }
    else
    {
      viewer->glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(flat_size/3));
    }
    VAOs[viewer]->release();

  }
}
