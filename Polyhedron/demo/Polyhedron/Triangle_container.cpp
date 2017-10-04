#include <CGAL/Three/Triangle_container.h>
#include <QOpenGLFramebufferObject>

typedef CGAL::Three::Viewer_interface VI;
using namespace CGAL::Three;
Triangle_container::Triangle_container(int program, bool indexed)
  : Primitive_container(program, indexed)
{
  VBOs.resize(NbOfVbos);
}

void Triangle_container::initGL(const CGAL::Three::Scene_item& item, CGAL::Three::Viewer_interface* viewer)const
{
  viewer->makeCurrent();
  if(indexed)
  {
    switch(program_id)
    {
    case VI::PROGRAM_WITHOUT_LIGHT:
    case VI::PROGRAM_WITH_LIGHT:
    {
      if(!VBOs[Smooth_vertices])
        VBOs[Smooth_vertices] =
            new CGAL::Three::Vbo("vertex");
      if(!VBOs[Vertex_indices])
        VBOs[Vertex_indices] =
            new CGAL::Three::Vbo("indices",
                                 QOpenGLBuffer::IndexBuffer);
      if(!VAOs[viewer])
        VAOs[viewer] = new CGAL::Three::Vao(item.getShaderProgram(program_id, viewer));
      VAOs[viewer]->addVbo(VBOs[Smooth_vertices]);
      VAOs[viewer]->addVbo(VBOs[Vertex_indices]);
    }
      break;
    default:
      Q_UNUSED(viewer);
      break;
    }
    switch(program_id)
    {
    case VI::PROGRAM_WITH_LIGHT:
    {
      if(!VBOs[Smooth_normals])
        VBOs[Smooth_normals] =
            new CGAL::Three::Vbo("normals");
      if(!VBOs[VColors])
        VBOs[VColors] =
            new CGAL::Three::Vbo("colors", QOpenGLBuffer::VertexBuffer, GL_FLOAT, 0, 4);
      VAOs[viewer]->addVbo(VBOs[Smooth_normals]);
      VAOs[viewer]->addVbo(VBOs[VColors]);
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

    case VI::PROGRAM_WITH_LIGHT:
    case VI::PROGRAM_C3T3:
    //case VI::PROGRAM_C3T3_TETS:
    case VI::PROGRAM_SPHERES:
    case VI::PROGRAM_CUTPLANE_SPHERES:
    {
      if(!VBOs[Flat_vertices])
      {
        VBOs[Flat_vertices] =
            new CGAL::Three::Vbo("vertex");
      }
      if(!VBOs[Flat_normals])
        VBOs[Flat_normals] =
            new CGAL::Three::Vbo("normals");
      if(!VBOs[FColors])
        VBOs[FColors] =
            new CGAL::Three::Vbo("colors", QOpenGLBuffer::VertexBuffer, GL_FLOAT, 0, 4);
      VAOs[viewer] = new CGAL::Three::Vao(item.getShaderProgram(program_id, viewer));
      VAOs[viewer]->addVbo(VBOs[Flat_vertices]);
      VAOs[viewer]->addVbo(VBOs[Flat_normals]);
      VAOs[viewer]->addVbo(VBOs[FColors]);

    }
      break;
    default:
      Q_UNUSED(viewer);
      break;
    }
    switch(program_id)
    {
    case VI::PROGRAM_C3T3:
      //case VI::PROGRAM_C3T3_TETS:
    case VI::PROGRAM_SPHERES:
    case VI::PROGRAM_CUTPLANE_SPHERES:
    {
      if(!VBOs[Facet_barycenters])
        VBOs[Facet_barycenters] =
            new CGAL::Three::Vbo("barycenter");
      VAOs[viewer]->addVbo(VBOs[Facet_barycenters]);
    }
      break;
    default:
      Q_UNUSED(viewer);
      break;
    }
    if(program_id == VI::PROGRAM_SPHERES
       || program_id == VI::PROGRAM_CUTPLANE_SPHERES)
    {
      if(!VBOs[Radius])
        VBOs[Radius] =
            new CGAL::Three::Vbo("radius",
                                 QOpenGLBuffer::VertexBuffer, GL_FLOAT, 0, 1);
      VAOs[viewer]->addVbo(VBOs[Radius]);
      viewer->glVertexAttribDivisor(VAOs[viewer]->program->attributeLocation("barycenter"), 1);
      viewer->glVertexAttribDivisor(VAOs[viewer]->program->attributeLocation("radius"), 1);
      viewer->glVertexAttribDivisor(VAOs[viewer]->program->attributeLocation("colors"), 1);
    }
  }
  is_gl_init[viewer] = true;
}

void Triangle_container::draw(const CGAL::Three::Scene_item& item,
                              CGAL::Three::Viewer_interface* viewer,
                              bool is_color_uniform,
                              QOpenGLFramebufferObject* fbo) const
{

  item.attribBuffers(viewer, program_id);

  if(indexed)
  {
    VAOs[viewer]->bind();
    VAOs[viewer]->program->setAttributeValue("is_selected", is_selected);
    if(is_color_uniform)
      VAOs[viewer]->program->setAttributeValue("colors", color);
    VBOs[Vertex_indices]->bind();
    if(program_id == VI::PROGRAM_WITH_LIGHT)
    {
      VAOs[viewer]->program->setUniformValue("comparing", comparing);
      VAOs[viewer]->program->setUniformValue("width", width);
      VAOs[viewer]->program->setUniformValue("height", height);
      VAOs[viewer]->program->setUniformValue("near", near);
      VAOs[viewer]->program->setUniformValue("far", far);
      VAOs[viewer]->program->setUniformValue("writing", writing);
      if( fbo)
        viewer->glBindTexture(GL_TEXTURE_2D, fbo->texture());
    }
    viewer->glDrawElements(GL_TRIANGLES, static_cast<unsigned int>(idx_size),
                           GL_UNSIGNED_INT, 0 );
    VBOs[Vertex_indices]->release();
    VAOs[viewer]->release();
  }
  else
  {
    VAOs[viewer]->bind();
    if(program_id == VI::PROGRAM_C3T3)
      VAOs[viewer]->program->setUniformValue("shrink_factor", shrink_factor);
    if(program_id == VI::PROGRAM_C3T3
       || program_id == VI::PROGRAM_CUTPLANE_SPHERES)
      VAOs[viewer]->program->setUniformValue("cutplane", plane);
    VAOs[viewer]->program->setAttributeValue("is_selected", is_selected);
    if(is_color_uniform)
      VAOs[viewer]->program->setAttributeValue("colors", color);
    if(program_id == VI::PROGRAM_SPHERES
       || program_id == VI::PROGRAM_CUTPLANE_SPHERES)
    {
      viewer->glDrawArraysInstanced(GL_TRIANGLES, 0,
                                    static_cast<GLsizei>(flat_size/3),
                                    static_cast<GLsizei>(center_size/3));
    }
    else
    {
      if(program_id == VI::PROGRAM_WITH_LIGHT)
      {
        VAOs[viewer]->program->setUniformValue("comparing", comparing);
        VAOs[viewer]->program->setUniformValue("width", width);
        VAOs[viewer]->program->setUniformValue("height", height);
        VAOs[viewer]->program->setUniformValue("near", near);
        VAOs[viewer]->program->setUniformValue("far", far);
        VAOs[viewer]->program->setUniformValue("writing", writing);
        if( fbo)
          viewer->glBindTexture(GL_TEXTURE_2D, fbo->texture());
      }
      viewer->glDrawArrays(GL_TRIANGLES,0,static_cast<GLsizei>(flat_size/3));
    }

    VAOs[viewer]->release();
  }
}
