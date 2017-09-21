#include <CGAL/Three/Triangle_container.h>

typedef CGAL::Three::Viewer_interface VI;
using namespace CGAL::Three;
Triangle_container::Triangle_container(int program, bool indexed)
  : Primitive_container(program, indexed)
{
  VBOs.resize(NbOfVbos);
}

void Triangle_container::initGL(CGAL::Three::Scene_item* item, CGAL::Three::Viewer_interface* viewer)const
{
  if(indexed)
  {
    switch(program_id)
    {
    case VI::PROGRAM_WITH_LIGHT:
    {
      VBOs[Smooth_vertices] =
          new CGAL::Three::Vbo("vertex");
      VBOs[Smooth_normals] =
          new CGAL::Three::Vbo("normals");
      VBOs[VColors] =
          new CGAL::Three::Vbo("colors");
      VBOs[Vertex_indices] =
          new CGAL::Three::Vbo("indices",
                               QOpenGLBuffer::IndexBuffer);
      VAO = new CGAL::Three::Vao(item->getShaderProgram(program_id, viewer));
      VAO->addVbo(VBOs[Smooth_vertices]);
      VAO->addVbo(VBOs[Vertex_indices]);
      VAO->addVbo(VBOs[Smooth_normals]);
      VAO->addVbo(VBOs[VColors]);
    }
      break;
    default:
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

    {
      VBOs[Flat_vertices] =
          new CGAL::Three::Vbo("vertex");
      VBOs[Flat_normals] =
          new CGAL::Three::Vbo("normals");
      VBOs[FColors] =
          new CGAL::Three::Vbo("colors");
      VAO = new CGAL::Three::Vao(item->getShaderProgram(program_id, viewer));
      VAO->addVbo(VBOs[Flat_vertices]);
      VAO->addVbo(VBOs[Flat_normals]);
      VAO->addVbo(VBOs[FColors]);

    }
      break;
    }
    switch(program_id)
    {
    case VI::PROGRAM_C3T3:
      //case VI::PROGRAM_C3T3_TETS:
    {
      VBOs[Facet_barycenters] =
          new CGAL::Three::Vbo("barycenter");
      VAO->addVbo(VBOs[Facet_barycenters]);
    }
      break;
    default:
      break;
    }
  }
  is_gl_init = true;
}
void Triangle_container::draw(const CGAL::Three::Scene_item& item,
                              CGAL::Three::Viewer_interface* viewer,
                              bool is_color_uniform ) const
{
  if(!is_init)
  {
    initializeBuffers();
  }
  item.attribBuffers(viewer, program_id);

  if(indexed)
  {
    VAO->bind();
    if(item.isSelected())
      VAO->program->setAttributeValue("is_selected", true);
    else
      VAO->program->setAttributeValue("is_selected", false);
    if(is_color_uniform)
      VAO->program->setAttributeValue("colors", color);
    VBOs[Vertex_indices]->bind();
    viewer->glDrawElements(GL_TRIANGLES, static_cast<unsigned int>(idx_size),
                           GL_UNSIGNED_INT, 0 );
    VBOs[Vertex_indices]->release();
    VAO->release();
  }
  else
  {
    VAO->bind();
    if(program_id == VI::PROGRAM_C3T3
       //|| program_id == VI::PROGRAM_C3T3_TETS
       )
      VAO->program->setUniformValue("shrink_factor", shrink_factor);
    if(program_id == VI::PROGRAM_C3T3)
      VAO->program->setUniformValue("cutplane", plane);
    if(item.isSelected())
      VAO->program->setAttributeValue("is_selected", true);
    else
      VAO->program->setAttributeValue("is_selected", false);
    if(is_color_uniform)
      VAO->program->setAttributeValue("colors", color);
    viewer->glDrawArrays(GL_TRIANGLES,0,static_cast<GLsizei>(flat_size/3));

    VAO->release();
  }
}
