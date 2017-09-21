#include <CGAL/Three/Edge_container.h>


typedef CGAL::Three::Viewer_interface VI;
using namespace CGAL::Three;
Edge_container::Edge_container(CGAL::Three::Scene_item* item, CGAL::Three::Viewer_interface* viewer, int program, bool indexed)
  :Primitive_container(program, indexed)
{
  VBOs.resize(NbOfVbos);
  if(indexed)
  {
    switch(program_id)
    {
    case VI::PROGRAM_WITHOUT_LIGHT:
    case VI::PROGRAM_NO_SELECTION:
    {
      VBOs[Vertices] =
          new CGAL::Three::Vbo("vertex");
      VBOs[Indices] =
          new CGAL::Three::Vbo("indices",
                               QOpenGLBuffer::IndexBuffer);
      VAO = new CGAL::Three::Vao(item->getShaderProgram(program_id, viewer));
      VAO->addVbo(VBOs[Vertices]);
      VAO->addVbo(VBOs[Indices]);
    }
      break;
    default:
      break;
    }
  }
}

void Edge_container::draw(const CGAL::Three::Scene_item& item,
                          CGAL::Three::Viewer_interface* viewer,
                          bool is_color_uniform) const
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
    VBOs[Indices]->bind();
    glDrawElements(GL_LINES, static_cast<GLuint>(idx_size),
                   GL_UNSIGNED_INT, 0);
    VBOs[Indices]->release();
    VAO->release();
  }
}
