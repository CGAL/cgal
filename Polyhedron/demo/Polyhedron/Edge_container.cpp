#include <CGAL/Three/Edge_container.h>


typedef CGAL::Three::Viewer_interface VI;
using namespace CGAL::Three;
Edge_container::Edge_container(int program, bool indexed)
  :Primitive_container(program, indexed)
{
  VBOs.resize(NbOfVbos);

}

void Edge_container::initGL(Scene_item *item, Viewer_interface *viewer) const
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
      VAOs[viewer] = new CGAL::Three::Vao(item->getShaderProgram(program_id, viewer));
      VAOs[viewer]->addVbo(VBOs[Vertices]);
      VAOs[viewer]->addVbo(VBOs[Indices]);
    }
      break;
    default:
      break;
    }
  }
  is_gl_init[viewer] = true;
}
void Edge_container::draw(const CGAL::Three::Scene_item& item,
                          CGAL::Three::Viewer_interface* viewer,
                          bool is_color_uniform) const
{
  if(!is_init[viewer])
  {
    initializeBuffers(viewer);
  }
  item.attribBuffers(viewer, program_id);

  if(indexed)
  {
    VAOs[viewer]->bind();
    if(item.isSelected())
      VAOs[viewer]->program->setAttributeValue("is_selected", true);
    else
      VAOs[viewer]->program->setAttributeValue("is_selected", false);
    if(is_color_uniform)
      VAOs[viewer]->program->setAttributeValue("colors", color);
    VBOs[Indices]->bind();
    viewer->glDrawElements(GL_LINES, static_cast<GLuint>(idx_size),
                   GL_UNSIGNED_INT, 0);
    VBOs[Indices]->release();
    VAOs[viewer]->release();
  }
}
