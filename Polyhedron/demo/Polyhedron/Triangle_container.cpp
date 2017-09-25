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
  viewer->makeCurrent();
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
      VAOs[viewer] = new CGAL::Three::Vao(item->getShaderProgram(program_id, viewer));
      VAOs[viewer]->addVbo(VBOs[Smooth_vertices]);
      VAOs[viewer]->addVbo(VBOs[Vertex_indices]);
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

    {
      if(!VBOs[Flat_vertices])
        VBOs[Flat_vertices] =
            new CGAL::Three::Vbo("vertex");
      if(!VBOs[Flat_normals])
        VBOs[Flat_normals] =
            new CGAL::Three::Vbo("normals");
      if(!VBOs[FColors])
        VBOs[FColors] =
            new CGAL::Three::Vbo("colors");
      VAOs[viewer] = new CGAL::Three::Vao(item->getShaderProgram(program_id, viewer));
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
  }
  is_gl_init[viewer] = true;
}

void Triangle_container::draw(const CGAL::Three::Scene_item& item,
                              CGAL::Three::Viewer_interface* viewer,
                              bool is_color_uniform ) const
{
//  if(!is_init[viewer])
//  {
//    initializeBuffers(viewer);
//  }
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
    VBOs[Vertex_indices]->bind();
    viewer->glDrawElements(GL_TRIANGLES, static_cast<unsigned int>(idx_size),
                           GL_UNSIGNED_INT, 0 );
    VBOs[Vertex_indices]->release();
    VAOs[viewer]->release();
  }
  else
  {
    VAOs[viewer]->bind();
    if(program_id == VI::PROGRAM_C3T3
       //|| program_id == VI::PROGRAM_C3T3_TETS
       )
      VAOs[viewer]->program->setUniformValue("shrink_factor", shrink_factor);
    if(program_id == VI::PROGRAM_C3T3)
      VAOs[viewer]->program->setUniformValue("cutplane", plane);
    if(item.isSelected())
      VAOs[viewer]->program->setAttributeValue("is_selected", true);
    else
      VAOs[viewer]->program->setAttributeValue("is_selected", false);
    if(is_color_uniform)
      VAOs[viewer]->program->setAttributeValue("colors", color);
    viewer->glDrawArrays(GL_TRIANGLES,0,static_cast<GLsizei>(flat_size/3));

    VAOs[viewer]->release();
  }
}
