#include "Scene_spheres_item.h"

void Scene_spheres_item::computeElements() const
{
  colors.clear();
  edges_colors.clear();
  centers.clear();
  radius.clear();
  Q_FOREACH(sphere_pair sp, spheres)
  {
    colors.push_back((float)sp.second.red()/255);
    colors.push_back((float)sp.second.green()/255);
    colors.push_back((float)sp.second.blue()/255);

    edges_colors.push_back((float)sp.second.red()/255);
    edges_colors.push_back((float)sp.second.green()/255);
    edges_colors.push_back((float)sp.second.blue()/255);

    centers.push_back(sp.first->center().x());
    centers.push_back(sp.first->center().y());
    centers.push_back(sp.first->center().z());

    radius.push_back(sp.first->squared_radius());

  }
}

void Scene_spheres_item::initializeBuffers(CGAL::Three::Viewer_interface *viewer) const
{
  if(has_plane)
  {
    program = getShaderProgram(PROGRAM_CUTPLANE_SPHERES, viewer);
    attrib_buffers(viewer, PROGRAM_CUTPLANE_SPHERES);
  }
  else
  {
    program = getShaderProgram(PROGRAM_SPHERES, viewer);
    attrib_buffers(viewer, PROGRAM_SPHERES);
  }

  program->bind();
  vaos[Facets]->bind();
  buffers[Vertices].bind();
  buffers[Vertices].allocate(vertices.data(),
                             static_cast<int>(vertices.size()*sizeof(float)));
  program->enableAttributeArray("vertex");
  program->setAttributeBuffer("vertex", GL_FLOAT, 0, 3);
  buffers[Vertices].release();

  buffers[Normals].bind();
  buffers[Normals].allocate(normals.data(),
                            static_cast<int>(normals.size()*sizeof(float)));
  program->enableAttributeArray("normals");
  program->setAttributeBuffer("normals", GL_FLOAT, 0, 3);
  buffers[Normals].release();

  buffers[Color].bind();
  buffers[Color].allocate(colors.data(),
                          static_cast<int>(colors.size()*sizeof(float)));
  program->enableAttributeArray("colors");
  program->setAttributeBuffer("colors", GL_FLOAT, 0, 3);
  buffers[Color].release();

  buffers[Radius].bind();
  buffers[Radius].allocate(radius.data(),
                           static_cast<int>(radius.size()*sizeof(float)));
  program->enableAttributeArray("radius");
  program->setAttributeBuffer("radius", GL_FLOAT, 0, 1);
  buffers[Radius].release();

  buffers[Center].bind();
  buffers[Center].allocate(centers.data(),
                           static_cast<int>(centers.size()*sizeof(float)));
  program->enableAttributeArray("center");
  program->setAttributeBuffer("center", GL_FLOAT, 0, 3);
  buffers[Center].release();

  viewer->glVertexAttribDivisor(program->attributeLocation("center"), 1);
  viewer->glVertexAttribDivisor(program->attributeLocation("radius"), 1);
  viewer->glVertexAttribDivisor(program->attributeLocation("colors"), 1);
  vaos[Facets]->release();


  vaos[Edges]->bind();
  buffers[Edge_vertices].bind();
  buffers[Edge_vertices].allocate(edges.data(),
                                  static_cast<int>(edges.size()*sizeof(float)));
  program->enableAttributeArray("vertex");
  program->setAttributeBuffer("vertex", GL_FLOAT, 0, 3);
  buffers[Edge_vertices].release();

  buffers[Normals].bind();
  program->enableAttributeArray("normals");
  program->setAttributeBuffer("normals", GL_FLOAT, 0, 3);
  buffers[Normals].release();

  buffers[Edge_color].bind();
  buffers[Edge_color].allocate(edges_colors.data(),
                               static_cast<int>(edges_colors.size()*sizeof(float)));
  program->enableAttributeArray("colors");
  program->setAttributeBuffer("colors", GL_FLOAT, 0, 3);
  buffers[Edge_color].release();

  buffers[Radius].bind();
  program->enableAttributeArray("radius");
  program->setAttributeBuffer("radius", GL_FLOAT, 0, 1);
  buffers[Radius].release();

  buffers[Center].bind();
  program->enableAttributeArray("center");
  program->setAttributeBuffer("center", GL_FLOAT, 0, 3);
  buffers[Center].release();

  viewer->glVertexAttribDivisor(program->attributeLocation("center"), 1);
  viewer->glVertexAttribDivisor(program->attributeLocation("radius"), 1);
  viewer->glVertexAttribDivisor(program->attributeLocation("colors"), 1);
  vaos[Edges]->release();

  program->release();

  nb_centers = centers.size();
  centers.clear();
  centers.swap(centers);
  colors.clear();
  colors.swap(colors);
  radius.clear();
  radius.swap(radius);
  edges_colors.clear();
  edges_colors.swap(edges_colors);

  are_buffers_filled = true;
}

void Scene_spheres_item::draw(Viewer_interface *viewer) const
{
  if (!are_buffers_filled)
  {
    computeElements();
    initializeBuffers(viewer);
  }
  vaos[Facets]->bind();
  if(has_plane)
  {
    program = getShaderProgram(PROGRAM_CUTPLANE_SPHERES, viewer);
    attrib_buffers(viewer, PROGRAM_CUTPLANE_SPHERES);
    program->bind();
    QVector4D cp(plane.a(),plane.b(),plane.c(),plane.d());
    program->setUniformValue("cutplane", cp);

  }
  else
  {
    program = getShaderProgram(PROGRAM_SPHERES, viewer);
    attrib_buffers(viewer, PROGRAM_SPHERES);
    program->bind();
  }
  viewer->glDrawArraysInstanced(GL_TRIANGLES, 0,
                                static_cast<GLsizei>(vertices.size()/3),
                                static_cast<GLsizei>(nb_centers));
  program->release();
  vaos[Facets]->release();
}
void Scene_spheres_item::draw_edges(Viewer_interface *viewer) const
{
  if (!are_buffers_filled)
  {
    computeElements();
    initializeBuffers(viewer);
  }
  vaos[Edges]->bind();
  if(has_plane)
  {
    program = getShaderProgram(PROGRAM_CUTPLANE_SPHERES, viewer);
    attrib_buffers(viewer, PROGRAM_CUTPLANE_SPHERES);
    program->bind();
    QVector4D cp(plane.a(),plane.b(),plane.c(),plane.d());
    program->setUniformValue("cutplane", cp);
  }
  else
  {
    program = getShaderProgram(PROGRAM_SPHERES, viewer);
    attrib_buffers(viewer, PROGRAM_SPHERES);
    program->bind();
  }
  viewer->glDrawArraysInstanced(GL_LINES, 0,
                                static_cast<GLsizei>(edges.size()/3),
                                static_cast<GLsizei>(nb_centers));
  program->release();
  vaos[Edges]->release();
}
void Scene_spheres_item::add_sphere(CGAL::Sphere_3<Kernel> *sphere, CGAL::Color color)
{
  sphere_pair pair_(sphere, color);
  spheres.append(pair_);
}

void Scene_spheres_item::remove_sphere(CGAL::Sphere_3<Kernel> *sphere)
{
  Q_FOREACH(sphere_pair pair_, spheres)
    if(pair_.first == sphere)
    {
      spheres.removeAll(pair_);
      break;
    }
}
