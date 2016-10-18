#include "Scene_edit_box_item.h"
#include <QApplication>
#include <CGAL/Three/Viewer_interface.h>

struct Scene_edit_box_item::vertex{
  short id;
  CGAL::Point_3<Kernel> position;
  //  face& face_;
};
struct Scene_edit_box_item::edge{
  short id;
  vertex source;
  vertex target;
  //face& face_;
};
struct Scene_edit_box_item::triangle{
  short id;
  vertex vertices[3];
  //face& face_;
};
struct Scene_edit_box_item::face{
  short id;
  triangle triangles[2];
  edge edges[4];
  vertex vertices[4];
};

struct Scene_edit_box_item_priv{
  typedef CGAL::Simple_cartesian<double>  Kernel;
  enum VAOs{
    Edges = 0,
    Spheres,
    S_Edges,
    S_Spheres,
    Arrow,
    NumberOfVaos
  };
  enum VBOs{
    VertexEdges = 0,
    VertexSpheres,
    NormalSpheres,
    VertexArrow,
    NormalArrow,
    NumberOfVbos
  };

  Scene_edit_box_item_priv(const CGAL::Three::Scene_interface *scene_interface, Scene_edit_box_item* ebi)
  {
    scene = scene_interface;
    item = ebi;
    //      5-----6
    //  .   |  .  |
    // 4------7   |
    // |    | |   |
    // |    1-|---2
    // | .    |.
    // 0------3

    //vertices
    for(short i = 0; i< 8; ++i)
    {
      double x,y,z;
      x = ((i/2)%2)? scene->bbox().min(0):scene->bbox().max(0);
      y = (i/4)? scene->bbox().min(1):scene->bbox().max(1);
      z = (((i+1)/2)%2)? scene->bbox().min(2):scene->bbox().max(2);
      vertices[i].position = Kernel::Point_3(x,y,z);
      vertices[i].id = i;
    }

    //      .--5--.
    //  4   |  6  |
    // .--7-1-.   2
    // |    | |   |
    // 0    .-39--.
    // | 8    |10
    // .--11--.
    //edges
    for(short i=0; i<12; ++i)
    {
      edges[i].id = i;
      if(i<4)
      {
        edges[i].source = vertices[i];
        edges[i].target = vertices[i+4];
      }
      else if(i<8)
      {
        edges[i].source = vertices[i];
        edges[i].target = vertices[(i+1)%4 +4];
      }
      else
      {
        edges[i].source = vertices[i%4];
        edges[i].target = vertices[(i+1) %4];
      }

      vertex_edges.resize(0);
    }


    mutable std::vector<float> vertex_edges;
    mutable std::vector<float> vertex_spheres;
    mutable std::vector<float> normal_spheres;
    mutable std::vector<float> vertex_arrow;
    mutable std::vector<float> normal_arrow;

    mutable Scene_edit_box_item::vertex vertices[8];
    mutable Scene_edit_box_item::edge edges[12];
    mutable Scene_edit_box_item::triangle triangles[12];
    mutable Scene_edit_box_item::face faces[6];

    mutable QOpenGLShaderProgram *program;
    void initializeBuffers(CGAL::Three::Viewer_interface *viewer)const;

    void computeElements() const;

    const CGAL::Three::Scene_interface* scene;
    Scene_edit_box_item* item;
  };


  Scene_edit_box_item::Scene_edit_box_item(const CGAL::Three::Scene_interface *scene_interface)
    :  Scene_item(NumberOfVbos,NumberOfVaos)

  {
    d = new Scene_edit_box_item_priv(scene_interface, this);

    are_buffers_filled = false;
  }
  QString Scene_edit_box_item::toolTip() const {

    return QString();
  }

  void Scene_edit_box_item::drawEdges(CGAL::Three::Viewer_interface* viewer) const
  {
    if(!are_buffers_filled)
    {
      d->computeElements();
      d->initializeBuffers(viewer);
    }
    vaos[Edges]->bind();
    d->program = getShaderProgram(PROGRAM_WITHOUT_LIGHT);
    attribBuffers(viewer, PROGRAM_WITHOUT_LIGHT);
    d->program->bind();
    d->program->setAttributeValue("colors", this->color());
    viewer->glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(d->vertex_edges.size()/3));
    vaos[Edges]->release();
    d->program->release();

  }

  void Scene_edit_box_item::compute_bbox() const
  {

  }


  void Scene_edit_box_item_priv::initializeBuffers(CGAL::Three::Viewer_interface *viewer)const
  {

    //vao containing the data for the lines
    {
      program = item->getShaderProgram(Scene_edit_box_item::PROGRAM_WITHOUT_LIGHT, viewer);
      program->bind();

      item->vaos[Edges]->bind();
      item->buffers[VertexEdges].bind();
      item->buffers[VertexEdges].allocate(vertex_edges.data(),
                                          static_cast<GLsizei>(vertex_edges.size()*sizeof(float)));
      program->enableAttributeArray("vertex");
      program->setAttributeBuffer("vertex",GL_FLOAT,0,3);
      item->buffers[VertexEdges].release();

      item->vaos[Edges]->release();
      program->release();

    }
    item->are_buffers_filled = true;
  }

  void Scene_edit_box_item_priv::computeElements() const
  {
    QApplication::setOverrideCursor(Qt::WaitCursor);
    //edges
    for(short i=0; i<12; ++i)
    {
      if(i<4)
      {
        vertex_edges.push_back(vertices[i].position.x());
        vertex_edges.push_back(vertices[i].position.y());
        vertex_edges.push_back(vertices[i].position.z());

        vertex_edges.push_back(vertices[i+4].position.x());
        vertex_edges.push_back(vertices[i+4].position.y());
        vertex_edges.push_back(vertices[i+4].position.z());
      }
      else if(i<8)
      {
        vertex_edges.push_back(vertices[i].position.x());
        vertex_edges.push_back(vertices[i].position.y());
        vertex_edges.push_back(vertices[i].position.z());

        vertex_edges.push_back(vertices[(i+1)%4 +4].position.x());
        vertex_edges.push_back(vertices[(i+1)%4 +4].position.y());
        vertex_edges.push_back(vertices[(i+1)%4 +4].position.z());
      }
      else
      {
        vertex_edges.push_back(vertices[i%4].position.x());
        vertex_edges.push_back(vertices[i%4].position.y());
        vertex_edges.push_back(vertices[i%4].position.z());

        vertex_edges.push_back(vertices[(i+1) %4].position.x());
        vertex_edges.push_back(vertices[(i+1) %4].position.y());
        vertex_edges.push_back(vertices[(i+1) %4].position.z());
      }
    }
    QApplication::restoreOverrideCursor();
  }

  Scene_edit_box_item::~Scene_edit_box_item()
  {
    delete d;
  }
