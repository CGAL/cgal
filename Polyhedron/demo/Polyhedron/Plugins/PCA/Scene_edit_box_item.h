#ifndef SCENE_EDIT_BOX_ITEM_H
#define SCENE_EDIT_BOX_ITEM_H

#include <CGAL/Three/Scene_item.h>
#include <CGAL/Simple_cartesian.h>
#include "create_sphere.h"
struct Scene_edit_box_item_priv;
class Q_DECL_EXPORT Scene_edit_box_item: public CGAL::Three::Scene_item
{
    Q_OBJECT
  public:
    typedef CGAL::Simple_cartesian<double>  Kernel;
    struct vertex;
    struct edge;
    struct triangle;
    struct face;
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
    Scene_edit_box_item(const CGAL::Three::Scene_interface* scene_interface);
    ~Scene_edit_box_item();
    bool isFinite() const { return true; }
    bool isEmpty() const { return true; }
    void compute_bbox() const;

    Scene_edit_box_item* clone() const {
      return 0;
    }

    QString toolTip() const;

    // Indicate if rendering mode is supported
    bool supportsRenderingMode(RenderingMode m) const {
      return (m == Wireframe);
    }

    void drawEdges(CGAL::Three::Viewer_interface* viewer) const;
    void invalidateOpenGLBuffers()
    {
      compute_bbox();
      are_buffers_filled = false;
    }

protected:
    friend struct Scene_edit_box_item_priv;
    Scene_edit_box_item_priv* d;
};
#endif // SCENE_EDIT_BOX_ITEM_H
