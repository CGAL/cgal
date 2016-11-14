#ifndef SCENE_EDIT_BOX_ITEM_H
#define SCENE_EDIT_BOX_ITEM_H

#include <CGAL/Three/Scene_item.h>
#include <CGAL/Simple_cartesian.h>
#include "create_sphere.h"
#include "Scene_edit_box_item_config.h"
struct Scene_edit_box_item_priv;
class SCENE_EDIT_BOX_ITEM_EXPORT Scene_edit_box_item: public CGAL::Three::Scene_item
{
    Q_OBJECT
  public:
    typedef CGAL::Simple_cartesian<double>  Kernel;
    struct vertex;
    struct edge;
    struct face;
    Scene_edit_box_item(const CGAL::Three::Scene_interface* scene_interface);
    ~Scene_edit_box_item();
    bool isFinite() const { return true; }
    bool isEmpty() const { return false; }
    void compute_bbox() const;

    bool manipulatable() const { return true; }
    ManipulatedFrame* manipulatedFrame();
    Scene_edit_box_item* clone() const {
      return 0;
    }

    QString toolTip() const;

    bool eventFilter(QObject *, QEvent *);
    // Indicate if rendering mode is supported
    bool supportsRenderingMode(RenderingMode m) const;
    void draw(CGAL::Three::Viewer_interface *) const;
    void drawHl(CGAL::Three::Viewer_interface *) const;
    void drawEdges(CGAL::Three::Viewer_interface* viewer) const;
    void drawSpheres(CGAL::Three::Viewer_interface* viewer, const QMatrix4x4 f_matrix) const;
    void invalidateOpenGLBuffers()
    {
      compute_bbox();
      are_buffers_filled = false;
    }
    double point(short i, short j) const;

public Q_SLOTS:
    void highlight();
    void clearHL();
protected:
    friend struct Scene_edit_box_item_priv;
    Scene_edit_box_item_priv* d;
};
#endif // SCENE_EDIT_BOX_ITEM_H
