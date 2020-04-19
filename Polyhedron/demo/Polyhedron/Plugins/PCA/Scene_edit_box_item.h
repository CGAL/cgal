#ifndef SCENE_EDIT_BOX_ITEM_H
#define SCENE_EDIT_BOX_ITEM_H

#include <CGAL/Three/Scene_item_rendering_helper.h>
#include <CGAL/Three/Scene_transparent_interface.h>
#include <CGAL/Simple_cartesian.h>
#include "create_sphere.h"
#include "Scene_edit_box_item_config.h"
struct Scene_edit_box_item_priv;
class SCENE_EDIT_BOX_ITEM_EXPORT Scene_edit_box_item:
    public CGAL::Three::Scene_item_rendering_helper,
    public CGAL::Three::Scene_transparent_interface
{
    Q_OBJECT
    Q_INTERFACES(CGAL::Three::Scene_transparent_interface)
    Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.TransparentInterface/1.0")
  public:
    typedef CGAL::Simple_cartesian<double>  Kernel;
    struct vertex;
    struct edge;
    struct face;
    Scene_edit_box_item();
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
    void drawTransparent(CGAL::Three::Viewer_interface*)const;
    void drawHl(CGAL::Three::Viewer_interface *) const;
    void drawEdges(CGAL::Three::Viewer_interface* viewer) const;
    void drawSpheres(CGAL::Three::Viewer_interface* viewer, const QMatrix4x4 f_matrix) const;
    void invalidateOpenGLBuffers();

    //      5-----6
    //  .   |  .  |
    // 4------7   |
    // |    | |   |
    // |    1-|---2
    // | .    |.
    // 0------3

    double point(short i, short j) const;
    void initializeBuffers(CGAL::Three::Viewer_interface *) const;
    void computeElements() const;
public Q_SLOTS:
    void highlight(CGAL::Three::Viewer_interface* viewer);
    void clearHL();
    void connectNewViewer(QObject* o);
protected:
    friend struct Scene_edit_box_item_priv;
    Scene_edit_box_item_priv* d;
};
#endif // SCENE_EDIT_BOX_ITEM_H
