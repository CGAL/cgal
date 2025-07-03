#ifndef VARIATIONAL_MEDIAL_AXIS_SKELETON_ITEM_H
#define VARIATIONAL_MEDIAL_AXIS_SKELETON_ITEM_H

#include <CGAL/Three/Scene_item_rendering_helper.h>
#include <CGAL/Simple_cartesian.h>
#include "Scene_surface_mesh_item.h"
#include "create_sphere.h"
#include "Variational_medial_axis_skeleton_item_config.h"

struct Variational_medial_axis_skeleton_item_priv;
class VARIATIONAL_MEDIAL_AXIS_SKELETON_ITEM_EXPORT Variational_medial_axis_skeleton_item:
    public CGAL::Three::Scene_item_rendering_helper
{
    Q_OBJECT
  public:
    typedef CGAL::Simple_cartesian<double>  Kernel;
    struct vertex;
    struct edge;
    struct face;

    Variational_medial_axis_skeleton_item(const CGAL::Three::Scene_interface* scene_interface,
                                          const Scene_surface_mesh_item* sm_item,
                                          std::size_t nb_pts);
    ~Variational_medial_axis_skeleton_item();
    bool isFinite() const { return true; }
    bool isEmpty() const { return false; }
    void compute_bbox() const;

    bool manipulatable() const { return true; }
    ManipulatedFrame* manipulatedFrame();
    Variational_medial_axis_skeleton_item* clone() const {
      return nullptr;
    }

    QString toolTip() const;
    QMenu* contextMenu();

    bool eventFilter(QObject *, QEvent *);
    // Indicate if rendering mode is supported
    bool supportsRenderingMode(RenderingMode m) const;
    void draw(CGAL::Three::Viewer_interface *) const;
    void drawHl(CGAL::Three::Viewer_interface *) const;
    //~ void drawEdges(CGAL::Three::Viewer_interface* viewer) const;
    void drawSpheres(CGAL::Three::Viewer_interface* viewer, const QMatrix4x4 f_matrix) const;
//    void drawPath() const;
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
    friend struct Variational_medial_axis_skeleton_item_priv;
    Variational_medial_axis_skeleton_item_priv* d;
};
#endif // VARIATIONAL_MEDIAL_AXIS_SKELETON
