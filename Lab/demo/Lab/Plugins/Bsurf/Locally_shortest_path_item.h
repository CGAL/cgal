#ifndef LOCALLY_SHORTEST_PATH_ITEM_H
#define LOCALLY_SHORTEST_PATH_ITEM_H

#include <CGAL/Three/Scene_item_rendering_helper.h>
#include <CGAL/Simple_cartesian.h>
#include "Scene_surface_mesh_item.h"
#include "Scene_polylines_item.h"
#include "create_sphere.h"
#include "Locally_shortest_path_item_config.h"

struct Locally_shortest_path_item_priv;
class LOCALLY_SHORTEST_PATH_ITEM_EXPORT Locally_shortest_path_item:
    public CGAL::Three::Scene_item_rendering_helper
{
    Q_OBJECT
  public:
    typedef CGAL::Simple_cartesian<double>  Kernel;
    struct vertex;
    struct edge;
    struct face;

    Locally_shortest_path_item(const CGAL::Three::Scene_interface* scene_interface,
                               const Scene_surface_mesh_item* sm_item,
                               Scene_polylines_item* polyline_item,
                               std::size_t nb_pts);
    ~Locally_shortest_path_item();
    bool isFinite() const { return true; }
    bool isEmpty() const { return false; }
    void compute_bbox() const;

    bool manipulatable() const { return true; }
    ManipulatedFrame* manipulatedFrame();
    Locally_shortest_path_item* clone() const {
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
    void drawPath() const;
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
    friend struct Locally_shortest_path_item_priv;
    Locally_shortest_path_item_priv* d;
};
#endif // SCENE_EDIT_BOX_ITEM_H
