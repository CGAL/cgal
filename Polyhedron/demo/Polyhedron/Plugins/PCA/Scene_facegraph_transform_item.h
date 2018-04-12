#ifndef Scene_facegraph_transform_item_H
#define Scene_facegraph_transform_item_H

#include "Scene_facegraph_transform_item_config.h"
#include "Kernel_type.h"

#ifdef CGAL_USE_SURFACE_MESH
#include <CGAL/Surface_mesh/Surface_mesh_fwd.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>
#include <CGAL/Three/Scene_interface.h>
#include <CGAL/Three/Scene_item.h>
#else
#include "Scene_polyhedron_item.h"
#endif

#include <CGAL/Qt/manipulatedFrame.h>
#include <CGAL/Qt/qglviewer.h>

#include <QKeyEvent>



using namespace CGAL::Three;
struct Scene_facegraph_transform_item_priv;
#ifdef CGAL_USE_SURFACE_MESH
typedef CGAL::Surface_mesh<Kernel::Point_3> FaceGraph;
#else
typedef Polyhedron FaceGraph;
#endif
// This class represents a polyhedron in the OpenGL scene
class SCENE_FACEGRAPH_TRANSFORM_ITEM_EXPORT Scene_facegraph_transform_item
        : public Scene_item {
    Q_OBJECT
    
public: 
    Scene_facegraph_transform_item(const CGAL::qglviewer::Vec& pos, FaceGraph *sm,
                                    const QString name);
    Scene_item* clone() const{return NULL;}
    QString toolTip() const;
    void drawEdges(CGAL::Three::Viewer_interface*) const;
    void compute_bbox() const;
    ~Scene_facegraph_transform_item();
    bool manipulatable() const;
    ManipulatedFrame* manipulatedFrame();
    void setManipulatable(bool);
    const CGAL::qglviewer::Vec& center() const;
    virtual bool supportsRenderingMode(RenderingMode m) const { return m==Wireframe ; }
    virtual void invalidateOpenGLBuffers();
    virtual bool keyPressEvent(QKeyEvent*);
    void setFMatrix(double matrix[16]);
    bool isEmpty() const {return false;}
    FaceGraph* getFaceGraph();

protected:
    friend struct Scene_facegraph_transform_item_priv;
    Scene_facegraph_transform_item_priv* d;

Q_SIGNALS:
    void stop();
    void killed();
}; // end class Scene_facegraph_transform_item

#endif // Scene_facegraph_transform_item_H
