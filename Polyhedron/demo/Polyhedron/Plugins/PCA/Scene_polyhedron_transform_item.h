#ifndef SCENE_POLYHEDRON_TRANSFORM_ITEM_H
#define SCENE_POLYHEDRON_TRANSFORM_ITEM_H

#include "Scene_polyhedron_item.h"
#include "Scene_polyhedron_transform_item_config.h"
#include <QGLViewer/manipulatedFrame.h>
#include <QGLViewer/qglviewer.h>
#include <QKeyEvent>
struct Scene_polyhedron_transform_item_priv;
// This class represents a polyhedron in the OpenGL scene
class SCENE_POLYHEDRON_TRANSFORM_ITEM_EXPORT Scene_polyhedron_transform_item 
        : public CGAL::Three::Scene_item {
    Q_OBJECT
    
    typedef Scene_polyhedron_item Base;
    
public: 
    Scene_polyhedron_transform_item(const qglviewer::Vec& pos,const Scene_polyhedron_item* poly_item,const CGAL::Three::Scene_interface* scene_interface);
    Scene_item* clone() const{return NULL;}
    QString toolTip() const;
    void drawEdges(CGAL::Three::Viewer_interface*) const;
    void compute_bbox() const;
    ~Scene_polyhedron_transform_item();
    bool manipulatable() const;
    ManipulatedFrame* manipulatedFrame();
    void setManipulatable(bool);
    const Scene_polyhedron_item* getBase() const;
    const qglviewer::Vec& center() const;
    virtual bool supportsRenderingMode(RenderingMode m) const { return m==Wireframe ; }
    virtual void invalidateOpenGLBuffers();
    virtual bool keyPressEvent(QKeyEvent*);
    void setFMatrix(double matrix[16]);

protected:
    friend struct Scene_polyhedron_transform_item_priv;
    Scene_polyhedron_transform_item_priv* d;

Q_SIGNALS:
    void stop();
    void killed();
}; // end class Scene_polyhedron_transform_item

#endif // SCENE_POLYHEDRON_TRANSFORM_ITEM_H
