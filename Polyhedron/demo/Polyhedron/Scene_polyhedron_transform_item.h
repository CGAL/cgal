#ifndef SCENE_POLYHEDRON_TRANSFORM_ITEM_H
#define SCENE_POLYHEDRON_TRANSFORM_ITEM_H

#include "Scene_polyhedron_item.h"
#include "Scene_polyhedron_transform_item_config.h"
#include <QGLViewer/manipulatedFrame.h>
#include <QGLViewer/qglviewer.h>
#include <QKeyEvent>

// This class represents a polyhedron in the OpenGL scene
class SCENE_POLYHEDRON_TRANSFORM_ITEM_EXPORT Scene_polyhedron_transform_item 
        : public Scene_item {
    Q_OBJECT
    
    typedef Scene_polyhedron_item Base;
    
public: 
    Scene_polyhedron_transform_item(const qglviewer::Vec& pos,const Scene_polyhedron_item* poly_item,const Scene_interface* scene_interface);
    Scene_item* clone() const{return NULL;}
    QString toolTip() const;
    void draw_edges(Viewer_interface*) const;
    Bbox bbox() const;
    ~Scene_polyhedron_transform_item() {delete frame; Q_EMIT killed();}
    bool manipulatable() const { return manipulable; }
    ManipulatedFrame* manipulatedFrame() { return frame; }
    void setManipulatable(bool b = true) { manipulable = b;}
    const Scene_polyhedron_item* getBase() const{ return poly_item;  };
    const qglviewer::Vec& center() const { return center_; }
    virtual bool supportsRenderingMode(RenderingMode m) const { return m==Wireframe ; }
    virtual void invalidate_buffers();
    virtual bool keyPressEvent(QKeyEvent*);

private:
    const Scene_polyhedron_item* poly_item;
    bool manipulable;
    qglviewer::ManipulatedFrame* frame;
    const Polyhedron* poly;
    qglviewer::Vec center_;
    mutable QOpenGLShaderProgram *program;
    mutable std::vector<float> positions_lines;
    mutable std::size_t nb_lines;
    using Scene_item::initialize_buffers;
    void initialize_buffers(Viewer_interface *viewer) const;
    void compute_elements();

Q_SIGNALS:
    void stop();
    void killed();
}; // end class Scene_polyhedron_transform_item

#endif // SCENE_POLYHEDRON_TRANSFORM_ITEM_H
