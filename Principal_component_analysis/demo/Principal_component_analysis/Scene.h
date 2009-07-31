#ifndef SCENE_H
#define SCENE_H

#include <QtOpenGL/qgl.h>
#include <iostream>
#include <cmath>

#include "types.h"


class Scene
{
public:
    Scene();
    ~Scene();
public:
    // types
    typedef CGAL::Bbox_3 Bbox;

public:
    void draw(); 
    void update_bbox();
    Bbox bbox() { return m_bbox; }

private:
    // member data
    Bbox m_bbox;
    Polyhedron *m_pPolyhedron;

private:

public:
    // file menu
    int open(QString filename);

    // toggle view options
    void toggle_view_poyhedron();

    // view options
    bool m_view_polyhedron;

    // refinement
    void refine_loop();

    // drawing
    void draw_polyhedron();
}; // end class Scene


#endif // SCENE_H
