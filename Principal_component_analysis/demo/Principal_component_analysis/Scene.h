#ifndef SCENE_H
#define SCENE_H


#include <iostream>
#include <cmath>
#include <QOpenGLShaderProgram>
#include <QOpenGLBuffer>
#include <QOpenGLVertexArrayObject>
#include "types.h"
class Viewer;

class Scene
{
public:
    Scene();
    ~Scene();
public:
    // types
    typedef CGAL::Bbox_3 Bbox;

public:
    void update_bbox();
    Bbox bbox() { return m_bbox; }

private:
    // member data
    Bbox m_bbox;
    Line m_line;
    Plane m_plane;
    Point m_centroid;
    Polyhedron *m_pPolyhedron;
    QOpenGLShaderProgram rendering_program;
    QOpenGLBuffer buffers[4];
    QOpenGLVertexArrayObject vao[4];
    bool is_gl_init;
    void gl_init();

    // view options
    bool m_view_polyhedron;

public:
    // file menu
    int open(QString filename);

    // toggle view options
    void toggle_view_poyhedron();

    // algorithms
    Vector normalize(const Vector& v);

    void refine_loop();
    void fit_edges();
    void fit_vertices();
    void fit_triangles();

    // rendering
    void draw(Viewer *viewer);
    void render_line(Viewer *viewer);
    void render_plane(Viewer* viewer);
    void render_centroid(Viewer* viewer);
    void render_polyhedron(Viewer* viewer);

private:




}; // end class Scene


#endif // SCENE_H
