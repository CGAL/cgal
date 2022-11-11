#ifndef GLWIDGET_H
#define GLWIDGET_H

// Qt
#include <QOpenGLWidget>
#include <QOpenGLFunctions_2_1>
#include <QPaintEvent>

// local
#include "scene.h"

class GlViewer : public QOpenGLWidget, public QOpenGLFunctions_2_1
{
    Q_OBJECT

private:
    Scene* m_scene;

    // toggles
    bool m_view_points;
    bool m_view_tolerance;
    bool m_view_vertices;
    bool m_view_edges;
    bool m_view_ghost_edges;
    bool m_view_edge_cost;
    bool m_view_edge_priority;
    bool m_view_bins;
    bool m_view_foot_points;
    bool m_view_relocation;
    bool m_view_edge_relevance;

    // interactive modes
    bool m_insert_points;

    // rendering options
    double m_line_thickness;
    double m_point_size;
    double m_vertex_size;

    // camera
    double m_scale;
    double m_center_x, m_center_y;

    // mouse
    QPoint m_mouse_click;
    QPoint m_mouse_move;
    Point m_mouse_pick;

public:
    GlViewer(QWidget *parent);

    void set_scene(Scene* pScene) { m_scene = pScene; }

    void set_camera(const double x, const double y, const double s)
    {
        m_center_x = x;
        m_center_y = y;
        m_scale = s;
    }

    // options
    double& line_thickness() { return m_line_thickness; }
    const double& line_thickness() const { return m_line_thickness; }

    double& point_size() { return m_point_size; }
    const double& point_size() const { return m_point_size; }

    double& vertex_size() { return m_vertex_size; }
    const double& vertex_size() const { return m_vertex_size; }

    // toggles
    void toggle_view_points() { m_view_points = !m_view_points; }

    void toggle_view_tolerance() { m_view_tolerance = !m_view_tolerance; }

    void toggle_view_vertices() { m_view_vertices = !m_view_vertices; }

    void toggle_view_edges() { m_view_edges = !m_view_edges; }

    void toggle_view_ghost_edges() { m_view_ghost_edges = !m_view_ghost_edges; }

    void toggle_view_edge_cost() { m_view_edge_cost = !m_view_edge_cost; }

    void toggle_view_edge_priority() {
      m_view_edge_priority = !m_view_edge_priority;
    }

    void toggle_view_bins () { m_view_bins = !m_view_bins; }

    void toggle_view_foot_points() { m_view_foot_points = !m_view_foot_points; }

    void toggle_view_relocation() { m_view_relocation = !m_view_relocation; }

    void toggle_view_edge_relevance() { m_view_edge_relevance = !m_view_edge_relevance; }

    void toggle_insert_points() { m_insert_points = !m_insert_points; }

protected:
    // GL
    void paintGL();
    void initializeGL();
    void resizeGL(int width, int height);

    // mouse
    void wheelEvent(QWheelEvent *event);
    void mouseMoveEvent(QMouseEvent *event);
    void mousePressEvent(QMouseEvent *event);
    void mouseReleaseEvent(QMouseEvent *event);

    void sample_mouse_path(const QPoint& point);
    void move_camera(const QPoint& p0, const QPoint& p1);
    void convert_to_world_space(const QPoint& point, double &x, double &y);
};

#endif
