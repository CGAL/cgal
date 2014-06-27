// Qt
#include <QtGui>

// local
#include "glviewer.h"

GlViewer::GlViewer(QWidget *pParent)
: QGLWidget(QGLFormat(QGL::SampleBuffers), pParent)
{
	m_scene = NULL;
    
	m_view_points = true;
	m_view_vertices      = true;
	m_view_edges         = false;
	m_view_ghost_edges   = false;
	m_view_edge_cost     = false;
	m_view_edge_priority = false;
	m_view_bins        = false;
	m_view_foot_points = false;
	m_view_relocation  = false;
	m_view_tolerance   = false;
	m_view_incolors       = false;
	m_view_edge_relevance = true;
    
	m_insert_points      = false;
	m_activate_simulation = 0;
	m_simulation_stage    = 0;

	m_line_thickness = 2.0;
	m_point_size = 2.0;
	m_vertex_size = 2.0;

	m_scale = 1.0;
	m_center_x = m_center_y = 0.5;
    
	setAutoFillBackground(false);
}

void GlViewer::resizeGL(int width, int height) 
{
	glViewport(0, 0, width, height);
	double aspect_ratio = double(height) / double(width);
    
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-1.0, 1.0, -aspect_ratio, aspect_ratio, -1.0, 1.0);
    
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

void GlViewer::initializeGL() 
{
	glClearColor(1., 1., 1., 0.);
	glDisable(GL_DEPTH_TEST);
	glEnable(GL_SMOOTH);
}

void GlViewer::paintGL() 
{
	glClear(GL_COLOR_BUFFER_BIT);
	if (!m_scene) return;
    
	glPushMatrix();
	glScaled(m_scale, m_scale, m_scale);
	glTranslated(-m_center_x, -m_center_y, 0.0);

	m_scene->render(m_view_points,
                    m_view_vertices, 
                    m_view_edges, 
                    m_view_ghost_edges,
                    m_view_edge_cost,
                    m_view_edge_priority,
                    m_view_bins,
                    m_view_foot_points,
                    m_view_relocation,                    
                    m_view_tolerance,
                    m_view_incolors,
                    m_view_edge_relevance,
                    float(m_point_size),
                    float(m_vertex_size),
                    float(m_line_thickness));
            
    if (m_activate_simulation == 2)
        m_scene->render_simulation(m_mouse_pick, 
                                   m_simulation_stage, 
                                   float(2*m_vertex_size), 
                                   float(m_line_thickness));
    
	glPopMatrix();
}

void GlViewer::wheelEvent(QWheelEvent *event) 
{
	if (!m_scene) return;
	m_scale += 0.05 * (event->delta() / 120);
    if (m_scale <= 0.0) m_scale = 0.0;
	updateGL();
}

void GlViewer::mousePressEvent(QMouseEvent *event) 
{
	if (!m_scene) return;
	m_mouse_click = event->pos();
    
	if (event->button() == Qt::LeftButton) 
    {
        if (m_activate_simulation != 0) m_activate_simulation = 1;
		setCursor(QCursor(Qt::PointingHandCursor));
		sample_mouse_path(m_mouse_click);
	} 
    else 
    {
		setCursor(QCursor(Qt::ClosedHandCursor));
	}
}

void GlViewer::mouseMoveEvent(QMouseEvent *event)
{
	if(!m_scene) return;    
	m_mouse_move = event->pos();
    
	if (event->buttons() == Qt::LeftButton)
	{
        if (m_mouse_move != m_mouse_click)
            sample_mouse_path(m_mouse_move);
	}
    else 
    {        
		move_camera(m_mouse_click, m_mouse_move);
    }
    
    m_mouse_click = m_mouse_move;
	updateGL();
}

void GlViewer::mouseReleaseEvent(QMouseEvent *event) 
{
	if (!m_scene) return;
	m_mouse_move = event->pos();
    
	if (event->button() == Qt::LeftButton) 
    {
        if (m_mouse_move != m_mouse_click)
            sample_mouse_path(m_mouse_move);
	} 
    else 
    {
		move_camera(m_mouse_click, m_mouse_move);
	}    
    
    m_mouse_click = m_mouse_move;
	setCursor(QCursor(Qt::ArrowCursor));
	updateGL();
}

void GlViewer::sample_mouse_path(const QPoint& point)
{
	double x, y;
	convert_to_world_space(point, x, y); 
    
	if (m_insert_points)
		m_scene->add_sample(Point(x, y));
    
    if (m_activate_simulation == 1)
    {
        m_activate_simulation = 2;
        m_mouse_pick = Point(x, y);
    }
}

void GlViewer::move_camera(const QPoint& p0, const QPoint& p1)
{
	m_center_x -= double(p1.x() - p0.x()) / double(width());
	m_center_y += double(p1.y() - p0.y()) / double(height()); 
}

void GlViewer::convert_to_world_space(const QPoint& point, double &x, double &y)
{
	double aspect_ratio = double(height()) / double(width());
    
	x = double(point.x()) / double(width());
	x = (2.0*x - 1.0) / m_scale;
	x += m_center_x;
    
	y = 1.0 - double(point.y()) / double(height());
	y = (2.0*y - 1.0) * aspect_ratio / m_scale;
	y += m_center_y;
}
