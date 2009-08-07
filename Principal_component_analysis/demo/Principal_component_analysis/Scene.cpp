#include "Scene.h"

#include <iostream>
#include <fstream>

#include <QString>
#include <QTextStream>
#include <QFileInfo>
#include <QInputDialog>

#include <CGAL/Timer.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Subdivision_method_3.h>

#include <CGAL/centroid.h>
#include <CGAL/linear_least_squares_fitting_3.h>


#include "render_edges.h"

Scene::Scene()
{
    m_pPolyhedron = NULL;

    // view options
    m_view_polyhedron = true;
}

Scene::~Scene()
{
    delete m_pPolyhedron;
}

int Scene::open(QString filename)
{
    QTextStream cerr(stderr);
    cerr << QString("Opening file \"%1\"\n").arg(filename);
    QApplication::setOverrideCursor(QCursor(::Qt::WaitCursor));

    QFileInfo fileinfo(filename);
    std::ifstream in(filename.toUtf8());

    if(!in || !fileinfo.isFile() || ! fileinfo.isReadable())
    {
        std::cerr << "unable to open file" << std::endl;
        QApplication::restoreOverrideCursor();
        return -1;
    }

    if(m_pPolyhedron != NULL)
        delete m_pPolyhedron;

    // allocate new polyhedron
    m_pPolyhedron = new Polyhedron;
    in >> *m_pPolyhedron;
    if(!in)
    {
        std::cerr << "invalid OFF file" << std::endl;
        QApplication::restoreOverrideCursor();

        delete m_pPolyhedron;
        m_pPolyhedron = NULL;

        return -1;
    }

    QApplication::restoreOverrideCursor();
    return 0;
}

void Scene::update_bbox()
{
    std::cout << "Compute bbox...";
    m_bbox = Bbox();

    if(m_pPolyhedron == NULL)
    {
        std::cout << "failed (no polyhedron)." << std::endl;
        return;
    }

    if(m_pPolyhedron->empty())
    {
        std::cout << "failed (empty polyhedron)." << std::endl;
        return;
    }

    Polyhedron::Point_iterator it = m_pPolyhedron->points_begin();
    m_bbox = (*it).bbox();
    for(; it != m_pPolyhedron->points_end();it++)
        m_bbox = m_bbox + (*it).bbox();
    std::cout << "done (" << m_pPolyhedron->size_of_facets()
        << " facets)" << std::endl;
}

void Scene::draw()
{
    if(m_view_polyhedron)
        render_polyhedron();

    render_line();
    render_plane();
    render_centroid();
}

void Scene::render_plane()
{
    ::glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
    ::glLineWidth(3.0f);
    ::glColor3ub(255,0,0);
    ::glBegin(GL_QUADS);
    Point o = m_plane.projection(m_centroid);
    Point a = o + normalize(m_plane.base1()) + normalize(m_plane.base2());
    Point b = o + normalize(m_plane.base1()) - normalize(m_plane.base2());
    Point c = o - normalize(m_plane.base1()) - normalize(m_plane.base2());
    Point d = o - normalize(m_plane.base1()) + normalize(m_plane.base2());
    ::glVertex3d(a.x(),a.y(),a.z());
    ::glVertex3d(b.x(),b.y(),b.z());
    ::glVertex3d(c.x(),c.y(),c.z());
    ::glVertex3d(d.x(),d.y(),d.z());
    ::glEnd();
}

void Scene::render_line()
{
    ::glLineWidth(3.0f);
    ::glColor3ub(0,0,255);
    ::glBegin(GL_LINES);
    Point o = m_line.projection(m_centroid);
    Point a = o + normalize(m_line.to_vector());
    Point b = o - normalize(m_line.to_vector());
    ::glVertex3d(a.x(),a.y(),a.z());
    ::glVertex3d(b.x(),b.y(),b.z());
    ::glEnd();
}

void Scene::render_centroid()
{
    ::glPointSize(10.0f);
    ::glColor3ub(0,128,0);
    ::glBegin(GL_POINTS);
    ::glVertex3d(m_centroid.x(),m_centroid.y(),m_centroid.z());
    ::glEnd();
}


Vector Scene::normalize(const Vector& v)
{
    return v / std::sqrt(v*v);
}

void Scene::render_polyhedron()
{
    // draw black edges
    if(m_pPolyhedron != NULL)
    {
        ::glDisable(GL_LIGHTING);
        ::glColor3ub(0,0,0);
        ::glLineWidth(1.0f);
        gl_render_edges(*m_pPolyhedron);
    }
}

void Scene::refine_loop()
{
    if(m_pPolyhedron == NULL)
    {
        std::cout << "Load polyhedron first." << std::endl;
        return;
    }
    std::cout << "Loop subdivision...";
    CGAL::Subdivision_method_3::Loop_subdivision(*m_pPolyhedron, 1);
    std::cout << "done (" << m_pPolyhedron->size_of_facets() << " facets)" << std::endl;
}

void Scene::toggle_view_poyhedron()
{
    m_view_polyhedron = !m_view_polyhedron;
}

void Scene::fit_triangles()
{
    std::cout << "Fit triangles...";

    std::list<Triangle> triangles;
    Polyhedron::Facet_iterator it;
    for(it = m_pPolyhedron->facets_begin();
        it != m_pPolyhedron->facets_end();
        it++)
    {
        Polyhedron::Halfedge_handle he = it->halfedge();
        const Point& a = he->vertex()->point();
        const Point& b = he->next()->vertex()->point();
        const Point& c = he->next()->next()->vertex()->point();
        Triangle triangle(a,b,c);
        triangles.push_back(triangle);
    }

    m_centroid = CGAL::centroid(triangles.begin(),triangles.end());
    CGAL::linear_least_squares_fitting_3(triangles.begin(),
        triangles.end(), m_line, CGAL::Dimension_tag<2>()); 
    CGAL::linear_least_squares_fitting_3(triangles.begin(),
        triangles.end(), m_plane, CGAL::Dimension_tag<2>()); 

    std::cout << "done" << std::endl;
}

void Scene::fit_edges()
{
    std::cout << "Fit edges...";

    std::list<Segment> segments;
    Polyhedron::Edge_iterator he;
    for(he = m_pPolyhedron->edges_begin();
        he != m_pPolyhedron->edges_end();
        he++)
    {
        const Point& a = he->vertex()->point();
        const Point& b = he->opposite()->vertex()->point();
        Segment segment(a,b);
        segments.push_back(segment);
    }
    
    m_centroid = CGAL::centroid(segments.begin(),segments.end());
    CGAL::linear_least_squares_fitting_3(segments.begin(),
        segments.end(), m_line, CGAL::Dimension_tag<1>()); 
    CGAL::linear_least_squares_fitting_3(segments.begin(),
        segments.end(), m_plane, CGAL::Dimension_tag<1>()); 

    std::cout << "done" << std::endl;
}

void Scene::fit_vertices()
{
    std::cout << "Fit vertices...";

    std::list<Point> points;
    Polyhedron::Vertex_iterator v;
    for(v = m_pPolyhedron->vertices_begin();
        v != m_pPolyhedron->vertices_end();
        v++)
    {
        const Point& p = v->point();
        points.push_back(p);
    }
    
    m_centroid = CGAL::centroid(points.begin(),points.end());
    CGAL::linear_least_squares_fitting_3(points.begin(),
        points.end(), m_line, CGAL::Dimension_tag<0>()); 
    CGAL::linear_least_squares_fitting_3(points.begin(),
        points.end(), m_plane, CGAL::Dimension_tag<0>()); 

    std::cout << "done" << std::endl;
}





