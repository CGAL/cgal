#include "Scene.h"

#include <iostream>
#include <fstream>

#include <QString>
#include <QTextStream>
#include <QFileInfo>
#include <QInputDialog>

#include "Refiner.h"
#include "render_edges.h"

#include <CGAL/Timer.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Subdivision_method_3.h>

Scene::Scene()
{
    m_pPolyhedron = NULL;

    // view options
    m_view_points = true;
    m_view_segments = true;
    m_view_polyhedron = true;
    m_view_distance_function = true;

    // distance function
    m_red_ramp.build_red();
    m_blue_ramp.build_blue();
    m_max_distance_function = (FT)0.0;
    m_signed_distance_function = false;
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
        draw_polyhedron();

    if(m_view_points)
        draw_points();

    if(m_view_segments)
        draw_segments();

    if(m_view_distance_function)
    {
        if(m_signed_distance_function)
            draw_signed_distance_function();
        else
            draw_unsigned_distance_function();
    }
}

void Scene::draw_polyhedron()
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

void Scene::draw_segments()
{
    if(m_segments.size() != 0)
    {
        ::glDisable(GL_LIGHTING);
        ::glColor3ub(0,100,0);
        ::glLineWidth(2.0f);
        ::glBegin(GL_LINES);
        std::list<Segment>::iterator it;
        for(it = m_segments.begin(); it != m_segments.end(); it++)
        {
            const Segment& s = *it;
            const Point& p = s.source();
            const Point& q = s.target();
            ::glVertex3d(p.x(),p.y(),p.z());
            ::glVertex3d(q.x(),q.y(),q.z());
        }
        ::glEnd();
    }
}

void Scene::draw_points()
{
    // draw red points
    if(m_points.size() != 0)
    {
        ::glDisable(GL_LIGHTING);
        ::glColor3ub(180,0,0);
        ::glPointSize(2.0f);
        ::glBegin(GL_POINTS);
        std::list<Point>::iterator it;
        for(it = m_points.begin(); it != m_points.end(); it++)
        {
            const Point& p = *it;
            ::glVertex3d(p.x(),p.y(),p.z());
        }
        ::glEnd();
    }
}

void Scene::draw_unsigned_distance_function()
{
    if(m_max_distance_function == (FT)0.0)
        return;

    ::glDisable(GL_LIGHTING);
    ::glShadeModel(GL_SMOOTH);
    ::glBegin(GL_QUADS);
    int i,j;
    const int nb_quads = 99;
    for(i=0;i<nb_quads;i++)
    {
        for(j=0;j<nb_quads;j++)
        {
            Point_distance& pd00 = m_distance_function[i][j];
            Point_distance& pd01 = m_distance_function[i][j+1];
            Point_distance& pd11 = m_distance_function[i+1][j+1];
            Point_distance& pd10 = m_distance_function[i+1][j];
            Point& p00 = pd00.first;
            Point& p01 = pd01.first;
            Point& p11 = pd11.first;
            Point& p10 = pd10.first;
            FT& d00 = pd00.second;
            FT& d01 = pd01.second;
            FT& d11 = pd11.second;
            FT& d10 = pd10.second;
            unsigned int i00 = 255-(unsigned int)(255.0 * d00 / m_max_distance_function);
            unsigned int i01 = 255-(unsigned int)(255.0 * d01 / m_max_distance_function);
            unsigned int i11 = 255-(unsigned int)(255.0 * d11 / m_max_distance_function);
            unsigned int i10 = 255-(unsigned int)(255.0 * d10 / m_max_distance_function);
            ::glColor3ub(m_thermal_ramp.r(i00),m_thermal_ramp.g(i00),m_thermal_ramp.b(i00));
            ::glVertex3d(p00.x(),p00.y(),p00.z());
            ::glColor3ub(m_thermal_ramp.r(i01),m_thermal_ramp.g(i01),m_thermal_ramp.b(i01));
            ::glVertex3d(p01.x(),p01.y(),p01.z());
            ::glColor3ub(m_thermal_ramp.r(i11),m_thermal_ramp.g(i11),m_thermal_ramp.b(i11));
            ::glVertex3d(p11.x(),p11.y(),p11.z());
            ::glColor3ub(m_thermal_ramp.r(i10),m_thermal_ramp.g(i10),m_thermal_ramp.b(i10));
            ::glVertex3d(p10.x(),p10.y(),p10.z());
        }
    }
    ::glEnd();
}

void Scene::draw_signed_distance_function()
{
    if(m_max_distance_function == (FT)0.0)
        return;

    ::glDisable(GL_LIGHTING);
    ::glShadeModel(GL_SMOOTH);
    ::glBegin(GL_QUADS);
    int i,j;
    const int nb_quads = 99;
    for(i=0;i<nb_quads;i++)
    {
        for(j=0;j<nb_quads;j++)
        {
            Point_distance& pd00 = m_distance_function[i][j];
            Point_distance& pd01 = m_distance_function[i][j+1];
            Point_distance& pd11 = m_distance_function[i+1][j+1];
            Point_distance& pd10 = m_distance_function[i+1][j];
            Point& p00 = pd00.first;
            Point& p01 = pd01.first;
            Point& p11 = pd11.first;
            Point& p10 = pd10.first;
            FT& d00 = pd00.second;
            FT& d01 = pd01.second;
            FT& d11 = pd11.second;
            FT& d10 = pd10.second;

            // determines grey level
            unsigned int i00 = 255-(unsigned)(255.0 * (double)std::fabs(d00) / m_max_distance_function);
            unsigned int i01 = 255-(unsigned)(255.0 * (double)std::fabs(d01) / m_max_distance_function);
            unsigned int i11 = 255-(unsigned)(255.0 * (double)std::fabs(d11) / m_max_distance_function);
            unsigned int i10 = 255-(unsigned)(255.0 * (double)std::fabs(d10) / m_max_distance_function);

            // assembles one quad
            if(d00 > 0.0)
                ::glColor3ub(m_red_ramp.r(i00),m_red_ramp.g(i00),m_red_ramp.b(i00));
            else
                ::glColor3ub(m_blue_ramp.r(i00),m_blue_ramp.g(i00),m_blue_ramp.b(i00));
            ::glVertex3d(p00.x(),p00.y(),p00.z());

            if(d01 > 0.0)
                ::glColor3ub(m_red_ramp.r(i01),m_red_ramp.g(i01),m_red_ramp.b(i01));
            else
                ::glColor3ub(m_blue_ramp.r(i01),m_blue_ramp.g(i01),m_blue_ramp.b(i01));
            ::glVertex3d(p01.x(),p01.y(),p01.z());

            if(d11 > 0)
                ::glColor3ub(m_red_ramp.r(i11),m_red_ramp.g(i11),m_red_ramp.b(i11));
            else
                ::glColor3ub(m_blue_ramp.r(i11),m_blue_ramp.g(i11),m_blue_ramp.b(i11));
            ::glVertex3d(p11.x(),p11.y(),p11.z());

            if(d10 > 0)
                ::glColor3ub(m_red_ramp.r(i10),m_red_ramp.g(i10),m_red_ramp.b(i10));
            else
                ::glColor3ub(m_blue_ramp.r(i10),m_blue_ramp.g(i10),m_blue_ramp.b(i10));
            ::glVertex3d(p10.x(),p10.y(),p10.z());
        }
    }
    ::glEnd();
}

FT Scene::random_in(const double a,
                    const double b)
{
    double r = rand() / (double)RAND_MAX;
    return (FT)(a + (b - a) * r);
}

Point Scene::random_point(const CGAL::Bbox_3& bbox)
{
    FT x = random_in(bbox.xmin(),bbox.xmax());
    FT y = random_in(bbox.ymin(),bbox.ymax());
    FT z = random_in(bbox.zmin(),bbox.zmax());
    return Point(x,y,z);
}

Vector Scene::random_vector()
{
    FT x = random_in(0.0,1.0);
    FT y = random_in(0.0,1.0);
    FT z = random_in(0.0,1.0);
    return Vector(x,y,z);
}

Ray Scene::random_ray(const CGAL::Bbox_3& bbox)
{
    Point p = random_point(bbox);
    Point q = random_point(bbox);
    return Ray(p,q);
}

Segment Scene::random_segment(const CGAL::Bbox_3& bbox)
{
    Point p = random_point(bbox);
    Point q = random_point(bbox);
    return Segment(p,q);
}

Line Scene::random_line(const CGAL::Bbox_3& bbox)
{
    Point p = random_point(bbox);
    Point q = random_point(bbox);
    return Line(p,q);
}

Plane Scene::random_plane(const CGAL::Bbox_3& bbox)
{
    Point p = random_point(bbox);
    Vector vec = random_vector();
    return Plane(p,vec);
}

void Scene::generate_points_in(const unsigned int nb_points,
                               const double min,
                               const double max)
{
    if(m_pPolyhedron == NULL)
    {
        std::cout << "Load polyhedron first." << std::endl;
        return;
    }

    typedef CGAL::AABB_polyhedron_triangle_primitive<Kernel,Polyhedron> Primitive;
    typedef CGAL::AABB_traits<Kernel, Primitive> Traits;
    typedef CGAL::AABB_tree<Traits> Tree;

    std::cout << "Construct AABB tree...";
    Tree tree(m_pPolyhedron->facets_begin(),m_pPolyhedron->facets_end());
    std::cout << "done." << std::endl;

    CGAL::Timer timer;
    timer.start();
    std::cout << "Generate " << nb_points << " points in interval ["
        << min << ";" << max << "]";

    unsigned int nb_trials = 0;
    Vector vec = random_vector();
    while(m_points.size() < nb_points)
    {
        Point p = random_point(tree.bbox());

        // measure distance
        FT signed_distance = std::sqrt(tree.squared_distance(p));

        // measure sign
        Ray ray(p,vec);
        int nb_intersections = (int)tree.number_of_intersected_primitives(ray);
        if(nb_intersections % 2 != 0)
            signed_distance *= -1.0;

        if(signed_distance >= min &&
            signed_distance <= max)
        {
            m_points.push_back(p);
            if(m_points.size()%(nb_points/10) == 0)
                std::cout << "."; // ASCII progress bar
        }
        nb_trials++;
    }
    double speed = (double)nb_trials / timer.time();
    std::cout << "done (" << nb_trials << " trials, "
        << timer.time() << " s, "
        << speed << " queries/s)" << std::endl;
}


void Scene::generate_inside_points(const unsigned int nb_points)
{
    if(m_pPolyhedron == NULL)
    {
        std::cout << "Load polyhedron first." << std::endl;
        return;
    }

    typedef CGAL::AABB_polyhedron_triangle_primitive<Kernel,Polyhedron> Primitive;
    typedef CGAL::AABB_traits<Kernel, Primitive> Traits;
    typedef CGAL::AABB_tree<Traits> Tree;

    std::cout << "Construct AABB tree...";
    Tree tree(m_pPolyhedron->facets_begin(),m_pPolyhedron->facets_end());
    std::cout << "done." << std::endl;

    CGAL::Timer timer;
    timer.start();
    std::cout << "Generate " << nb_points << " inside points";

    unsigned int nb_trials = 0;
    Vector vec = random_vector();
    while(m_points.size() < nb_points)
    {
        Point p = random_point(tree.bbox());
        Ray ray(p,vec);
        int nb_intersections = (int)tree.number_of_intersected_primitives(ray);
        if(nb_intersections % 2 != 0)
        {
            m_points.push_back(p);
            if(m_points.size()%(nb_points/10) == 0)
                std::cout << "."; // ASCII progress bar
        }
        nb_trials++;
    }
    double speed = (double)nb_trials / timer.time();
    std::cout << "done (" << nb_trials << " trials, "
        << timer.time() << " s, "
        << speed << " queries/s)" << std::endl;
}

void Scene::generate_boundary_segments(const unsigned int nb_slices)
{
    if(m_pPolyhedron == NULL)
    {
        std::cout << "Load polyhedron first." << std::endl;
        return;
    }

    typedef CGAL::AABB_polyhedron_triangle_primitive<Kernel,Polyhedron> Primitive;
    typedef CGAL::AABB_traits<Kernel, Primitive> Traits;
    typedef CGAL::AABB_tree<Traits> Tree;
    typedef Tree::Object_and_primitive_id Object_and_primitive_id;

    std::cout << "Construct AABB tree...";
    Tree tree(m_pPolyhedron->facets_begin(),m_pPolyhedron->facets_end());
    std::cout << "done." << std::endl;

    CGAL::Timer timer;
    timer.start();
    std::cout << "Generate boundary segments from " << nb_slices << " slices: ";

    Vector normal((FT)0.0,(FT)0.0,(FT)1.0);
    unsigned int i;

    const double dz = m_bbox.zmax() - m_bbox.zmin();
    for(i=0;i<nb_slices;i++)
    {
        FT z = m_bbox.zmin() + (FT)i / (FT)nb_slices * dz;
        Point p((FT)0.0, (FT)0.0, z);
        Plane plane(p,normal);

        std::list<Object_and_primitive_id> intersections;
        tree.all_intersections(plane,std::back_inserter(intersections));

        std::list<Object_and_primitive_id>::iterator it;
        for(it = intersections.begin();
            it != intersections.end();
            it++)
        {
            Object_and_primitive_id op = *it;
            CGAL::Object object = op.first;
            Segment segment;
            if(CGAL::assign(segment,object))
                m_segments.push_back(segment);
        }
    }
    std::cout << m_segments.size() << " segments, " << timer.time() << " s." << std::endl;
}

void Scene::generate_boundary_points(const unsigned int nb_points)
{
    if(m_pPolyhedron == NULL)
    {
        std::cout << "Load polyhedron first." << std::endl;
        return;
    }

    typedef CGAL::AABB_polyhedron_triangle_primitive<Kernel,Polyhedron> Primitive;
    typedef CGAL::AABB_traits<Kernel, Primitive> Traits;
    typedef CGAL::AABB_tree<Traits> Tree;
    typedef Tree::Object_and_primitive_id Object_and_primitive_id;

    std::cout << "Construct AABB tree...";
    Tree tree(m_pPolyhedron->facets_begin(),m_pPolyhedron->facets_end());
    std::cout << "done." << std::endl;

    CGAL::Timer timer;
    timer.start();
    std::cout << "Generate boundary points: ";

    unsigned int nb = 0;
    unsigned int nb_lines = 0;
    while(nb < nb_points)
    {
        Line line = random_line(tree.bbox());

        std::list<Object_and_primitive_id> intersections;
        tree.all_intersections(line,std::back_inserter(intersections));
        nb_lines++;

        std::list<Object_and_primitive_id>::iterator it;
        for(it = intersections.begin();
            it != intersections.end();
            it++)
        {
            Object_and_primitive_id op = *it;
            CGAL::Object object = op.first;
            Point point;
            if(CGAL::assign(point,object))
            {
                m_points.push_back(point);
                nb++;
            }
        }
    }
    std::cout << nb_lines << " line queries, " << timer.time() << " s." << std::endl;
}

void Scene::generate_edge_points(const unsigned int nb_points)
{
    if(m_pPolyhedron == NULL)
    {
        std::cout << "Load polyhedron first." << std::endl;
        return;
    }

    typedef CGAL::AABB_polyhedron_segment_primitive<Kernel,Polyhedron> Primitive;
    typedef CGAL::AABB_traits<Kernel, Primitive> Traits;
    typedef CGAL::AABB_tree<Traits> Tree;
    typedef Tree::Object_and_primitive_id Object_and_primitive_id;

    std::cout << "Construct AABB tree...";
    Tree tree(m_pPolyhedron->edges_begin(),m_pPolyhedron->edges_end());
    std::cout << "done." << std::endl;

    CGAL::Timer timer;
    timer.start();
    std::cout << "Generate edge points: ";

    unsigned int nb = 0;
    unsigned int nb_planes = 0;
    while(nb < nb_points)
    {
        Plane plane = random_plane(tree.bbox());

        std::list<Object_and_primitive_id> intersections;
        tree.all_intersections(plane,std::back_inserter(intersections));
        nb_planes++;

        std::list<Object_and_primitive_id>::iterator it;
        for(it = intersections.begin();
            it != intersections.end();
            it++)
        {
            Object_and_primitive_id op = *it;
            CGAL::Object object = op.first;
            Point point;
            if(CGAL::assign(point,object))
            {
                m_points.push_back(point);
                nb++;
            }
        }
    }
    std::cout << nb_planes << " plane queries, " << timer.time() << " s." << std::endl;
}

void Scene::unsigned_distance_function()
{
    if(m_pPolyhedron == NULL)
    {
        std::cout << "Load polyhedron first." << std::endl;
        return;
    }

    CGAL::Timer timer;
    timer.start();
    std::cout << "Construct AABB tree...";
    Facet_tree tree(m_pPolyhedron->facets_begin(),m_pPolyhedron->facets_end());
    tree.accelerate_distance_queries();
    std::cout << "done (" << timer.time() << " s)" << std::endl;

    m_max_distance_function = (FT)0.0;
    int i,j;
    const double dx = m_bbox.xmax() - m_bbox.xmin();
    const double dy = m_bbox.ymax() - m_bbox.ymin();
    const double z = 0.5 * (m_bbox.zmax() + m_bbox.zmin());
    for(i=0;i<100;i++)
    {
        FT x = m_bbox.xmin() + (FT)((double)i/100.0 * dx);
        for(j=0;j<100;j++)
        {
            FT y = m_bbox.ymin() + (FT)((double)j/100.0 * dy);
            Point query(x,y,z);
            FT sq_distance = tree.squared_distance(query);
            FT distance = std::sqrt(sq_distance);
            m_distance_function[i][j] = Point_distance(query,distance);
            m_max_distance_function = distance > m_max_distance_function ?
distance : m_max_distance_function;
        }
    }
    m_signed_distance_function = false;
}

void Scene::unsigned_distance_function_to_edges()
{
    if(m_pPolyhedron == NULL)
    {
        std::cout << "Load polyhedron first." << std::endl;
        return;
    }

    typedef CGAL::AABB_polyhedron_segment_primitive<Kernel,Polyhedron> Primitive;
    typedef CGAL::AABB_traits<Kernel, Primitive> Traits;
    typedef CGAL::AABB_tree<Traits> Edge_tree;

    CGAL::Timer timer;
    timer.start();
    std::cout << "Construct AABB tree from edges...";
    Edge_tree tree(m_pPolyhedron->edges_begin(),m_pPolyhedron->edges_end());
    tree.accelerate_distance_queries();
    std::cout << "done (" << timer.time() << " s)" << std::endl;

    m_max_distance_function = (FT)0.0;
    const double dx = m_bbox.xmax() - m_bbox.xmin();
    const double dy = m_bbox.ymax() - m_bbox.ymin();
    const double z = 0.5 * (m_bbox.zmax() + m_bbox.zmin());
    int i,j;
    for(i=0;i<100;i++)
    {
        FT x = m_bbox.xmin() + (FT)((double)i/100.0 * dx);
        for(j=0;j<100;j++)
        {
            FT y = m_bbox.ymin() + (FT)((double)j/100.0 * dy);
            Point query(x,y,z);
            FT sq_distance = tree.squared_distance(query);
            FT distance = std::sqrt(sq_distance);
            m_distance_function[i][j] = Point_distance(query,distance);
            m_max_distance_function = distance > m_max_distance_function ?
distance : m_max_distance_function;
        }
    }
    m_signed_distance_function = false;
}

void Scene::signed_distance_function()
{
    if(m_pPolyhedron == NULL)
    {
        std::cout << "Load polyhedron first." << std::endl;
        return;
    }

    CGAL::Timer timer;
    timer.start();
    std::cout << "Construct AABB tree...";
    Facet_tree tree(m_pPolyhedron->facets_begin(),m_pPolyhedron->facets_end());
    tree.accelerate_distance_queries();
    std::cout << "done (" << timer.time() << " s)" << std::endl;

    m_max_distance_function = (FT)0.0;
    Vector vec = random_vector();

    const double dx = m_bbox.xmax() - m_bbox.xmin();
    const double dy = m_bbox.ymax() - m_bbox.ymin();
    const double z = 0.5 * (m_bbox.zmax() + m_bbox.zmin());
    int i,j;
    for(i=0;i<100;i++)
    {
        FT x = m_bbox.xmin() + (FT)((double)i/100.0 * dx);
        for(j=0;j<100;j++)
        {
            FT y = m_bbox.ymin() + (FT)((double)j/100.0 * dy);
            Point query(x,y,z);
            FT sq_distance = tree.squared_distance(query);
            FT unsigned_distance = std::sqrt(sq_distance);

            // get sign through ray casting (random vector)
            Ray ray(query,vec);
            unsigned int nbi = tree.number_of_intersected_primitives(ray);
            FT sign = nbi%2 == 0 ? (FT)1.0 : (FT)-1.0;
            FT signed_distance = sign * unsigned_distance;

            m_distance_function[i][j] = Point_distance(query,signed_distance);
            m_max_distance_function = unsigned_distance > m_max_distance_function ?
unsigned_distance : m_max_distance_function;
        }
    }
    m_signed_distance_function = true;
}

void Scene::toggle_view_poyhedron()
{
    m_view_polyhedron = !m_view_polyhedron;
}

void Scene::toggle_view_segments()
{
    m_view_segments = !m_view_segments;
}

void Scene::toggle_view_points()
{
    m_view_points = !m_view_points;
}

void Scene::toggle_view_distance_function()
{
    m_view_distance_function = !m_view_distance_function;
}

void Scene::refine_bisection(const FT max_sqlen)
{
    if(m_pPolyhedron == NULL)
    {
        std::cout << "Load polyhedron first." << std::endl;
        return;
    }
    std::cout << "Refine through recursive longest edge bisection...";
    Refiner<Kernel,Polyhedron> refiner(m_pPolyhedron);
    refiner(max_sqlen);
    std::cout << "done (" << m_pPolyhedron->size_of_facets() << " facets)" << std::endl;
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
