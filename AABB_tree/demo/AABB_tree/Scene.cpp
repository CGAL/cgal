#include "Scene.h"

#include <iostream>
#include <fstream>

#include <QtCore/qglobal.h>
#include <QString>
#include <QTextStream>
#include <QFileInfo>
#include <QInputDialog>

#include "Refiner.h"
#include "render_edges.h"

#include <CGAL/Timer.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Subdivision_method_3.h>


// constants
const int slow_distance_grid_size = 100;
const int fast_distance_grid_size = 20;

Scene::Scene()
  : m_frame (new ManipulatedFrame())
  , m_view_plane(false)
  , m_grid_size(slow_distance_grid_size)
  , m_cut_plane(NONE)
{
    m_pPolyhedron = NULL;

    // view options
    m_view_points = true;
    m_view_segments = true;
    m_view_polyhedron = true;

    // distance function
    m_red_ramp.build_red();
    m_blue_ramp.build_blue();
    m_max_distance_function = (FT)0.0;
}

Scene::~Scene()
{
    delete m_pPolyhedron;
    delete m_frame;
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
    
    // clear tree
    clear_internal_data();

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
    if(m_view_plane)
        ::glEnable(GL_DEPTH_TEST);
    else
        ::glDisable(GL_DEPTH_TEST);
  
    if(m_view_polyhedron)
        draw_polyhedron();

    if(m_view_points)
        draw_points();

    if(m_view_segments)
        draw_segments();

    if (m_view_plane)
    {
        switch( m_cut_plane )
        {
          case UNSIGNED_EDGES:
          case UNSIGNED_FACETS:
              draw_distance_function(m_thermal_ramp, m_thermal_ramp);
              break;
          case SIGNED_FACETS:
              draw_distance_function(m_red_ramp, m_blue_ramp);
              break;
          case CUT_SEGMENTS:
              draw_cut_segment_plane();
              break;
          case NONE: // do nothing
              break;
        }
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

void Scene::draw_distance_function(const Color_ramp& ramp_pos,
                                   const Color_ramp& ramp_neg) const
{
    ::glDisable(GL_LIGHTING);
    if ( m_fast_distance ) { ::glShadeModel(GL_FLAT); }
    else { ::glShadeModel(GL_SMOOTH); }
    
    ::glBegin(GL_QUADS);
    int i,j;
    const int nb_quads = m_grid_size-1;
    for(i=0;i<nb_quads;i++)
    {
        for(j=0;j<nb_quads;j++)
        {
            const Point_distance& pd00 = m_distance_function[i][j];
            const Point_distance& pd01 = m_distance_function[i][j+1];
            const Point_distance& pd11 = m_distance_function[i+1][j+1];
            const Point_distance& pd10 = m_distance_function[i+1][j];
            const Point& p00 = pd00.first;
            const Point& p01 = pd01.first;
            const Point& p11 = pd11.first;
            const Point& p10 = pd10.first;
            const FT& d00 = pd00.second;
            const FT& d01 = pd01.second;
            const FT& d11 = pd11.second;
            const FT& d10 = pd10.second;
            
            // determines grey level
            unsigned int i00 = 255-(unsigned)(255.0 * (double)std::fabs(d00) / m_max_distance_function);
            unsigned int i01 = 255-(unsigned)(255.0 * (double)std::fabs(d01) / m_max_distance_function);
            unsigned int i11 = 255-(unsigned)(255.0 * (double)std::fabs(d11) / m_max_distance_function);
            unsigned int i10 = 255-(unsigned)(255.0 * (double)std::fabs(d10) / m_max_distance_function);
            
            // assembles one quad
            if(d00 > 0.0)
                ::glColor3ub(ramp_pos.r(i00),ramp_pos.g(i00),ramp_pos.b(i00));
            else
                ::glColor3ub(ramp_neg.r(i00),ramp_neg.g(i00),ramp_neg.b(i00));
            ::glVertex3d(p00.x(),p00.y(),p00.z());
            
            if(d01 > 0.0)
                ::glColor3ub(ramp_pos.r(i01),ramp_pos.g(i01),ramp_pos.b(i01));
            else
                ::glColor3ub(ramp_neg.r(i01),ramp_neg.g(i01),ramp_neg.b(i01));
            ::glVertex3d(p01.x(),p01.y(),p01.z());
            
            if(d11 > 0)
                ::glColor3ub(ramp_pos.r(i11),ramp_pos.g(i11),ramp_pos.b(i11));
            else
                ::glColor3ub(ramp_neg.r(i11),ramp_neg.g(i11),ramp_neg.b(i11));
            ::glVertex3d(p11.x(),p11.y(),p11.z());
            
            if(d10 > 0)
                ::glColor3ub(ramp_pos.r(i10),ramp_pos.g(i10),ramp_pos.b(i10));
            else
                ::glColor3ub(ramp_neg.r(i10),ramp_neg.g(i10),ramp_neg.b(i10));
            ::glVertex3d(p10.x(),p10.y(),p10.z());
        }
    }
    ::glEnd();
}

void Scene::draw_cut_segment_plane() const
{
    float diag = .6f * float(bbox_diag());

    ::glDisable(GL_LIGHTING);
    ::glLineWidth(1.0f);
    ::glColor3f(.6f, .6f, .6f);

    // draw grid
    ::glPushMatrix();
    ::glMultMatrixd(m_frame->matrix());
    QGLViewer::drawGrid(diag);
    ::glPopMatrix();

    // draw cut segments
    ::glLineWidth(2.0f);
    ::glColor3f(1.f, 0.f, 0.f);
    ::glBegin(GL_LINES);
    for ( std::vector<Segment>::const_iterator it = m_cut_segments.begin(), 
          end = m_cut_segments.end() ; it != end ; ++it )
    {
        const Point& a = it->source();
        const Point& b = it->target();
      
        ::glVertex3d(a.x(), a.y(), a.z());
        ::glVertex3d(b.x(), b.y(), b.z());
    }
    ::glEnd();
  
    // fill grid with transparent blue
    ::glPushMatrix();
    ::glMultMatrixd(m_frame->matrix());
    ::glColor4f(.6f, .85f, 1.f, .65f);

    ::glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA); 
    ::glEnable(GL_BLEND);
    ::glBegin(GL_QUADS);
    ::glVertex3d(-diag, -diag, 0.);
    ::glVertex3d(-diag,  diag, 0.);
    ::glVertex3d( diag,  diag, 0.);
    ::glVertex3d( diag, -diag, 0.);
    ::glEnd();
    ::glDisable(GL_BLEND);
  
    ::glPopMatrix();
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

Plane Scene::frame_plane() const
{
    const qglviewer::Vec& pos = m_frame->position();
    const qglviewer::Vec& n = m_frame->inverseTransformOf(qglviewer::Vec(0.f, 0.f, 1.f));

    return Plane(n[0], n[1],  n[2], - n * pos);
}

Aff_transformation Scene::frame_transformation() const
{
    const ::GLdouble* m = m_frame->matrix();
  
    // OpenGL matrices are row-major matrices
    return Aff_transformation (m[0], m[4], m[8], m[12],
                               m[1], m[5], m[9], m[13],
                               m[2], m[6], m[10], m[14]);
}

FT Scene::bbox_diag() const
{
  double dx = m_bbox.xmax()-m_bbox.xmin();
  double dy = m_bbox.ymax()-m_bbox.ymin();
  double dz = m_bbox.zmax()-m_bbox.zmin();
  
  return FT(std::sqrt(dx*dx + dy*dy + dz*dz));
}

void Scene::build_facet_tree()
{
    if ( NULL == m_pPolyhedron )
    {
        std::cerr << "Build facet tree failed: load polyhedron first." << std::endl;
        return;
    }

    // Don't rebuild tree if it is already built
    if ( !m_facet_tree.empty() ) { return; }
  
    // build tree
    CGAL::Timer timer;
    timer.start();
    std::cout << "Construct Facet AABB tree...";
    m_facet_tree.rebuild(m_pPolyhedron->facets_begin(),m_pPolyhedron->facets_end());
    m_facet_tree.accelerate_distance_queries();
    std::cout << "done (" << timer.time() << " s)" << std::endl;
}

void Scene::build_edge_tree()
{
    if ( NULL == m_pPolyhedron )
    {
        std::cerr << "Build edge tree failed: load polyhedron first." << std::endl;
        return;
    }
  
    // Don't rebuild tree if it is already built
    if ( !m_edge_tree.empty() ) { return; }
    
    // build tree
    CGAL::Timer timer;
    timer.start();
    std::cout << "Construct Edge AABB tree...";
    m_edge_tree.rebuild(m_pPolyhedron->edges_begin(),m_pPolyhedron->edges_end());
    m_edge_tree.accelerate_distance_queries();
    std::cout << "done (" << timer.time() << " s)" << std::endl;
}

void Scene::clear_internal_data()
{
    m_facet_tree.clear();
    m_edge_tree.clear();

    clear_points();
    clear_segments();
    clear_cutting_plane();
}

void Scene::clear_cutting_plane()
{
    m_cut_segments.clear();
    m_cut_plane = NONE;
  
    deactivate_cutting_plane();
}

void Scene::update_grid_size()
{
    m_grid_size = m_fast_distance ? fast_distance_grid_size
                                  : slow_distance_grid_size;
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


template <typename Tree>
void Scene::compute_distance_function(const Tree& tree)
{
    // Get transformation
    Aff_transformation t = frame_transformation();
    
    m_max_distance_function = FT(0);
    FT diag = bbox_diag();
    
    const FT dx = diag;
    const FT dy = diag;
    const FT z (0);
    
    for(int i=0 ; i<m_grid_size ; ++i)
    {
        FT x = -diag/FT(2) + FT(i)/FT(m_grid_size) * dx;
        
        for(int j=0 ; j<m_grid_size ; ++j)
        {
            FT y = -diag/FT(2) + FT(j)/FT(m_grid_size) * dy;
            
            Point query = t( Point(x,y,z) );
            FT dist = CGAL::sqrt( tree.squared_distance(query) );
            
            m_distance_function[i][j] = Point_distance(query,dist);
            m_max_distance_function = (std::max)(dist, m_max_distance_function);
        }
    }
}

template <typename Tree>
void Scene::sign_distance_function(const Tree& tree)
{
    typedef typename Tree::size_type size_type;
    Vector random_vec = random_vector();
    
    for(int i=0 ; i<m_grid_size ; ++i)
    {
        for(int j=0 ; j<m_grid_size ; ++j)
        {
            const Point& p = m_distance_function[i][j].first;
            const FT unsigned_distance = m_distance_function[i][j].second;
            
            // get sign through ray casting (random vector)
            Ray ray(p, random_vec);
            size_type nbi = tree.number_of_intersected_primitives(ray);
            
            FT sign ( (nbi&1) == 0 ? 1 : -1);
            m_distance_function[i][j].second = sign * unsigned_distance;
        }
    }
}


void Scene::unsigned_distance_function()
{
    // Build tree (if build fail, exit)
    build_facet_tree();
    if ( m_facet_tree.empty() ) { return; }
  
    compute_distance_function(m_facet_tree);
    
    m_cut_plane = UNSIGNED_FACETS;
}


void Scene::unsigned_distance_function_to_edges()
{
    // Build tree (if build fail, exit)
    build_edge_tree();
    if ( m_edge_tree.empty() ) { return; }
    
    compute_distance_function(m_edge_tree);
    
    m_cut_plane = UNSIGNED_EDGES;
}


void Scene::signed_distance_function()
{
    // Build tree (if build fail, exit)
    build_facet_tree();
    if ( m_facet_tree.empty() ) { return; }
    
    compute_distance_function(m_facet_tree);
    sign_distance_function(m_facet_tree);

    m_cut_plane = SIGNED_FACETS;
}


void Scene::cut_segment_plane()
{
    // Build tree (if build fail, exit)
    build_facet_tree();
    if ( m_facet_tree.empty() ) { return; }
    
    Plane plane = frame_plane();
    
    // Compute intersections
    typedef std::vector<Facet_tree::Object_and_primitive_id> Intersections;
    Intersections intersections;
    m_facet_tree.all_intersections(plane, std::back_inserter(intersections));
    
    // Fill data structure
    m_cut_segments.clear();
    for ( Intersections::iterator it = intersections.begin(),
         end = intersections.end() ; it != end ; ++it )
    {
        const Segment* inter_seg = CGAL::object_cast<Segment>(&(it->first));
        
        if ( NULL != inter_seg )
        {
            m_cut_segments.push_back(*inter_seg);
        }
    }
    
    m_cut_plane = CUT_SEGMENTS;
}

void Scene::cutting_plane()
{
    switch( m_cut_plane )
    {
      case UNSIGNED_FACETS:
          return unsigned_distance_function();
      case SIGNED_FACETS:
          return signed_distance_function();
      case UNSIGNED_EDGES:
          return unsigned_distance_function_to_edges();
      case CUT_SEGMENTS:
          return cut_segment_plane();
      case NONE: // do nothing 
          return;
    }
    
    // Should not be here
    std::cerr << "Unknown cut_plane type" << std::endl;
    CGAL_assertion(false);
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

void Scene::toggle_view_plane()
{
    m_view_plane = !m_view_plane;
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
  
    clear_internal_data();
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
  
    clear_internal_data();
}


void Scene::activate_cutting_plane()
{
    connect(m_frame, SIGNAL(modified()), this, SLOT(cutting_plane()));
    m_view_plane = true;
}

void Scene::deactivate_cutting_plane()
{
    disconnect(m_frame, SIGNAL(modified()), this, SLOT(cutting_plane()));
    m_view_plane = false;
}
