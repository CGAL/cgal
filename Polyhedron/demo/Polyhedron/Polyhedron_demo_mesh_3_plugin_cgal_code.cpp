#include <CGAL/AABB_intersections.h>
#include <CGAL/AABB_tree.h>

#include "Polyhedron_type.h"

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
#include <CGAL/make_mesh_3.h>

#include <Scene_polyhedron_item.h>
#include <Scene_polygon_soup_item.h>
#include <fstream>
#include <sstream>

#include <CGAL/Timer.h>

#include <QMenu>

typedef CGAL::Polyhedral_mesh_domain_with_features_3<Kernel,
Polyhedron> Mesh_domain;

// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;

// 3D complex
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

// Mesh Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
typedef Mesh_criteria::Edge_criteria Edge_criteria;
typedef Mesh_criteria::Facet_criteria Facet_criteria;
typedef Mesh_criteria::Cell_criteria Cell_criteria;

typedef Tr::Point Point_3;

#include "Scene_item.h"
#include <QtCore/qglobal.h>
#include <CGAL/gl.h>
#include <QGLViewer/manipulatedFrame.h>
#include <QGLViewer/qglviewer.h>


class Q_DECL_EXPORT Scene_c3t3_item : public Scene_item
{
    Q_OBJECT
public:
    typedef qglviewer::ManipulatedFrame ManipulatedFrame;

    Scene_c3t3_item(const C3t3& c3t3)
        : Scene_item(7,3), c3t3_(c3t3), frame(new ManipulatedFrame()), last_known_scene(NULL)
    {
        positions_lines.resize(0);
        positions_poly.resize(0);
        color_lines.resize(0);
        color_poly.resize(0);
        color_grid.resize(0);
        normals.resize(0);
        //Generates an integer which will be used as ID for each buffer
    }

    ~Scene_c3t3_item()
    {
        delete frame;
    }

    void invalidate_buffers()
    {
        compute_elements();
        are_buffers_filled = false;
    }

    void contextual_changed()
    {
        if(frame->isManipulated()||frame->isSpinning())
            invalidate_buffers();
    }
    const C3t3& c3t3() const {
        return c3t3_;
    }

    bool manipulatable() const {
        return true;
    }
    ManipulatedFrame* manipulatedFrame() {
        return frame;
    }

    void setPosition(float x, float y, float z) {
        frame->setPosition(x, y, z);
    }

    void setNormal(float x, float y, float z) {
        frame->setOrientation(x, y, z, 0.f);
    }

    Kernel::Plane_3 plane() const {
        const qglviewer::Vec& pos = frame->position();
        const qglviewer::Vec& n =
                frame->inverseTransformOf(qglviewer::Vec(0.f, 0.f, 1.f));
        return Kernel::Plane_3(n[0], n[1],  n[2], - n * pos);
    }

    bool isFinite() const { return true; }
    bool isEmpty() const {
        return c3t3().triangulation().number_of_vertices() == 0;
    }

    Bbox bbox() const {
        if(isEmpty())
            return Bbox();
        else {
            CGAL::Bbox_3 result = c3t3().triangulation().finite_vertices_begin()->point().bbox();
            for(Tr::Finite_vertices_iterator
                vit = ++c3t3().triangulation().finite_vertices_begin(),
                end = c3t3().triangulation().finite_vertices_end();
                vit != end; ++vit)
            {
                result = result + vit->point().bbox();
            }
            return Bbox(result.xmin(), result.ymin(), result.zmin(),
                        result.xmax(), result.ymax(), result.zmax());
        }
    }

    Scene_c3t3_item* clone() const {
        return 0;
    }

    QString toolTip() const {
        int number_of_tets = 0;
        for(Tr::Finite_cells_iterator
            cit = c3t3().triangulation().finite_cells_begin(),
            end = c3t3().triangulation().finite_cells_end();
            cit != end; ++cit)
        {
            if( c3t3().is_in_complex(cit) )
                ++number_of_tets;
        }
        return tr("<p><b>3D complex in a 3D triangulation</b></p>"
                  "<p>Number of vertices: %1<br />"
                  "Number of surface facets: %2<br />"
                  "Number of volume tetrahedra: %3</p>")
                .arg(c3t3().triangulation().number_of_vertices())
                .arg(c3t3().number_of_facets())
                .arg(number_of_tets);
    }

    // Indicate if rendering mode is supported
    bool supportsRenderingMode(RenderingMode m) const {
        return (m != Gouraud && m!=PointsPlusNormals && m!=Splatting); // CHECK THIS!
    }

    void draw(Viewer_interface* viewer) const {
        if(!are_buffers_filled)
            initialize_buffers(viewer);
        vaos[0]->bind();
        program = getShaderProgram(PROGRAM_WITH_LIGHT);
        attrib_buffers(viewer, PROGRAM_WITH_LIGHT);
        program->bind();
        viewer->glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(positions_poly.size()/3));
        program->release();
        vaos[0]->release();


    }
    void draw_edges(Viewer_interface* viewer) const {
        if(!are_buffers_filled)
            initialize_buffers(viewer);
        vaos[2]->bind();
        program = getShaderProgram(PROGRAM_WITHOUT_LIGHT);
        attrib_buffers(viewer, PROGRAM_WITHOUT_LIGHT);
        program->bind();
        QMatrix4x4 f_mat;
        for(int i=0; i<16; i++)
            f_mat.data()[i]=frame->matrix()[i];
        program->setUniformValue("f_matrix",f_mat);
        viewer->glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(positions_grid.size()/3));
        program->release();
        vaos[2]->release();

        vaos[1]->bind();
        program = getShaderProgram(PROGRAM_WITHOUT_LIGHT);
        attrib_buffers(viewer, PROGRAM_WITHOUT_LIGHT);
        program->bind();
        viewer->glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(positions_lines.size()/3));
        program->release();
        vaos[1]->release();

    }
    void draw_points(Viewer_interface * viewer) const
    {
        if(!are_buffers_filled)
            initialize_buffers(viewer);
        vaos[1]->bind();
        program = getShaderProgram(PROGRAM_WITHOUT_LIGHT);
        attrib_buffers(viewer, PROGRAM_WITHOUT_LIGHT);
        program->bind();
        viewer->glDrawArrays(GL_POINTS, 0, static_cast<GLsizei>(positions_lines.size()/3));
       vaos[1]->release();
       program->release();

       vaos[2]->bind();
       program = getShaderProgram(PROGRAM_WITHOUT_LIGHT);
       attrib_buffers(viewer, PROGRAM_WITHOUT_LIGHT);
       program->bind();
       QMatrix4x4 f_mat;
       for(int i=0; i<16; i++)
           f_mat.data()[i]=frame->matrix()[i];
       program->setUniformValue("f_matrix",f_mat);
       viewer->glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(positions_grid.size()/3));
       program->release();
       vaos[2]->release();
    }
private:
    void draw_triangle(const Kernel::Point_3& pa,
                       const Kernel::Point_3& pb,
                       const Kernel::Point_3& pc, bool is_cut) {

#undef darker
        Kernel::Vector_3 n = cross_product(pb - pa, pc - pa);
        n = n / CGAL::sqrt(n*n);

if(!is_cut)
{
    for(int i=0; i<3; i++)
    {

        color_poly.push_back(this->color().redF());
        color_poly.push_back(this->color().greenF());
        color_poly.push_back(this->color().blueF());
    }
}
else
{
    for(int i=0; i<3; i++)
    {
      color_poly.push_back(this->color().darker(150).redF());
      color_poly.push_back(this->color().darker(150).greenF());
      color_poly.push_back(this->color().darker(150).blueF());
    }
}
        for(int i=0; i<3; i++)
        {
            normals.push_back(n.x());
            normals.push_back(n.y());
            normals.push_back(n.z());
        }
        positions_poly.push_back(pa.x());
        positions_poly.push_back(pa.y());
        positions_poly.push_back(pa.z());

        positions_poly.push_back(pb.x());
        positions_poly.push_back(pb.y());
        positions_poly.push_back(pb.z());

        positions_poly.push_back(pc.x());
        positions_poly.push_back(pc.y());
        positions_poly.push_back(pc.z());

    }

    void draw_triangle_edges(const Kernel::Point_3& pa,
                       const Kernel::Point_3& pb,
                       const Kernel::Point_3& pc) {

#undef darker
        Kernel::Vector_3 n = cross_product(pb - pa, pc - pa);
        n = n / CGAL::sqrt(n*n);
        for(int i=0; i<6; i++)
        {
            color_lines.push_back(0.0);
            color_lines.push_back(0.0);
            color_lines.push_back(0.0);
        }
        positions_lines.push_back(pa.x());
        positions_lines.push_back(pa.y());
        positions_lines.push_back(pa.z());

        positions_lines.push_back(pb.x());
        positions_lines.push_back(pb.y());
        positions_lines.push_back(pb.z());

        positions_lines.push_back(pb.x());
        positions_lines.push_back(pb.y());
        positions_lines.push_back(pb.z());

        positions_lines.push_back(pc.x());
        positions_lines.push_back(pc.y());
        positions_lines.push_back(pc.z());

        positions_lines.push_back(pc.x());
        positions_lines.push_back(pc.y());
        positions_lines.push_back(pc.z());

        positions_lines.push_back(pa.x());
        positions_lines.push_back(pa.y());
        positions_lines.push_back(pa.z());

    }



    double complex_diag() const {
        const Bbox& bbox = this->bbox();
        const double& xdelta = bbox.xmax-bbox.xmin;
        const double& ydelta = bbox.ymax-bbox.ymin;
        const double& zdelta = bbox.zmax-bbox.zmin;
        const double diag = std::sqrt(xdelta*xdelta +
                                      ydelta*ydelta +
                                      zdelta*zdelta);
        return diag * 0.7;
    }

public Q_SLOTS:
    void export_facets_in_complex()
    {
        std::stringstream off_sstream;
        c3t3().output_facets_in_complex_to_off(off_sstream);
        std::string backup = off_sstream.str();
        // Try to read .off in a polyhedron
        Scene_polyhedron_item* item = new Scene_polyhedron_item();
        if(!item->load(off_sstream))
        {
            delete item;
            off_sstream.str(backup);

            // Try to read .off in a polygon soup
            Scene_polygon_soup_item* soup_item = new Scene_polygon_soup_item;

            if(!soup_item->load(off_sstream)) {
                delete soup_item;
                return;
            }

            soup_item->setName(QString("%1_%2").arg(this->name()).arg("facets"));
            last_known_scene->addItem(soup_item);
        }
        else{
            item->setName(QString("%1_%2").arg(this->name()).arg("facets"));
            last_known_scene->addItem(item);
        }
    }

public:

    QMenu* contextMenu()
    {
        const char* prop_name = "Menu modified by Scene_c3t3_item.";

        QMenu* menu = Scene_item::contextMenu();

        // Use dynamic properties:
        // http://doc.qt.io/qt-5/qobject.html#property
        bool menuChanged = menu->property(prop_name).toBool();

        if(!menuChanged) {
            QAction* actionExportFacetsInComplex =
                    menu->addAction(tr("Export facets in complex"));
            actionExportFacetsInComplex->setObjectName("actionExportFacetsInComplex");
            connect(actionExportFacetsInComplex,
                    SIGNAL(triggered()),this,
                    SLOT(export_facets_in_complex()));
        }
        return menu;
    }

    void set_scene(Scene_interface* scene){ last_known_scene=scene; }

private:
    C3t3 c3t3_;
    qglviewer::ManipulatedFrame* frame;
    Scene_interface* last_known_scene;


    std::vector<float> positions_lines;
    std::vector<float> positions_grid;
    std::vector<float> positions_poly;
    std::vector<float> normals;
    std::vector<float> color_lines;
    std::vector<float> color_poly;
    std::vector<float> color_grid;

    mutable QOpenGLShaderProgram *program;

    using Scene_item::initialize_buffers;
    void initialize_buffers(Viewer_interface *viewer)const
    {
        //vao containing the data for the facets
        {
            program = getShaderProgram(PROGRAM_WITH_LIGHT, viewer);
            program->bind();

            vaos[0]->bind();
            buffers[0].bind();
            buffers[0].allocate(positions_poly.data(),
                                static_cast<int>(positions_poly.size()*sizeof(float)));
            program->enableAttributeArray("vertex");
            program->setAttributeBuffer("vertex",GL_FLOAT,0,3);
            buffers[0].release();

            buffers[1].bind();
            buffers[1].allocate(normals.data(),
                                static_cast<int>(normals.size()*sizeof(float)));
            program->enableAttributeArray("normals");
            program->setAttributeBuffer("normals",GL_FLOAT,0,3);
            buffers[1].release();

            buffers[2].bind();
            buffers[2].allocate(color_poly.data(),
                                static_cast<int>(color_poly.size()*sizeof(float)));
            program->enableAttributeArray("colors");
            program->setAttributeBuffer("colors",GL_FLOAT,0,3);
            buffers[2].release();
            vaos[0]->release();
            program->release();

        }

        //vao containing the data for the lines
        {
            program = getShaderProgram(PROGRAM_WITHOUT_LIGHT, viewer);
            program->bind();

            vaos[1]->bind();
            buffers[3].bind();
            buffers[3].allocate(positions_lines.data(),
                                static_cast<int>(positions_lines.size()*sizeof(float)));
            program->enableAttributeArray("vertex");
            program->setAttributeBuffer("vertex",GL_FLOAT,0,3);
            buffers[3].release();

            buffers[4].bind();
            buffers[4].allocate(color_lines.data(),
                                static_cast<int>(color_lines.size()*sizeof(float)));
            program->enableAttributeArray("colors");
            program->setAttributeBuffer("colors",GL_FLOAT,0,3);
            buffers[4].release();
            vaos[1]->release();
            program->release();

        }

        //vao containing the data for the grid
        {
            program = getShaderProgram(PROGRAM_WITHOUT_LIGHT, viewer);
            program->bind();

            vaos[2]->bind();
            buffers[5].bind();
            buffers[5].allocate(positions_grid.data(),
                                static_cast<int>(positions_grid.size()*sizeof(float)));
            program->enableAttributeArray("vertex");
            program->setAttributeBuffer("vertex",GL_FLOAT,0,3);
            buffers[5].release();

            buffers[6].bind();
            buffers[6].allocate(color_grid.data(),
                                static_cast<int>(color_grid.size()*sizeof(float)));
            program->enableAttributeArray("colors");
            program->setAttributeBuffer("colors",GL_FLOAT,0,3);
            buffers[6].release();
            vaos[2]->release();
            program->release();
        }
        are_buffers_filled = true;
    }
    void compute_elements()
    {
        positions_lines.clear();
        positions_poly.clear();
        color_lines.clear();
        color_grid.clear();
        color_poly.clear();
        normals.clear();

        //The grid
        {
            float x = (2*(float)complex_diag())/10.0;
            float y = (2*(float)complex_diag())/10.0;
            for(int u = 0; u < 11; u++)
            {

                positions_grid.push_back(-(float)complex_diag() + x* u);
                positions_grid.push_back(-(float)complex_diag());
                positions_grid.push_back(0.0);

                positions_grid.push_back(-(float)complex_diag() + x* u);
                positions_grid.push_back((float)complex_diag());
                positions_grid.push_back(0.0);
            }
            for(int v=0; v<11; v++)
            {

                positions_grid.push_back(-(float)complex_diag());
                positions_grid.push_back(-(float)complex_diag() + v * y);
                positions_grid.push_back(0.0);

                positions_grid.push_back((float)complex_diag());
                positions_grid.push_back(-(float)complex_diag() + v * y);
                positions_grid.push_back(0.0);
            }
            float colors[3];
            colors[0] = this->color().redF();
            colors[1] = this->color().greenF();
            colors[2] = this->color().blueF();

            for(int i=0; i< 132; i++)
            {
                color_grid.push_back(colors[i%3]);
            }
        }

        //The facets
        {
            if(isEmpty())
                return;

            const Kernel::Plane_3& plane = this->plane();
            GLdouble clip_plane[4];
            clip_plane[0] = -plane.a();
            clip_plane[1] = -plane.b();
            clip_plane[2] = -plane.c();
            clip_plane[3] = -plane.d();



            for(C3t3::Facet_iterator
                fit = c3t3().facets_begin(),
                end = c3t3().facets_end();
                fit != end; ++fit)
            {
                const Tr::Cell_handle& cell = fit->first;
                const int& index = fit->second;
                const Kernel::Point_3& pa = cell->vertex((index+1)&3)->point();
                const Kernel::Point_3& pb = cell->vertex((index+2)&3)->point();
                const Kernel::Point_3& pc = cell->vertex((index+3)&3)->point();
                typedef Kernel::Oriented_side Side;
                using CGAL::ON_ORIENTED_BOUNDARY;
                const Side sa = plane.oriented_side(pa);
                const Side sb = plane.oriented_side(pb);
                const Side sc = plane.oriented_side(pc);
                bool is_showned = false;
                if(pa.x() * clip_plane[0]  + pa.y() * clip_plane[1]  + pa.z() * clip_plane[2] + clip_plane[3]  > 0
                        && pb.x() * clip_plane[0]  + pb.y() * clip_plane[1]  + pb.z() * clip_plane[2] + clip_plane[3]  > 0
                        && pc.x() * clip_plane[0]  + pc.y() * clip_plane[1]  + pc.z() * clip_plane[2] + clip_plane[3]  > 0)
                    is_showned = true;

                if(is_showned && sa != ON_ORIENTED_BOUNDARY &&
                        sb != ON_ORIENTED_BOUNDARY &&
                        sc != ON_ORIENTED_BOUNDARY &&
                        sb == sa && sc == sa )
                {
                    if ( (index%2 == 1) == c3t3().is_in_complex(cell)) draw_triangle(pb, pa, pc, false);
                    else draw_triangle(pa, pb, pc, false);
                    draw_triangle_edges(pa, pb, pc);
                }

            }


            for(Tr::Finite_cells_iterator
                cit = c3t3().triangulation().finite_cells_begin(),
                end = c3t3().triangulation().finite_cells_end();
                cit != end; ++cit)
            {
                if(! c3t3().is_in_complex(cit) )
                    continue;

                const Kernel::Point_3& pa = cit->vertex(0)->point();
                const Kernel::Point_3& pb = cit->vertex(1)->point();
                const Kernel::Point_3& pc = cit->vertex(2)->point();
                const Kernel::Point_3& pd = cit->vertex(3)->point();
                typedef Kernel::Oriented_side Side;
                using CGAL::ON_ORIENTED_BOUNDARY;
                const Side sa = plane.oriented_side(pa);
                const Side sb = plane.oriented_side(pb);
                const Side sc = plane.oriented_side(pc);
                const Side sd = plane.oriented_side(pd);

                if( sa == ON_ORIENTED_BOUNDARY ||
                        sb == ON_ORIENTED_BOUNDARY ||
                        sc == ON_ORIENTED_BOUNDARY ||
                        sd == ON_ORIENTED_BOUNDARY ||
                        sb != sa || sc != sa || sd != sa)
                {
                    draw_triangle(pb,pa,pc, true);
                    draw_triangle(pa,pb,pd, true);
                    draw_triangle(pa,pd,pc, true);
                    draw_triangle(pb,pc,pd, true);

                    draw_triangle_edges(pa,pb,pc);
                    draw_triangle_edges(pa,pb,pd);
                    draw_triangle_edges(pa,pc,pd);
                    draw_triangle_edges(pb,pc,pd);
                }
            }
        }
    }

};

Scene_item* cgal_code_mesh_3(const Polyhedron* pMesh,
                             QString filename,
                             const double angle,
                             const double facet_sizing,
                             const double approx,
                             const double tet_sizing,
                             const double tet_shape,
                             const bool protect_features,
                             Scene_interface* scene)
{
  if(!pMesh) return 0;

  // remesh

  // Set mesh criteria
  Edge_criteria edge_criteria(facet_sizing);
  Facet_criteria facet_criteria(angle, facet_sizing, approx); // angle, size, approximation
  Cell_criteria cell_criteria(tet_shape, tet_sizing); // radius-edge ratio, size
  Mesh_criteria criteria(edge_criteria, facet_criteria, cell_criteria);

  CGAL::Timer timer;
  timer.start();
  std::cerr << "Meshing file \"" << qPrintable(filename) << "\"\n";
  std::cerr << "  angle: " << angle << std::endl
            << "  facets size bound: " << facet_sizing << std::endl
            << "  approximation bound: " << approx << std::endl
            << "  tetrahedra size bound: " << tet_sizing << std::endl;
  std::cerr << "Build AABB tree...";
  // Create domain
  Mesh_domain domain(*pMesh);
  if(protect_features) {
      domain.detect_features();
  }
  std::cerr << "done (" << timer.time() << " ms)" << std::endl;

  // Meshing
  std::cerr << "Mesh...";
  CGAL::parameters::internal::Features_options features =
          protect_features ?
              CGAL::parameters::features(domain) :
              CGAL::parameters::no_features();

  Scene_c3t3_item* new_item =
          new Scene_c3t3_item(CGAL::make_mesh_3<C3t3>(domain, criteria, features));
  new_item->set_scene(scene);
  std::cerr << "done (" << timer.time() << " ms, " << new_item->c3t3().triangulation().number_of_vertices() << " vertices)" << std::endl;

  if(new_item->c3t3().triangulation().number_of_vertices() > 0)
  {
    std::ofstream medit_out("out.mesh");
    new_item->c3t3().output_to_medit(medit_out);

    const Scene_item::Bbox& bbox = new_item->bbox();
    new_item->setPosition((float)(bbox.xmin + bbox.xmax)/2.f,
                          (float)(bbox.ymin + bbox.ymax)/2.f,
                          (float)(bbox.zmin + bbox.zmax)/2.f);
    return new_item;
  }
  else {
    delete new_item;
    return 0;
  }
}

#include "Polyhedron_demo_mesh_3_plugin_cgal_code.moc"
//#include "Scene_c3t3_item.moc" //Check this one, it's strange moc include.

