#include "Scene_polylines_item.h"
#include "create_sphere.h"

#include <CGAL/bounding_box.h>
#include <CGAL/gl.h>
#include <QMenu>
#include <QAction>

#include <QInputDialog>

class Scene_polylines_item_private {
public:
    typedef Scene_polylines_item::K K;
    typedef K::Point_3 Point_3;

    Scene_polylines_item_private() :
        draw_extremities(false),
        spheres_drawn_radius(0)
    {}

    void draw_sphere(const K::Point_3&, double) const;
    void draw_spheres(const Scene_polylines_item*) const;

    bool draw_extremities;
    double spheres_drawn_radius;
};

void
Scene_polylines_item::create_Sphere(float R) const
{
  create_flat_and_wire_sphere(R, positions_spheres, normals_spheres, positions_wire_spheres);
}

void
Scene_polylines_item::initialize_buffers(CGAL::Three::Viewer_interface *viewer = 0) const
{
    QOpenGLShaderProgram *program;
   //vao for the lines
    {
        program = getShaderProgram(PROGRAM_NO_SELECTION, viewer);
        program->bind();

        vaos[Edges]->bind();
        buffers[Edges_Vertices].bind();
        buffers[Edges_Vertices].allocate(positions_lines.data(),
                            static_cast<int>(positions_lines.size()*sizeof(float)));
        program->enableAttributeArray("vertex");
        program->setAttributeBuffer("vertex",GL_FLOAT,0,4);
        buffers[Edges_Vertices].release();
        vaos[Edges]->release();
        program->release();
    }
   //vao for the spheres
    {
        if(viewer->extension_is_found)
        {
            program = getShaderProgram(PROGRAM_INSTANCED, viewer);
            program->bind();

            vaos[Spheres]->bind();
            buffers[Spheres_Vertices].bind();
            buffers[Spheres_Vertices].allocate(positions_spheres.data(),
                                static_cast<int>(positions_spheres.size()*sizeof(float)));
            program->enableAttributeArray("vertex");
            program->setAttributeBuffer("vertex",GL_FLOAT,0,3);
            buffers[Spheres_Vertices].release();

            buffers[Spheres_Normals].bind();
            buffers[Spheres_Normals].allocate(normals_spheres.data(),
                                static_cast<int>(normals_spheres.size()*sizeof(float)));
            program->enableAttributeArray("normals");
            program->setAttributeBuffer("normals",GL_FLOAT,0,3);
            buffers[Spheres_Normals].release();

            buffers[Spheres_Colors].bind();
            buffers[Spheres_Colors].allocate(color_spheres.data(),
                                static_cast<int>(color_spheres.size()*sizeof(float)));
            program->enableAttributeArray("colors");
            program->setAttributeBuffer("colors",GL_FLOAT,0,3);
            buffers[Spheres_Colors].release();

            buffers[Spheres_Center].bind();
            buffers[Spheres_Center].allocate(positions_center.data(),
                                static_cast<int>(positions_center.size()*sizeof(float)));
            program->enableAttributeArray("center");
            program->setAttributeBuffer("center",GL_FLOAT,0,3);
            buffers[Spheres_Center].release();

            viewer->glVertexAttribDivisor(program->attributeLocation("center"), 1);
            viewer->glVertexAttribDivisor(program->attributeLocation("colors"), 1);

        }
        else
        {
            program = getShaderProgram(PROGRAM_NO_SELECTION, viewer);
            program->bind();

            vaos[Spheres]->bind();
            buffers[Spheres_Vertices].bind();
            buffers[Spheres_Vertices].allocate(positions_center.data(),
                                static_cast<int>(positions_center.size()*sizeof(float)));
            program->enableAttributeArray("vertex");
            program->setAttributeBuffer("vertex",GL_FLOAT,0,3);
            buffers[Spheres_Vertices].release();

            buffers[Spheres_Normals].bind();
            buffers[Spheres_Normals].allocate(color_spheres.data(),
                                static_cast<int>(color_spheres.size()*sizeof(float)));
            program->enableAttributeArray("colors");
            program->setAttributeBuffer("colors",GL_FLOAT,0,3);
            buffers[Spheres_Normals].release();
        }
        vaos[Spheres]->release();

        program->release();
    }

//vao for the wired spheres
    {
        if(viewer->extension_is_found)
        {
            program = getShaderProgram(PROGRAM_INSTANCED_WIRE, viewer);
            program->bind();

            vaos[Wired_Spheres]->bind();
            buffers[Wired_Spheres_Vertices].bind();
            buffers[Wired_Spheres_Vertices].allocate(positions_wire_spheres.data(),
                                static_cast<int>(positions_wire_spheres.size()*sizeof(float)));
            program->enableAttributeArray("vertex");
            program->setAttributeBuffer("vertex",GL_FLOAT,0,3);
            buffers[Wired_Spheres_Vertices].release();

            buffers[Spheres_Colors].bind();
            program->enableAttributeArray("colors");
            program->setAttributeBuffer("colors",GL_FLOAT,0,3);
            buffers[Spheres_Colors].release();

            buffers[Spheres_Normals].bind();
            program->enableAttributeArray("normals");
            program->setAttributeBuffer("normals",GL_FLOAT,0,3);
            buffers[Spheres_Normals].release();


            buffers[Spheres_Center].bind();
            program->enableAttributeArray("center");
            program->setAttributeBuffer("center",GL_FLOAT,0,3);
            buffers[Spheres_Center].release();


            viewer->glVertexAttribDivisor(program->attributeLocation("center"), 1);
            viewer->glVertexAttribDivisor(program->attributeLocation("colors"), 1);

            vaos[Wired_Spheres]->release();
            program->release();

            nb_lines = positions_lines.size();
            positions_lines.resize(0);
            std::vector<float>(positions_lines).swap(positions_lines);
            nb_spheres = positions_spheres.size();
            positions_spheres.resize(0);
            std::vector<float>(positions_spheres).swap(positions_spheres);
            normals_spheres.resize(0);
            std::vector<float>(normals_spheres).swap(normals_spheres);
            color_spheres.resize(0);
            std::vector<float>(color_spheres).swap(color_spheres);
            nb_centers = positions_center.size();
            positions_center.resize(0);
            std::vector<float>(positions_center).swap(positions_center);
            nb_wire = positions_wire_spheres.size();
            positions_wire_spheres.resize(0);
            std::vector<float>(positions_wire_spheres).swap(positions_wire_spheres);
        }
    }

   are_buffers_filled = true;

}
void
Scene_polylines_item::compute_elements() const
{
    positions_spheres.resize(0);
    positions_wire_spheres.resize(0);
    positions_lines.resize(0);
    color_spheres.resize(0);
    normals_spheres.resize(0);
    positions_center.resize(0);
    nbSpheres = 0;

    //Fills the VBO with the lines
    for(std::list<std::vector<Point_3> >::const_iterator it = polylines.begin();
        it != polylines.end();
        ++it){
        if(it->empty()) continue;
        for(size_t i = 0, end = it->size()-1;
            i < end; ++i)
        {
            const Point_3& a = (*it)[i];
            const Point_3& b = (*it)[i+1];
            positions_lines.push_back(a.x());
            positions_lines.push_back(a.y());
            positions_lines.push_back(a.z());
            positions_lines.push_back(1.0);

            positions_lines.push_back(b.x());
            positions_lines.push_back(b.y());
            positions_lines.push_back(b.z());
            positions_lines.push_back(1.0);

        }

    }
    //Fills the VBO with the spheres
    if(d->draw_extremities)
    {

        // FIRST, count the number of incident cycles and polylines
        // for all extremities.
        typedef std::map<Point_3, int> Point_to_int_map;
        typedef Point_to_int_map::iterator iterator;
        Point_to_int_map corner_polyline_nb;

        { // scope to fill corner_polyline_nb'
            Point_to_int_map corner_cycles_nb;

            for(std::list<std::vector<Point_3> >::const_iterator
                it = this->polylines.begin(),
                end = this->polylines.end();
                it != end; ++it)
            {
                const K::Point_3& a = *it->begin();
                const K::Point_3& b = *it->rbegin();
                if(a == b) {
                    if ( it->size()>1 )
                        ++corner_cycles_nb[a];
                    else
                        ++corner_polyline_nb[a];
                }
                else {
                    ++corner_polyline_nb[a];
                    ++corner_polyline_nb[b];
                }
            }
            // THEN, ignore points that are incident to one cycle only.
            for(iterator
                c_it = corner_cycles_nb.begin(),
                end = corner_cycles_nb.end();
                c_it != end; ++c_it)
            {
                const Point_3& a = c_it->first;

                iterator p_it = corner_polyline_nb.find(a);

                // If the point 'a'=c_it->first has only incident cycles...
                if(p_it == corner_polyline_nb.end()) {
                    // ...then count it as a corner only if it has two incident cycles
                    // or more.
                    if(c_it->second > 1) {
                        corner_polyline_nb[a] = c_it->second;
                    }
                } else {
                    // else add the number of cycles.
                    p_it->second += c_it->second;
                }
            }
        }
        // At this point, 'corner_polyline_nb' gives the multiplicity of all
        // corners.
        //Finds the centers of the spheres and their color
        for(iterator
            p_it = corner_polyline_nb.begin(),
            end = corner_polyline_nb.end();
            p_it != end; ++p_it)
        {
            nbSpheres++;
            const K::Point_3& centre = p_it->first;
            positions_center.push_back(centre.x());
            positions_center.push_back(centre.y());
            positions_center.push_back(centre.z());

            float colors[3];
            switch(p_it->second) {
            case 1:
                colors[0] = 0.0; // black
                colors[1] = 0.0;
                colors[2] = 0.0;
                break;
            case 2:
                colors[0] = 0.0; // green
                colors[1] = 0.8f;
                colors[2] = 0.0;
                break;
            case 3:
                colors[0] = 0.0; // blue
                colors[1] = 0.0;
                colors[2] = 0.8f;
                break;
            case 4:
                colors[0] = 0.8f; //red
                colors[1] = 0.0;
                colors[2] = 0.0;
                break;
            default:
                colors[0] = 0.8f; //fuschia
                colors[1] = 0.0;
                colors[2] = 0.8f;
            }

            color_spheres.push_back(colors[0]);
            color_spheres.push_back(colors[1]);
            color_spheres.push_back(colors[2]);
            color_wire_spheres.push_back(colors[0]);
            color_wire_spheres.push_back(colors[1]);
            color_wire_spheres.push_back(colors[2]);
            color_wire_spheres.push_back(colors[0]);
            color_wire_spheres.push_back(colors[1]);
            color_wire_spheres.push_back(colors[2]);
        }
        create_Sphere(d->spheres_drawn_radius);

    }
}


Scene_polylines_item::Scene_polylines_item() 
    :CGAL::Three::Scene_item(NbOfVbos,NbOfVaos)
    ,d(new Scene_polylines_item_private())
    ,nbSpheres(0)
{
    setRenderingMode(FlatPlusEdges);
    nb_spheres = 0;
    nb_wire = 0;
    nb_centers = 0;
    nb_lines = 0;
    invalidateOpenGLBuffers();

}

Scene_polylines_item::~Scene_polylines_item()
{
    delete d;

}

bool
Scene_polylines_item::isEmpty() const {
    return polylines.empty();
}

void
Scene_polylines_item::compute_bbox() const {
    typedef K::Iso_cuboid_3 Iso_cuboid_3;

    if(isEmpty())
    {
        _bbox =Bbox();
        return;
    }
    std::list<Point_3> boxes;
    for(std::list<std::vector<Point_3> >::const_iterator it = polylines.begin();
        it != polylines.end();
        ++it){
        if(it->begin() != it->end()) {
            Iso_cuboid_3 cub = CGAL::bounding_box(it->begin(), it->end());
            boxes.push_back((cub.min)());
            boxes.push_back((cub.max)());
        }
    }
    Iso_cuboid_3 bbox =
            boxes.begin() != boxes.end() ?
                CGAL::bounding_box(boxes.begin(), boxes.end()) :
                Iso_cuboid_3();

    _bbox = Bbox(bbox.xmin(),
                bbox.ymin(),
                bbox.zmin(),
                bbox.xmax(),
                bbox.ymax(),
                bbox.zmax());
}

Scene_polylines_item* 
Scene_polylines_item::clone() const {
    Scene_polylines_item* item = new Scene_polylines_item;
    item->polylines = polylines;
    QVariant metadata_variant = property("polylines metadata");
    if(metadata_variant.type() == QVariant::StringList)
    {
        item->setProperty("polylines metadata", metadata_variant);
    }
    return item;
}

QString
Scene_polylines_item::toolTip() const {
    QString s =
            tr("<p><b>%1</b> (mode: %2, color: %3)<br />"
               "<i>Polylines</i></p>"
               "<p>Number of polylines: %4</p>")
            .arg(this->name())
            .arg(this->renderingModeName())
            .arg(this->color().name())
            .arg(polylines.size());
    if(d->draw_extremities) {
        s += tr("<p>Legende of endpoints colors: <ul>"
                "<li>black: one incident polyline</li>"
                "<li>green: two incident polylines</li>"
                "<li>blue: three incident polylines</li>"
                "<li>red: four incident polylines</li>"
                "<li>fuchsia: five or more incident polylines</li>"
                "</ul></p>");
    }
    return s;
}

bool
Scene_polylines_item::supportsRenderingMode(RenderingMode m) const {
    return (m == Wireframe ||
            m == FlatPlusEdges ||
            m == Points);
}

// Shaded OpenGL drawing: only draw spheres
void
Scene_polylines_item::draw(CGAL::Three::Viewer_interface* viewer) const {

    if(!are_buffers_filled)
    {
        compute_elements();
        initialize_buffers(viewer);
    }
    if(d->draw_extremities)
    {
        if(viewer->extension_is_found)
        {
            vaos[Spheres]->bind();
            QOpenGLShaderProgram* program = getShaderProgram(PROGRAM_INSTANCED);
            attrib_buffers(viewer, PROGRAM_INSTANCED);
            program->bind();
            viewer->glDrawArraysInstanced(GL_TRIANGLES, 0,
                                          static_cast<GLsizei>(nb_spheres/3), nbSpheres);
            program->release();
            vaos[Spheres]->release();
        }
        else
        {
            vaos[Spheres]->bind();
            QOpenGLShaderProgram* program = getShaderProgram(PROGRAM_NO_SELECTION);
            attrib_buffers(viewer, PROGRAM_NO_SELECTION);
            glPointSize(8.0f);
            glEnable(GL_POINT_SMOOTH);
            program->bind();
            viewer->glDrawArrays(GL_POINTS, 0, static_cast<GLsizei>(nb_centers/3));
            glDisable(GL_POINT_SMOOTH);
            program->release();
            vaos[Spheres]->release();
        }
    }
}

// Wireframe OpenGL drawing
void 
Scene_polylines_item::draw_edges(CGAL::Three::Viewer_interface* viewer) const {
    if(!are_buffers_filled)
    {
        compute_elements();
        initialize_buffers(viewer);
    }

    vaos[Edges]->bind();
    attrib_buffers(viewer, PROGRAM_NO_SELECTION);
    QOpenGLShaderProgram *program = getShaderProgram(PROGRAM_NO_SELECTION);
    program->bind();
    program->setAttributeValue("colors", this->color());
    viewer->glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(nb_lines/4));
    program->release();
    vaos[Edges]->release();
    if(d->draw_extremities)
    {
        if(viewer->extension_is_found)
        {
            vaos[Wired_Spheres]->bind();
            attrib_buffers(viewer, PROGRAM_INSTANCED_WIRE);
            program = getShaderProgram(PROGRAM_INSTANCED_WIRE);
            program->bind();
            viewer->glDrawArraysInstanced(GL_LINES, 0,
                                          static_cast<GLsizei>(nb_wire/3), nbSpheres);
            program->release();
            vaos[Wired_Spheres]->release();
        }
    }

}

void 
Scene_polylines_item::draw_points(CGAL::Three::Viewer_interface* viewer) const {
    if(!are_buffers_filled)
    {
        compute_elements();
        initialize_buffers(viewer);
    }

    vaos[Edges]->bind();
    attrib_buffers(viewer, PROGRAM_NO_SELECTION);
    QOpenGLShaderProgram *program = getShaderProgram(PROGRAM_NO_SELECTION);
    program->bind();
    QColor temp = this->color();
    program->setAttributeValue("colors", temp);
    viewer->glDrawArrays(GL_POINTS, 0, static_cast<GLsizei>(nb_lines/4));
    // Clean-up
   vaos[Edges]->release();
   program->release();
}

QMenu* Scene_polylines_item::contextMenu() 
{
    const char* prop_name = "Menu modified by Scene_polylines_item.";

    QMenu* menu = Scene_item::contextMenu();

    // Use dynamic properties:
    // http://doc.qt.io/qt-5/qobject.html#property
    bool menuChanged = menu->property(prop_name).toBool();

    if(!menuChanged) {
        menu->addSeparator();
        // TODO: add actions to display corners
        QAction* action = menu->addAction(tr("Display corners with radius..."));
        connect(action, SIGNAL(triggered()),
                this, SLOT(change_corner_radii()));

        QAction* actionSmoothPolylines =
                menu->addAction(tr("Smooth polylines"));
        actionSmoothPolylines->setObjectName("actionSmoothPolylines");
        connect(actionSmoothPolylines, SIGNAL(triggered()),this, SLOT(smooth()));
        menu->setProperty(prop_name, true);
    }
    return menu;
}

void Scene_polylines_item::invalidateOpenGLBuffers()
{
    are_buffers_filled = false;
    compute_bbox();


}

void Scene_polylines_item::change_corner_radii() {
    bool ok = true;
    double proposed_radius = d->spheres_drawn_radius;
    if(proposed_radius == 0) {
        CGAL::Three::Scene_interface::Bbox b = bbox();
        proposed_radius = (std::max)(b.xmax - b.xmin,
                                     proposed_radius);
        proposed_radius = (std::max)(b.ymax - b.ymin,
                                     proposed_radius);
        proposed_radius = (std::max)(b.zmax - b.zmin,
                                     proposed_radius);
        proposed_radius /= 100;
    }
    double r = QInputDialog::getDouble(NULL,
                                       tr("Display corners with new radius..."),
                                       tr("Radius:"),
                                       proposed_radius, // value
                                       0.,          // min
                                       2147483647., // max
                                       10,          // decimals
                                       &ok);
    if(ok) {
        change_corner_radii(r);
    }
}

void Scene_polylines_item::change_corner_radii(double r) {
    if(r >= 0) {
        d->spheres_drawn_radius = r;
        d->draw_extremities = (r > 0);
        this->invalidateOpenGLBuffers();
    Q_EMIT itemChanged();
    }
}

void Scene_polylines_item::split_at_sharp_angles()
{
    typedef Polylines_container Bare_polyline_container;
    typedef Polyline Bare_polyline;
    Polylines_container& bare_polylines = polylines;

    int counter = 0;
    for(Bare_polyline_container::iterator
        bare_polyline_it = bare_polylines.begin();
        bare_polyline_it != bare_polylines.end(); // the end changes
        // during the loop
        ++counter /* bare_polyline_it is incremented in the loop */)
    {
        Bare_polyline_container::iterator current_polyline_it =
                bare_polyline_it;
        Bare_polyline& bare_polyline = *bare_polyline_it;
        Bare_polyline::iterator it = boost::next(bare_polyline.begin());

        if(boost::next(bare_polyline.begin()) == bare_polyline.end())
        {
            std::cerr << "WARNING: Isolated point in polylines\n";
            bare_polyline_it = bare_polylines.erase(bare_polyline_it);
            continue;
        }
        else
            ++bare_polyline_it;
        if(it != bare_polyline.end()) {
            for(; it != boost::prior(bare_polyline.end()); ++it) {
                const Point_3 pv = *it;
                const Point_3 pa = *boost::prior(it);
                const Point_3 pb = *boost::next(it);
                const K::Vector_3 av = pv - pa;
                const K::Vector_3 bv = pv - pb;
                const K::FT sc_prod = av * bv;
                if( sc_prod >= 0 ||
                        (sc_prod < 0 &&
                         CGAL::square(sc_prod) < (av * av) * (bv * bv) / 4 ) )
                {
#ifdef PROTECTION_DEBUG
                    std::cerr << "Split polyline (small angle) "
                              <<  std::acos(sqrt(CGAL::square(sc_prod) /
                                                 ((av*av) * (bv*bv)))) * 180 /CGAL_PI
                               << " degres\n";
#endif
                    Bare_polyline new_polyline;
                    std::copy(it, bare_polyline.end(),
                              std::back_inserter(new_polyline));

                    if(*bare_polyline.begin() == *bare_polyline.rbegin()) {
                        // if the polyline is a cycle, test if its beginning is a sharp
                        // angle...
                        const Point_3 pv = *bare_polyline.begin();
                        const Point_3 pa = *boost::prior(boost::prior(bare_polyline.end()));
                        const Point_3 pb = *boost::next(bare_polyline.begin());
                        const K::Vector_3 av = pv - pa;
                        const K::Vector_3 bv = pv - pb;
                        const K::FT sc_prod = av * bv;
                        if( sc_prod >= 0 ||
                                (sc_prod < 0 &&
                                 CGAL::square(sc_prod) < (av * av) * (bv * bv) / 4 ) )
                        {
                            // if its beginning is a sharp angle, then split
                            bare_polyline.erase(boost::next(it), bare_polyline.end());
                        }
                        else {
                            // ...if not, modifies its beginning
                            std::copy(boost::next(bare_polyline.begin()),
                                      boost::next(it),
                                      std::back_inserter(new_polyline));
                            bare_polylines.erase(current_polyline_it);
                        }
                    }
                    else {
                        bare_polyline.erase(boost::next(it), bare_polyline.end());
                    }
                    bare_polylines.push_back(new_polyline);
                    break;
                }
            }
        }
    }
  Q_EMIT itemChanged();
}

void
Scene_polylines_item::merge(Scene_polylines_item* other_item) {
    if(other_item == 0) return;
    std::copy(other_item->polylines.begin(),
              other_item->polylines.end(),
              std::back_inserter(polylines));
    QVariant other_metadata_variant = other_item->property("polylines metadata");
    if(other_metadata_variant.type() == QVariant::StringList)
    {
        QStringList metadata = property("polylines metadata").toStringList();
        metadata.append(other_metadata_variant.toStringList());
        setProperty("polylines metadata", metadata);
    }
    invalidateOpenGLBuffers();
}

