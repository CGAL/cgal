#include "Scene_polyhedron_transform_item.h"
#include "Kernel_type.h"
#include "Polyhedron_type.h"

Scene_polyhedron_transform_item::Scene_polyhedron_transform_item(const qglviewer::Vec& pos,const Scene_polyhedron_item* poly_item_,const Scene_interface*):
    Scene_item(1,1),
    poly_item(poly_item_),
    manipulable(false),
    frame(new ManipulatedFrame()),
    poly(poly_item->polyhedron()),
    center_(pos)

{
    frame->setPosition(pos);
    nb_lines = 0;
}

void Scene_polyhedron_transform_item::initialize_buffers(Viewer_interface *viewer =0) const
{
    //vao for the edges
    {
        program = getShaderProgram(PROGRAM_WITHOUT_LIGHT, viewer);
        program->bind();

        vaos[0]->bind();
        buffers[0].bind();
        buffers[0].allocate(positions_lines.data(),
                            static_cast<int>(positions_lines.size()*sizeof(float)));
        program->enableAttributeArray("vertex");
        program->setAttributeBuffer("vertex",GL_FLOAT,0,3);
        buffers[0].release();

        QColor color = this->color();
        program->setAttributeValue("colors",color);
        vaos[0]->release();
        program->release();
    }
    nb_lines = positions_lines.size();
    positions_lines.resize(0);
    std::vector<float>(positions_lines).swap(positions_lines);

    are_buffers_filled = true;
}

void Scene_polyhedron_transform_item::compute_elements()
{
     positions_lines.resize(0);
    typedef Kernel::Point_3		Point;
    typedef Polyhedron::Edge_const_iterator	Edge_iterator;

    Edge_iterator he;
    for(he = poly->edges_begin();
        he != poly->edges_end();
        he++)
    {
        const Point& a = he->vertex()->point();
        const Point& b = he->opposite()->vertex()->point();
        positions_lines.push_back(a.x()-center_.x);
        positions_lines.push_back(a.y()-center_.y);
        positions_lines.push_back(a.z()-center_.z);

        positions_lines.push_back(b.x()-center_.x);
        positions_lines.push_back(b.y()-center_.y);
        positions_lines.push_back(b.z()-center_.z);

    }

}

void Scene_polyhedron_transform_item::draw_edges(Viewer_interface* viewer) const
{
    if(!are_buffers_filled)
        initialize_buffers(viewer);
    vaos[0]->bind();
    program = getShaderProgram(PROGRAM_WITHOUT_LIGHT);
    attrib_buffers(viewer,PROGRAM_WITHOUT_LIGHT);
    program->bind();
    QMatrix4x4 f_matrix;
    for (int i=0; i<16; ++i){
        f_matrix.data()[i] = (float)frame->matrix()[i];
    }
    program->setUniformValue("f_matrix", f_matrix);
    viewer->glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(nb_lines/3));
    vaos[0]->release();
    program->release();

}

QString Scene_polyhedron_transform_item::toolTip() const {
    return QObject::tr("<p>Affine transformation of <b>%1</b></p>"
                       "<p>Keep <b>Ctrl</b> pressed and use the arcball to define an affine transformation.<br />"
                       "Press <b>S</b> to apply the affine transformation to a copy of <b>%1</b>.</p>")
            .arg(getBase()->name());
}
bool Scene_polyhedron_transform_item::keyPressEvent(QKeyEvent* e){
    if (e->key()==Qt::Key_S){
    Q_EMIT stop();
        return true;
    }
    return false;
}

Scene_polyhedron_transform_item::Bbox
Scene_polyhedron_transform_item::bbox() const {
    const Kernel::Point_3& p = *(poly->points_begin());
    CGAL::Bbox_3 bbox(p.x(), p.y(), p.z(), p.x(), p.y(), p.z());
    for(Polyhedron::Point_const_iterator it = poly->points_begin();
        it != poly->points_end();
        ++it) {
        bbox = bbox + it->bbox();
    }
    return Bbox(bbox.xmin(),bbox.ymin(),bbox.zmin(),
                bbox.xmax(),bbox.ymax(),bbox.zmax());
}


void Scene_polyhedron_transform_item::invalidate_buffers()
{
    compute_elements();
    are_buffers_filled = false;
}

