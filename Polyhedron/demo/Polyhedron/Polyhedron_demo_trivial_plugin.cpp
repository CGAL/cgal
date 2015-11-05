#include <QtCore/qglobal.h>

#include "Scene_item.h"
#include "Scene_interface.h"
#include <CGAL/gl.h>

#include "Viewer_interface.h"
#include <QAction>
#include <QMainWindow>

class Q_DECL_EXPORT Scene_bbox_item : public Scene_item
{
    Q_OBJECT
public:
    Scene_bbox_item(const Scene_interface* scene_interface)
        :  Scene_item(1,1), scene(scene_interface)

    {

        positions_lines.resize(0);
        //Generates an integer which will be used as ID for each buffer
    }
    ~Scene_bbox_item()
    {
    }
    bool isFinite() const { return true; }
    bool isEmpty() const { return true; }
    Bbox bbox() const { return Bbox(); }

    Scene_bbox_item* clone() const {
        return 0;
    }

    QString toolTip() const {
        const Bbox& bb = scene->bbox();
        return QString("<p><b>Scene bounding box</b></p>"
                       "<p>x range: (%1, %2)<br />"
                       "y range: (%3, %4)<br />"
                       "z range: (%5, %6)</p>")
                .arg(bb.xmin).arg(bb.xmax)
                .arg(bb.ymin).arg(bb.ymax)
                .arg(bb.zmin).arg(bb.zmax);
    }

    // Indicate if rendering mode is supported
    bool supportsRenderingMode(RenderingMode m) const {
        return (m == Wireframe);
    }

    void draw_edges(Viewer_interface* viewer) const
    {
        if(!are_buffers_filled)
            initialize_buffers(viewer);
        vaos[0]->bind();
        program = getShaderProgram(PROGRAM_WITHOUT_LIGHT);
        attrib_buffers(viewer, PROGRAM_WITHOUT_LIGHT);
        program->bind();
        program->setAttributeValue("colors", this->color());
        viewer->glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(positions_lines.size()/3));
        vaos[0]->release();
        program->release();

    }

    void invalidate_buffers()
    {
        compute_elements();
        are_buffers_filled = false;
    }

private:

    std::vector<float> positions_lines;
    mutable QOpenGLShaderProgram *program;
    using Scene_item::initialize_buffers;
    void initialize_buffers(Viewer_interface *viewer)const
    {

        //vao containing the data for the lines
        {
            program = getShaderProgram(PROGRAM_WITHOUT_LIGHT, viewer);
            program->bind();

            vaos[0]->bind();
            buffers[0].bind();
            buffers[0].allocate(positions_lines.data(),
                                static_cast<GLsizei>(positions_lines.size()*sizeof(float)));
            program->enableAttributeArray("vertex");
            program->setAttributeBuffer("vertex",GL_FLOAT,0,3);
            buffers[0].release();

            vaos[0]->release();
            program->release();

        }
        are_buffers_filled = true;
    }

    void compute_elements()
    {
        positions_lines.clear();
        const Bbox& bb = scene->bbox();
        positions_lines.push_back(bb.xmin); positions_lines.push_back(bb.ymin); positions_lines.push_back(bb.zmin);
        positions_lines.push_back(bb.xmax); positions_lines.push_back(bb.ymin); positions_lines.push_back(bb.zmin);
        positions_lines.push_back(bb.xmin); positions_lines.push_back(bb.ymin); positions_lines.push_back(bb.zmin);
        positions_lines.push_back(bb.xmin); positions_lines.push_back(bb.ymax); positions_lines.push_back(bb.zmin);
        positions_lines.push_back(bb.xmin); positions_lines.push_back(bb.ymin); positions_lines.push_back(bb.zmin);
        positions_lines.push_back(bb.xmin); positions_lines.push_back(bb.ymin); positions_lines.push_back(bb.zmax);

        positions_lines.push_back(bb.xmax); positions_lines.push_back(bb.ymin); positions_lines.push_back(bb.zmin);
        positions_lines.push_back(bb.xmax); positions_lines.push_back(bb.ymax); positions_lines.push_back(bb.zmin);
        positions_lines.push_back(bb.xmax); positions_lines.push_back(bb.ymin); positions_lines.push_back(bb.zmin);
        positions_lines.push_back(bb.xmax); positions_lines.push_back(bb.ymin); positions_lines.push_back(bb.zmax);

        positions_lines.push_back(bb.xmin); positions_lines.push_back(bb.ymax); positions_lines.push_back(bb.zmin);
        positions_lines.push_back(bb.xmax); positions_lines.push_back(bb.ymax); positions_lines.push_back(bb.zmin);
        positions_lines.push_back(bb.xmin); positions_lines.push_back(bb.ymax); positions_lines.push_back(bb.zmin);
        positions_lines.push_back(bb.xmin); positions_lines.push_back(bb.ymax); positions_lines.push_back(bb.zmax);

        positions_lines.push_back(bb.xmin); positions_lines.push_back(bb.ymin); positions_lines.push_back(bb.zmax);
        positions_lines.push_back(bb.xmax); positions_lines.push_back(bb.ymin); positions_lines.push_back(bb.zmax);
        positions_lines.push_back(bb.xmin); positions_lines.push_back(bb.ymin); positions_lines.push_back(bb.zmax);
        positions_lines.push_back(bb.xmin); positions_lines.push_back(bb.ymax); positions_lines.push_back(bb.zmax);

        positions_lines.push_back(bb.xmax); positions_lines.push_back(bb.ymax); positions_lines.push_back(bb.zmax);
        positions_lines.push_back(bb.xmin); positions_lines.push_back(bb.ymax); positions_lines.push_back(bb.zmax);
        positions_lines.push_back(bb.xmax); positions_lines.push_back(bb.ymax); positions_lines.push_back(bb.zmax);
        positions_lines.push_back(bb.xmax); positions_lines.push_back(bb.ymin); positions_lines.push_back(bb.zmax);
        positions_lines.push_back(bb.xmax); positions_lines.push_back(bb.ymax); positions_lines.push_back(bb.zmax);
        positions_lines.push_back(bb.xmax); positions_lines.push_back(bb.ymax); positions_lines.push_back(bb.zmin);
    }

    const Scene_interface* scene;
};

#include "Polyhedron_demo_plugin_interface.h"

class Polyhedron_demo_trivial_plugin : 
        public QObject,
        public Polyhedron_demo_plugin_interface
{
    Q_OBJECT
    Q_INTERFACES(Polyhedron_demo_plugin_interface)
    Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

public:
    void init(QMainWindow* mainWindow, Scene_interface* scene_interface);
    QList<QAction*> actions() const {
        return QList<QAction*>() << actionBbox;
    }

    bool applicable(QAction*) const {
        return true;
    }
public Q_SLOTS:

    void bbox();
    void enableAction();

private:
    Scene_interface* scene;
    QAction* actionBbox;


}; // end Polyhedron_demo_trivial_plugin

void Polyhedron_demo_trivial_plugin::init(QMainWindow* mainWindow, Scene_interface* scene_interface)
{
    scene = scene_interface;
    actionBbox = new QAction(tr("Create bbox"), mainWindow);
    connect(actionBbox, SIGNAL(triggered()),
            this, SLOT(bbox()));
}

void Polyhedron_demo_trivial_plugin::bbox()
{
    for(int i = 0, end = scene->numberOfEntries();
        i < end; ++i)
    {
        if(qobject_cast<Scene_bbox_item*>(scene->item(i)))
            return;
    }
    Scene_item* item = new Scene_bbox_item(scene);
    connect(item, SIGNAL(destroyed()),
            this, SLOT(enableAction()));
    item->setName("Scene bbox");
    item->setColor(Qt::black);
    item->setRenderingMode(Wireframe);
    scene->addItem(item);
    actionBbox->setEnabled(false);
}

void Polyhedron_demo_trivial_plugin::enableAction() {
    actionBbox->setEnabled(true);
}

#include "Polyhedron_demo_trivial_plugin.moc"
