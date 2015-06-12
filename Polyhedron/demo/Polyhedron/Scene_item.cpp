#include "Scene_item.h"
#include "Scene_interface.h"
#include <QMenu>
#include <iostream>
#include <QDebug>
#include "Viewer_interface.h"
const QColor Scene_item::defaultColor = QColor(100, 100, 255);

Scene_item::~Scene_item() {
    delete defaultContextMenu;
}

void Scene_item::itemAboutToBeDestroyed(Scene_item* item) {
    if(this == item)
        emit aboutToBeDestroyed();
}


QString modeName(RenderingMode mode) {
    switch(mode)
    {
    case Points:
        return QObject::tr("points");
    case Wireframe:
        return QObject::tr("wire");
    case Flat:
        return QObject::tr("flat");
    case FlatPlusEdges:
        return QObject::tr("flat+edges");
    case Gouraud:
        return QObject::tr("Gouraud");
    case PointsPlusNormals:
        return QObject::tr("pts+normals");
    case Splatting:
        return QObject::tr("splats");
    default:
        Q_ASSERT(false);
        return QObject::tr("unknown");
    }
}

const char* slotName(RenderingMode mode) {
    switch(mode)
    {
    case Points:
        return SLOT(setPointsMode());
    case Wireframe:
        return SLOT(setWireframeMode());
    case Flat:
        return SLOT(setFlatMode());
    case FlatPlusEdges:
        return SLOT(setFlatPlusEdgesMode());
    case Gouraud:
        return SLOT(setGouraudMode());
    case PointsPlusNormals:
        return SLOT(setPointsPlusNormalsMode());
    case Splatting:
        return SLOT(setSplattingMode());
    default:
        Q_ASSERT(false);
        return "";
    }
}

// Rendering mode as a human readable string
QString Scene_item::renderingModeName() const
{
    return modeName(renderingMode());
} 
QMenu* Scene_item::contextMenu()
{
    if(defaultContextMenu) {
        defaultContextMenu->setTitle(name());
        return defaultContextMenu;
    }

    defaultContextMenu = new QMenu(name());
    // defaultContextMenu->addAction(name());
    // defaultContextMenu->addSeparator();
    // QMenu* modeMenu = new QMenu(QObject::tr("Rendering mode"),
    //                             defaultContextMenu);
    for(unsigned int mode = 0; mode < NumberOfRenderingMode;
        ++mode)
    {
        if(!supportsRenderingMode(RenderingMode(mode))) continue;
        QString mName = modeName(RenderingMode(mode));
        QAction* action =
                defaultContextMenu->addAction(tr("Set %1 mode")
                                              .arg(mName),
                                              this,
                                              slotName(RenderingMode(mode)));
        QObject::connect(action, SIGNAL(triggered()),
                         this, SIGNAL(itemChanged()));
    }
    // defaultContextMenu->addAction(modeMenu->menuAction());
    return defaultContextMenu;
}

void Scene_item::changed() {
    // emit itemChanged();
}

void Scene_item::selection_changed(bool) {
    // emit itemChanged();
}


void Scene_item::select(double /*orig_x*/,
                        double /*orig_y*/,
                        double /*orig_z*/,
                        double /*dir_x*/,
                        double /*dir_y*/,
                        double /*dir_z*/)
{
}

void Scene_item::compile_shaders()
{
    qFunc.initializeOpenGLFunctions();
    for(int i=0; i<10; i++)
    {
     if(!buffers[i].create())
         qDebug()<<"ERROR";
     if(!vaos[i].create())
         qDebug()<<"ERROR";
    }
    //fill the vertex shader
    const char vertex_shader_source[] =
    {"attribute highp vec4 vertex;\n"
     "attribute highp vec3 normals;\n"
     "attribute highp vec3 colors;\n"
     "uniform highp mat4 mvp_matrix;\n"
     "uniform highp mat4 mv_matrix; \n"
     "varying highp vec4 fP; \n"
     "varying highp vec3 fN; \n"
     "varying highp vec4 color; \n"
     "void main(void)\n"
     "{\n"
     "   color = vec4(colors, 1.0); \n"
     "   fP = mv_matrix * vertex; \n"
     "   fN = mat3(mv_matrix)* normals; \n"
     "   gl_Position = mvp_matrix * vertex; \n"
     "}"

    };
    //fill the fragment shader
    const char fragment_shader_source[]=
    {
         "varying highp vec4 color; \n"
         "varying highp vec4 fP; \n"
         "varying highp vec3 fN; \n"
         "uniform highp vec4 light_pos;  \n"
         "uniform highp vec4 light_diff; \n"
         "uniform highp vec4 light_spec; \n"
         "uniform highp vec4 light_amb;  \n"
         "uniform highp float spec_power ; \n"
         "uniform int is_two_side; \n"

         "void main(void) { \n"

         "   highp vec3 L = light_pos.xyz - fP.xyz; \n"
         "   highp vec3 V = -fP.xyz; \n"
         "   highp vec3 N; \n"
         "   if(fN == highp vec3(0.0,0.0,0.0)) \n"
         "       N = highp vec3(0.0,0.0,0.0); \n"
         "   else \n"
         "       N = normalize(fN); \n"
         "   L = normalize(L); \n"
         "   V = normalize(V); \n"
         "   highp vec3 R = reflect(-L, N); \n"
            "vec4 diffuse; \n"
         "   if(is_two_side == 1) \n"
         "       diffuse = abs(dot(N,L)) * light_diff * color; \n"
         "   else \n"
         "       diffuse = max(dot(N,L), 0.0) * light_diff * color; \n"
         "   highp vec4 specular = pow(max(dot(R,V), 0.0), spec_power) * light_spec; \n"

         "gl_FragColor = vec4((color*light_amb).xyz + diffuse.xyz + specular.xyz,1); \n"
         "} \n"
         "\n"
    };

    QOpenGLShader *vertex_shader = new QOpenGLShader(QOpenGLShader::Vertex);
    if(!vertex_shader->compileSourceCode(vertex_shader_source))
    {
        std::cerr<<"Compiling vertex source FAILED"<<std::endl;
    }

    QOpenGLShader *fragment_shader= new QOpenGLShader(QOpenGLShader::Fragment);
    if(!fragment_shader->compileSourceCode(fragment_shader_source))
    {
        std::cerr<<"Compiling fragmentsource FAILED"<<std::endl;
    }

    if(!rendering_program_with_light.addShader(vertex_shader))
    {
        std::cerr<<"adding vertex shader FAILED"<<std::endl;
    }
    if(!rendering_program_with_light.addShader(fragment_shader))
    {
        std::cerr<<"adding fragment shader FAILED"<<std::endl;
    }
    if(!rendering_program_with_light.link())
    {
        std::cerr<<"linking Program FAILED"<<std::endl;
    }
    rendering_program_with_light.bind();


    //For the edges

    if(!rendering_program_without_light.addShader(vertex_shader))
    {
        std::cerr<<"adding vertex shader FAILED"<<std::endl;
    }
    if(!rendering_program_without_light.addShader(fragment_shader))
    {
        std::cerr<<"adding fragment shader FAILED"<<std::endl;
    }
    if(!rendering_program_without_light.link())
    {
        std::cerr<<"linking Program FAILED"<<std::endl;
    }
    rendering_program_without_light.bind();
}
void Scene_item::attrib_buffers(Viewer_interface* viewer) const
{

    GLint is_both_sides = 0;
    QMatrix4x4 mvp_mat;
    QMatrix4x4 mv_mat;

    //fills the MVP and MV matrices.

    GLdouble d_mat[16];
    viewer->camera()->getModelViewProjectionMatrix(d_mat);
    //Convert the GLdoubles matrices in GLfloats
    for (int i=0; i<16; ++i){
        mvp_mat.data()[i] = GLfloat(d_mat[i]);
    }

    viewer->camera()->getModelViewMatrix(d_mat);
    for (int i=0; i<16; ++i)
        mv_mat.data()[i] = GLfloat(d_mat[i]);

    qFunc.glGetIntegerv(GL_LIGHT_MODEL_TWO_SIDE, &is_both_sides);


    QVector4D position(0.0f,0.0f,1.0f, 1.0f );
    QVector4D ambient(0.4f, 0.4f, 0.4f, 0.4f);
    // Diffuse
    QVector4D diffuse(1.0f, 1.0f, 1.0f, 1.0f);
    // Specular
    QVector4D specular(0.0f, 0.0f, 0.0f, 1.0f);


    rendering_program_with_light.bind();

    rendering_program_with_light.setUniformValue("mvp_matrix", mvp_mat);
    rendering_program_with_light.setUniformValue("mv_matrix", mv_mat);

    rendering_program_with_light.setUniformValue("light_pos", position);
    rendering_program_with_light.setUniformValue("light_diff",diffuse);
    rendering_program_with_light.setUniformValue("light_spec", specular);
    rendering_program_with_light.setUniformValue("light_amb", ambient);
    rendering_program_with_light.setUniformValue("spec_power", 51.8f);
    rendering_program_with_light.setUniformValue("is_two_side", is_both_sides);

    rendering_program_with_light.release();

    rendering_program_without_light.bind();
    rendering_program_without_light.setUniformValue("mvp_matrix", mvp_mat);
    rendering_program_without_light.setUniformValue("mv_matrix", mv_mat);

    rendering_program_without_light.setUniformValue("light_pos", position);
    rendering_program_without_light.setUniformValue("light_diff", diffuse);
    rendering_program_without_light.setUniformValue("light_spec", specular);
    rendering_program_without_light.setUniformValue("light_amb", ambient);
    rendering_program_without_light.setUniformValue("spec_power", 51.8f);
    rendering_program_without_light.setUniformValue("is_two_side", is_both_sides);
    rendering_program_without_light.setAttributeValue("normals", 0.0,0.0,0.0);

    rendering_program_without_light.release();
}
#include "Scene_item.moc"

