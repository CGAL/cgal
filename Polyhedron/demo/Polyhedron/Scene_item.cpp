#include "Scene_item.h"
#include "Scene_interface.h"
#include <QMenu>
#include <iostream>
#include <QDebug>
#include "Viewer_interface.h"
const QColor Scene_item::defaultColor = QColor(100, 100, 255);

Scene_item::~Scene_item() {
    delete defaultContextMenu;
    for(int i=0; i<10; i++)
    {
        vaos[i].destroy();
        buffers[i].destroy();
    }
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
//set-up the shader programs
void Scene_item::compile_shaders()
{
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


    if(!rendering_program_with_texture.addShaderFromSourceFile(QOpenGLShader::Vertex,":/cgal/Polyhedron_3/resources/shader_with_texture.v"))
    {
        std::cerr<<"adding vertex shader FAILED"<<std::endl;
    }
    if(!rendering_program_with_texture.addShaderFromSourceFile(QOpenGLShader::Fragment,":/cgal/Polyhedron_3/resources/shader_with_texture.f"))
    {
        std::cerr<<"adding fragment shader FAILED"<<std::endl;
    }
    if(!rendering_program_with_texture.link())
    {
        std::cerr<<"linking Program FAILED"<<std::endl;
    }
    rendering_program_with_texture.bind();

}

// set-up the uniform attributes of the shader programs.
void Scene_item::attrib_buffers(Viewer_interface* viewer, int program_name) const
{
    GLint is_both_sides = 0;
    QMatrix4x4 mvp_mat;
    QMatrix4x4 mv_mat;
    QMatrix4x4 f_mat;
    f_mat.setToIdentity();
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
    QColor temp = this->color();
    switch(program_name)
    {
    case PROGRAM_WITH_LIGHT:
        shader_programs[PROGRAM_WITH_LIGHT]->bind();
        shader_programs[PROGRAM_WITH_LIGHT]->setUniformValue("mvp_matrix", mvp_mat);
        shader_programs[PROGRAM_WITH_LIGHT]->setUniformValue("mv_matrix", mv_mat);
        shader_programs[PROGRAM_WITH_LIGHT]->setUniformValue("light_pos", position);
        shader_programs[PROGRAM_WITH_LIGHT]->setUniformValue("light_diff",diffuse);
        shader_programs[PROGRAM_WITH_LIGHT]->setUniformValue("light_spec", specular);
        shader_programs[PROGRAM_WITH_LIGHT]->setUniformValue("light_amb", ambient);
        shader_programs[PROGRAM_WITH_LIGHT]->setUniformValue("spec_power", 51.8f);
        shader_programs[PROGRAM_WITH_LIGHT]->setUniformValue("is_two_side", is_both_sides);
        shader_programs[PROGRAM_WITH_LIGHT]->release();
        break;
    case PROGRAM_WITHOUT_LIGHT:
        shader_programs[PROGRAM_WITHOUT_LIGHT]->bind();
        shader_programs[PROGRAM_WITHOUT_LIGHT]->setUniformValue("mvp_matrix", mvp_mat);
        shader_programs[PROGRAM_WITHOUT_LIGHT]->setUniformValue("mv_matrix", mv_mat);

        shader_programs[PROGRAM_WITHOUT_LIGHT]->setUniformValue("light_pos", position);
        shader_programs[PROGRAM_WITHOUT_LIGHT]->setUniformValue("light_diff", diffuse);
        shader_programs[PROGRAM_WITHOUT_LIGHT]->setUniformValue("light_spec", specular);
        shader_programs[PROGRAM_WITHOUT_LIGHT]->setUniformValue("light_amb", ambient);
        shader_programs[PROGRAM_WITHOUT_LIGHT]->setUniformValue("spec_power", 51.8f);
        shader_programs[PROGRAM_WITHOUT_LIGHT]->setUniformValue("is_two_side", is_both_sides);
        shader_programs[PROGRAM_WITHOUT_LIGHT]->setAttributeValue("normals", 0.0,0.0,0.0);
        shader_programs[PROGRAM_WITHOUT_LIGHT]->setUniformValue("f_matrix",f_mat);


        shader_programs[PROGRAM_WITHOUT_LIGHT]->release();
        break;
    case PROGRAM_WITH_TEXTURE:
        if(is_selected)
        {

            shader_programs[PROGRAM_WITH_TEXTURE]->setAttributeValue("color_facets", temp.lighter(120).redF(),temp.lighter(120).greenF(), temp.lighter(120).blueF());
        }
        else
        {
            shader_programs[PROGRAM_WITH_TEXTURE]->setAttributeValue("color_facets", temp.redF(),temp.greenF(), temp.blueF());
        }

        shader_programs[PROGRAM_WITH_TEXTURE]->bind();
        shader_programs[PROGRAM_WITH_TEXTURE]->setUniformValue("mvp_matrix", mvp_mat);
        shader_programs[PROGRAM_WITH_TEXTURE]->setUniformValue("mv_matrix", mv_mat);
        shader_programs[PROGRAM_WITH_TEXTURE]->setUniformValue("light_pos", position);
        shader_programs[PROGRAM_WITH_TEXTURE]->setUniformValue("light_diff",diffuse);
        shader_programs[PROGRAM_WITH_TEXTURE]->setUniformValue("light_spec", specular);
        shader_programs[PROGRAM_WITH_TEXTURE]->setUniformValue("light_amb", ambient);
        shader_programs[PROGRAM_WITH_TEXTURE]->setUniformValue("spec_power", 51.8f);
        shader_programs[PROGRAM_WITH_TEXTURE]->setUniformValue("s_texture",0);
        shader_programs[PROGRAM_WITH_TEXTURE]->setUniformValue("f_matrix",f_mat);


        shader_programs[PROGRAM_WITH_TEXTURE]->release();
        break;
    case PROGRAM_WITH_TEXTURED_EDGES:
        shader_programs[PROGRAM_WITH_TEXTURED_EDGES]->bind();
        if(is_selected)
        {
            shader_programs[PROGRAM_WITH_TEXTURED_EDGES]->setUniformValue("color_lines",QVector3D(0.0,0.0,0.0));
        }
        else
        {
            shader_programs[PROGRAM_WITH_TEXTURED_EDGES]->setUniformValue("color_lines", QVector3D(temp.lighter(50).redF(), temp.lighter(50).greenF(), temp.lighter(50).blueF()));

        }

        shader_programs[PROGRAM_WITH_TEXTURED_EDGES]->setUniformValue("mvp_matrix", mvp_mat);
        shader_programs[PROGRAM_WITH_TEXTURED_EDGES]->setUniformValue("s_texture",0);
        shader_programs[PROGRAM_WITH_TEXTURED_EDGES]->release();
        break;
    case PROGRAM_INSTANCED:

        shader_programs[PROGRAM_INSTANCED]->bind();
        shader_programs[PROGRAM_INSTANCED]->setUniformValue("mvp_matrix", mvp_mat);
        shader_programs[PROGRAM_INSTANCED]->setUniformValue("mv_matrix", mv_mat);

        shader_programs[PROGRAM_INSTANCED]->setUniformValue("light_pos", position);
        shader_programs[PROGRAM_INSTANCED]->setUniformValue("light_diff",diffuse);
        shader_programs[PROGRAM_INSTANCED]->setUniformValue("light_spec", specular);
        shader_programs[PROGRAM_INSTANCED]->setUniformValue("light_amb", ambient);
        shader_programs[PROGRAM_INSTANCED]->setUniformValue("spec_power", 51.8f);
        shader_programs[PROGRAM_INSTANCED]->setUniformValue("is_two_side", is_both_sides);
        shader_programs[PROGRAM_INSTANCED]->release();

        break;
    case PROGRAM_INSTANCED_WIRE:
        shader_programs[PROGRAM_INSTANCED_WIRE]->bind();
        shader_programs[PROGRAM_INSTANCED_WIRE]->setUniformValue("mvp_matrix", mvp_mat);
        shader_programs[PROGRAM_INSTANCED_WIRE]->release();
        break;
    }
}


QOpenGLShaderProgram* Scene_item::getShaderProgram(int name, Viewer_interface * viewer) const
{
    switch(name)
    {
    case PROGRAM_WITH_LIGHT:
        if(shader_programs[PROGRAM_WITH_LIGHT])
        {
            return shader_programs[PROGRAM_WITH_LIGHT];
        }

        else
        {

            QOpenGLShaderProgram *program = new QOpenGLShaderProgram(viewer);
            if(!program->addShaderFromSourceFile(QOpenGLShader::Vertex,":/cgal/Polyhedron_3/resources/shader_with_light.v"))
            {
                std::cerr<<"adding vertex shader FAILED"<<std::endl;
            }
            if(!program->addShaderFromSourceFile(QOpenGLShader::Fragment,":/cgal/Polyhedron_3/resources/shader_with_light.f"))
            {
                std::cerr<<"adding fragment shader FAILED"<<std::endl;
            }
            program->link();
            shader_programs[PROGRAM_WITH_LIGHT] = program;
            return program;
        }
        break;
    case PROGRAM_WITHOUT_LIGHT:
        if( shader_programs[PROGRAM_WITHOUT_LIGHT])
        {
            return shader_programs[PROGRAM_WITHOUT_LIGHT];
        }
        else
        {
            QOpenGLShaderProgram *program = new QOpenGLShaderProgram(viewer);
            if(!program->addShaderFromSourceFile(QOpenGLShader::Vertex,":/cgal/Polyhedron_3/resources/shader_without_light.v"))
            {
                std::cerr<<"adding vertex shader FAILED"<<std::endl;
            }
            if(!program->addShaderFromSourceFile(QOpenGLShader::Fragment,":/cgal/Polyhedron_3/resources/shader_without_light.f"))
            {
                std::cerr<<"adding fragment shader FAILED"<<std::endl;
            }
            program->link();
            shader_programs[PROGRAM_WITHOUT_LIGHT] = program;
            return program;
        }
        break;
    case PROGRAM_WITH_TEXTURE:
        if( shader_programs[PROGRAM_WITH_TEXTURE])
        {
            return shader_programs[PROGRAM_WITH_TEXTURE];
        }
        else
        {
            QOpenGLShaderProgram *program = new QOpenGLShaderProgram(viewer);
            if(!program->addShaderFromSourceFile(QOpenGLShader::Vertex,":/cgal/Polyhedron_3/resources/shader_with_texture.v"))
            {
                std::cerr<<"adding vertex shader FAILED"<<std::endl;
            }
            if(!program->addShaderFromSourceFile(QOpenGLShader::Fragment,":/cgal/Polyhedron_3/resources/shader_with_texture.f"))
            {
                std::cerr<<"adding fragment shader FAILED"<<std::endl;
            }
            program->link();
            shader_programs[PROGRAM_WITH_TEXTURE] = program;
            return program;
        }
        break;
    case PROGRAM_WITH_TEXTURED_EDGES:
        if( shader_programs[PROGRAM_WITH_TEXTURED_EDGES])
        {
            return shader_programs[PROGRAM_WITH_TEXTURED_EDGES];
        }
        else
        {
            QOpenGLShaderProgram *program = new QOpenGLShaderProgram(viewer);
            if(!program->addShaderFromSourceFile(QOpenGLShader::Vertex,":/cgal/Polyhedron_3/resources/shader_with_textured_edges.v" ))
            {
                std::cerr<<"adding vertex shader FAILED"<<std::endl;
            }
            if(!program->addShaderFromSourceFile(QOpenGLShader::Fragment,":/cgal/Polyhedron_3/resources/shader_with_textured_edges.f" ))
            {
                std::cerr<<"adding fragment shader FAILED"<<std::endl;
            }
            program->link();
            shader_programs[PROGRAM_WITH_TEXTURED_EDGES] = program;
            return program;

        }
        break;
    case PROGRAM_INSTANCED:
        if( shader_programs[PROGRAM_INSTANCED])
        {
            return shader_programs[PROGRAM_INSTANCED];
        }
        else
        {
            QOpenGLShaderProgram *program = new QOpenGLShaderProgram(viewer);
            if(!program->addShaderFromSourceFile(QOpenGLShader::Vertex,":/cgal/Polyhedron_3/resources/shader_instanced.v" ))
            {
                std::cerr<<"adding vertex shader FAILED"<<std::endl;
            }
            if(!program->addShaderFromSourceFile(QOpenGLShader::Fragment,":/cgal/Polyhedron_3/resources/shader_with_light.f" ))
            {
                std::cerr<<"adding fragment shader FAILED"<<std::endl;
            }
            program->link();
            shader_programs[PROGRAM_INSTANCED] = program;
            return program;

        }
        break;
    case PROGRAM_INSTANCED_WIRE:
        if( shader_programs[PROGRAM_INSTANCED_WIRE])
        {
            return shader_programs[PROGRAM_INSTANCED_WIRE];
        }
        else
        {
            QOpenGLShaderProgram *program = new QOpenGLShaderProgram(viewer);
            if(!program->addShaderFromSourceFile(QOpenGLShader::Vertex,":/cgal/Polyhedron_3/resources/shader_instanced.v" ))
            {
                std::cerr<<"adding vertex shader FAILED"<<std::endl;
            }
            if(!program->addShaderFromSourceFile(QOpenGLShader::Fragment,":/cgal/Polyhedron_3/resources/shader_without_light.f" ))
            {
                std::cerr<<"adding fragment shader FAILED"<<std::endl;
            }
            program->link();
            shader_programs[PROGRAM_INSTANCED_WIRE] = program;
            return program;

        }
        break;
    default:
        std::cerr<<"ERROR : Program not found."<<std::endl;
    }
}
#include "Scene_item.moc"

