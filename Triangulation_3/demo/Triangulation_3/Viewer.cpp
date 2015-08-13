#include <boost/config.hpp>
#include <QDebug>

#if defined(BOOST_MSVC)
// Avoid warning concerning spatial_sort(QList::begin(), QList.end() QT "bug"
#  pragma warning(disable: 4267 )
#  pragma warning(disable: 4244 )
#endif

#include "Viewer.h"




using namespace std;


void Viewer::init()
{
    initializeOpenGLFunctions();
    glDrawArraysInstanced = (PFNGLDRAWARRAYSINSTANCEDARBPROC)this->context()->getProcAddress("glDrawArraysInstancedARB");
    if(!glDrawArraysInstanced)
    {
        qDebug()<<"glDrawArraysInstancedARB : extension not found. Spheres will be displayed as points.";
        extension_is_found = false;
    }
    else
        extension_is_found = true;

    glVertexAttribDivisor = (PFNGLVERTEXATTRIBDIVISORARBPROC)this->context()->getProcAddress("glVertexAttribDivisorARB");
    if(!glDrawArraysInstanced)
    {
        qDebug()<<"glVertexAttribDivisorARB : extension not found. Spheres will be displayed as points.";
        extension_is_found = false;
    }
    else
        extension_is_found = true;
//    extension_is_found = false;

    /* Initial timer for playing incremental construction */
    m_pTimer = new QTimer(this);
    connect(m_pTimer, SIGNAL(timeout()), this, SLOT(incremental_insert()));

    /* Scene inits */
    setBackgroundColor(::Qt::white);
    // scene are defined by a sphere of 2.0, camera at the center, i.e. (0, 0, 0)
    setSceneCenter( qglviewer::Vec(-0.,-0.,-0.) );
    setSceneRadius( 2. );
    // show text message
    setTextIsEnabled(true);
    setForegroundColor(::Qt::red);
    setFont(QFont("Arial Black", 16, QFont::Bold));

    /* OpenGL inits */
    // Increase the material shininess, so that the difference between
    // the two versions of the spiral is more visible.
    glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 50.0);
    GLfloat specular_color[4] = { 0.8f, 0.8f, 0.8f, 1.0 };
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR,  specular_color);
    // Set Smooth Shading
    ::glShadeModel(GL_SMOOTH);

    // depth buffer setup
    ::glClearDepth(1.0f);
    ::glEnable(GL_DEPTH_TEST);
    ::glDepthFunc(GL_LEQUAL);
    ::glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);

    // enable semi-transparent culling planes
    ::glEnable(GL_BLEND);
    ::glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    // anti-aliasing, i.e. reduce jaggedness (if the OpenGL driver permits that)
    ::glEnable(GL_POINT_SMOOTH);
    ::glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
    ::glEnable(GL_LINE_SMOOTH);
    ::glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);

    /* Add mouse and key description */
    setKeyDescription( Qt::CTRL + Qt::Key_G, tr("Generate points") );
    setKeyDescription( Qt::CTRL + Qt::Key_O, tr("Load points") );
    setKeyDescription( Qt::CTRL + Qt::Key_S, tr("Save points") );
    setKeyDescription( Qt::CTRL + Qt::Key_Comma, tr("Preference") );
    setKeyDescription( Qt::CTRL + Qt::Key_H, tr("Hide Kernel Demo") );
    setKeyDescription( Qt::CTRL + Qt::Key_Q, tr("Quit Kernel Demo") );
    setKeyDescription( Qt::Key_Return,
                       tr("Insert new point to triangulation in <u>Input-Point</u> mode") );
    setKeyDescription( Qt::Key_Escape,
                       tr("Cancel insertion in <u>Input-Point</u> mode;<br>")
                       + tr("Cancel current selection in <u>Select</u> mode") );
    setKeyDescription( Qt::Key_Delete, tr("Delete selected vertices in <u>Select</u> mode") );
#if QGLVIEWER_VERSION >= 0x020500
    setMouseBindingDescription(Qt::NoModifier, Qt::LeftButton,
                               tr("Hold to move new point in <u>Input-Point</u> mode;<br>")
                               + tr("Hold to move a vertex in <u>Move</u> mode") );
    setMouseBindingDescription(Qt::ShiftModifier, Qt::LeftButton,
                               tr("Click to insert a vertex in <u>Input-Vertex</u> mode;<br>")
                               + tr("Click to insert a point in <u>Input-Point</u> mode;<br>")
                               + tr("Click or Drag to select multiple points in <u>Select</u> mode;<br>")
                               + tr("Click to place a query point in <u>Find-Nearest-Neighbor</u> mode;<br>")
                               + tr("Click to place a query point in <u>Show-Empty-Sphere</u> mode") );
    setMouseBindingDescription(Qt::ControlModifier, Qt::LeftButton,
                               tr("Drag to add vertices to current selection in <u>Select</u> mode") );
#else
    setMouseBindingDescription( Qt::LeftButton,
                                tr("Hold to move new point in <u>Input-Point</u> mode;<br>")
                                + tr("Hold to move a vertex in <u>Move</u> mode") );
    setMouseBindingDescription( Qt::SHIFT + Qt::LeftButton,
                                tr("Click to insert a vertex in <u>Input-Vertex</u> mode;<br>")
                                + tr("Click to insert a point in <u>Input-Point</u> mode;<br>")
                                + tr("Click or Drag to select multiple points in <u>Select</u> mode;<br>")
                                + tr("Click to place a query point in <u>Find-Nearest-Neighbor</u> mode;<br>")
                                + tr("Click to place a query point in <u>Show-Empty-Sphere</u> mode") );
    setMouseBindingDescription( Qt::CTRL + Qt::LeftButton,
                                tr("Drag to add vertices to current selection in <u>Select</u> mode") );
#endif
    compile_shaders();
    are_buffers_initialized = false;
}

QString Viewer::helpString() const
{
    QString text("<h1>3D Triangulation Demo</h1>");

    text += "This example illustrates a generic interactive demo for 3D Triangulation in CGAL. ";
    text += "This demo could be used as a simple skeleton ";
    text += "for potential demos of other 3D packages or for teaching CGAL.<br><br>";

    text += "The key feature is to edit vertices/points with mouse.";
    text += "There are several modes:<br><br>";

    text += " - <u>Normal Mode</u>: ";
    text += "Rotate, zoom, or translate camera using mouse.<br>";
    text += " - <u>Insert Vertex</u>: ";
    text += "Insert a vertex on the surface of the trackball ";
    text += "and the triangulation will be updated correspondingly.<br>";
    text += " - <u>Insert Point</u>: ";
    text += "Insert a point on the surface of the trackball. ";
    text += "Its conflict region will be highlighted. ";
    text += "When the new point is moving, ";
    text += "its conflict region will be updated correspondingly.<br>";
    text += " - <u>Select</u>: ";
    text += "Click or drag mouse left button to select multiple points.<br>";
    text += " - <u>Move</u>: Hold mouse left button to move a vertex ";
    text += "and the triangulation will be updated correspondingly.<br>";
    text += " - <u>Find Nearest Neighbor</u>: ";
    text += "Place a query point and its nearest neighbor will be highlighted.<br>";
    text += " - <u>Show Empty Sphere</u>: ";
    text += "Place a query point, locate the point in a cell ";
    text += "and then show the empty sphere of that cell. ";
    text += "An empty sphere of a cell is a sphere ";
    text += "with all four vertices of the cell lying on it ";
    text += "and no other vertices inside it.<br><br>";
    text += "<b>Shift+Wheel</b> to resize the trackball when it exists. ";
    text += "See <b>Mouse</b> page for more details.<br><br>";

    text += "Other basic features include:<br>";
    text += " - Randomly generate points,<br>";
    text += " - Read/Write files,<br>";
    text += " - Show vertices, Voronoi edges, Delaunay edges, and/or facets,<br>";
    text += " - Incremental Construct: ";
    text += "Re-construct the current triangulation incrementally. ";
    text += "If no triangulation exists yet, randomly generate 100 points ";
    text += "and construct a Delaunay triangulation of those points.<br>";

    return text;
}
/*************************************************************/
/*  Shaders functions */

void Viewer::compile_shaders()
{

    for(int i=0; i< vboSize; i++)
        buffers[i].create();
    for(int i=0; i< vaoSize; i++)
        vao[i].create();

    draw_cylinder(m_fSizeDEdge,25,points_cylinder,normals_cylinder);
    draw_sphere(m_fSizeVertex, 15, points_sphere, normals_sphere);
    //Vertex source code
    const char vertex_source[] =
    {
        "#version 120 \n"
        "attribute highp vec4 vertex;\n"
        "uniform highp mat4 mvp_matrix;\n"
        "void main(void)\n"
        "{\n"
        "   gl_Position = mvp_matrix * vertex; \n"
        "}"
    };
    //Fragment source code
    const char fragment_source[] =
    {
        "#version 120 \n"
        "uniform highp vec4 color; \n"
        "void main(void) { \n"
        "gl_FragColor = color; \n"
        "} \n"
        "\n"
    };
    QOpenGLShader *vertex_shader = new QOpenGLShader(QOpenGLShader::Vertex);
    if(!vertex_shader->compileSourceCode(vertex_source))
    {
        std::cerr<<"Compiling vertex source FAILED"<<std::endl;
    }

    QOpenGLShader *fragment_shader= new QOpenGLShader(QOpenGLShader::Fragment);
    if(!fragment_shader->compileSourceCode(fragment_source))
    {
        std::cerr<<"Compiling fragmentsource FAILED"<<std::endl;
    }

    if(!rendering_program.addShader(vertex_shader))
    {
        std::cerr<<"adding vertex shader FAILED"<<std::endl;
    }
    if(!rendering_program.addShader(fragment_shader))
    {
        std::cerr<<"adding fragment shader FAILED"<<std::endl;
    }
    if(!rendering_program.link())
    {
        std::cerr<<"linking Program FAILED"<<std::endl;
    }
    rendering_program.bind();

    // sphere program
    //Vertex source code
    const char vertex_source_spheres[] =
    {
        "#version 120 \n"
        "attribute highp vec4 vertex;\n"
        "attribute highp vec3 normal;\n"
        "attribute highp vec4 center;\n"

        "uniform highp mat4 mvp_matrix;\n"
        "uniform highp mat4 mv_matrix; \n"
        "varying highp vec4 fP; \n"
        "varying highp vec3 fN; \n"
        "void main(void)\n"
        "{\n"
        "   fP = mv_matrix * vertex; \n"
        "   fN = mat3(mv_matrix)* normal; \n"
        "   gl_Position = mvp_matrix * (vertex+vec4(center.xyz, 0.0)); \n"
        "}"
    };
    //Fragment source code
    const char fragment_source_spheres[] =
    {
        "#version 120 \n"
        "varying highp vec4 fP; \n"
        "varying highp vec3 fN; \n"
        "uniform vec4 color; \n"
        "uniform highp vec4 light_pos;  \n"
        "uniform highp vec4 light_diff; \n"
        "uniform highp vec4 light_spec; \n"
        "uniform highp vec4 light_amb;  \n"
        "uniform float spec_power ; \n"

        "void main(void) { \n"

        "   vec3 L = light_pos.xyz - fP.xyz; \n"
        "   vec3 V = -fP.xyz; \n"

        "   vec3 N = normalize(fN); \n"
        "   L = normalize(L); \n"
        "   V = normalize(V); \n"

        "   vec3 R = reflect(-L, N); \n"
        "   vec4 diffuse = abs(dot(N,L)) * light_diff*color; \n"
        "   vec4 specular = pow(max(dot(R,V), 0.0), spec_power) * light_spec; \n"

        "gl_FragColor = color*light_amb + diffuse + specular; \n"
        "} \n"
        "\n"
    };
    QOpenGLShader *vertex_shader_spheres = new QOpenGLShader(QOpenGLShader::Vertex);
    if(!vertex_shader_spheres->compileSourceCode(vertex_source_spheres))
    {
        std::cerr<<"Compiling vertex source FAILED"<<std::endl;
    }


    QOpenGLShader *fragment_shader_spheres= new QOpenGLShader(QOpenGLShader::Fragment);
    if(!fragment_shader_spheres->compileSourceCode(fragment_source_spheres))
    {
        std::cerr<<"Compiling fragmentsource FAILED"<<std::endl;
    }


    if(!rendering_program_spheres.addShader(vertex_shader_spheres))
    {
        std::cerr<<"adding vertex shader FAILED"<<std::endl;
    }
    if(!rendering_program_spheres.addShader(fragment_shader_spheres))
    {
        std::cerr<<"adding fragment shader FAILED"<<std::endl;
    }
    if(!rendering_program_spheres.link())
    {
        std::cerr<<"linking Program FAILED"<<std::endl;
    }
    rendering_program_spheres.bind();

    // cylinder program
    //Vertex source code
    const char vertex_source_cylinders[] =
    {
        "#version 120 \n"
        "attribute highp vec4 vertex;\n"
        "attribute highp vec3 normal;\n"
        "attribute highp vec4 transfo1;\n"
        "attribute highp vec4 transfo2;\n"
        "attribute highp vec4 transfo3;\n"
        "attribute highp vec4 transfo4;\n"
        "mat4 transfo = mat4(transfo1, transfo2, transfo3, transfo4); \n"
        "uniform highp mat4 mvp_matrix;\n"
        "uniform highp mat4 mv_matrix; \n"
        "varying highp vec4 fP; \n"
        "varying highp vec3 fN; \n"
        "void main(void)\n"
        "{\n"
        "   fP = mv_matrix * vertex; \n"
        "   vec4 TN = transfo*vec4(normal,1.0); \n"
        "   fN = mat3(mv_matrix)* TN.xyz; \n"
        "   gl_Position =  mvp_matrix * transfo* vertex; \n"
        "}"
    };



    QOpenGLShader *vertex_shader_cylinders = new QOpenGLShader(QOpenGLShader::Vertex);
    if(!vertex_shader_cylinders->compileSourceCode(vertex_source_cylinders))
    {
        std::cerr<<"Compiling vertex source FAILED"<<std::endl;
    }
    if(!rendering_program_cylinders.addShader(vertex_shader_cylinders))
    {
        std::cerr<<"adding vertex shader FAILED"<<std::endl;
    }
    if(!rendering_program_cylinders.addShader(fragment_shader_spheres))
    {
        std::cerr<<"adding fragment shader FAILED"<<std::endl;
    }
    if(!rendering_program_cylinders.link())
    {
        std::cerr<<"linking Program FAILED"<<std::endl;
    }
    rendering_program_cylinders.bind();

}

void Viewer::compute_elements()
{
    pos_points->resize(0);
    pos_delaunay->resize(0);
    pos_voronoi->resize(0);
    pos_facets->resize(0);
    pos_newPoint->resize(0);
    pos_newFacet->resize(0);
    pos_selectedVertex->resize(0);
    pos_movingPoint->resize(0);
    pos_queryPoint->resize(0);
    pos_nearest_neighbor->resize(0);
    pos_emptyFacet->resize(0);
    transfo1_delaunay->resize(0);
    transfo2_delaunay->resize(0);
    transfo3_delaunay->resize(0);
    transfo4_delaunay->resize(0);

    pos_trackBall->resize(0);
    draw_sphere(m_fRadius, 35,points_trackBall, normals_trackBall);
    for(int i=0; i<3; i++)
        pos_trackBall->push_back(0.0);
    pos_emptySphere->resize(0);
    // Draw vertices
    if ( m_pScene->m_dt.number_of_vertices()>0 ) {

        for(QList<Vertex_handle>::iterator vit = m_pScene->m_vhArray.begin();
            vit < m_pScene->m_vhArray.end(); ++vit) {
            if( m_curMode == SELECT && (*vit)->isSeled() )  continue;
            if( (*vit) == m_nearestNb ) continue;
            drawVertex( (*vit)->point(), pos_points );
        }//end-for-points
    }//end-if-points
    //Draw Delaunay Edges
    for(edges_iterator eit = m_pScene->m_dt.finite_edges_begin();
        eit != m_pScene->m_dt.finite_edges_end(); ++eit) {
        Segment_3 seg = m_pScene->m_dt.segment(*eit);
        drawEdge( seg.vertex(0), seg.vertex(1), pos_delaunay );

        Vector_3 v = seg.vertex(1) - seg.vertex(0);
        float length = sqrt( v.squared_length() );
        // normalize
        v = v / length;
        // compute the angle: cos theta = v.z/1.0
        GLfloat angle = acos( v.y() ) / M_PI * 180;
        QMatrix4x4 matrix;
        matrix.setToIdentity();
        // move to "from" point
        matrix.translate( seg.vertex(0).x(), seg.vertex(0).y(), seg.vertex(0).z() );
        // rotate from z-axis to from-->to
        //  axis: cross product of z-axis and from-->to
        matrix.rotate( angle, v.z(), 0.0f,-v.x());
        matrix.scale(1,length,1);
        // stock the transformation
        for(int i=0; i<4; i++)
            transfo1_delaunay->push_back((float)matrix.data()[i]);
        for(int i=4; i<8; i++)
            transfo2_delaunay->push_back((float)matrix.data()[i]);
        for(int i=8; i<12; i++)
            transfo3_delaunay->push_back((float)matrix.data()[i]);
        for(int i=12; i<16; i++)
            transfo4_delaunay->push_back((float)matrix.data()[i]);
    }//end-for-edges
    // Draw Voronoi edges
    for(facets_iterator fit = m_pScene->m_dt.finite_facets_begin();
        fit != m_pScene->m_dt.finite_facets_end(); ++fit) {
        Object_3 o = m_pScene->m_dt.dual(*fit);
        if (const Segment_3 *s = CGAL::object_cast<Segment_3>(&o)) {
            drawEdge( s->vertex(0), s->vertex(1), pos_voronoi );

            Vector_3 v = s->vertex(1) - s->vertex(0);
            float length = sqrt( v.squared_length() );
            // normalize
            v = v / length;
            // compute the angle: cos theta = v.z/1.0
            GLfloat angle = acos( v.y() ) / M_PI * 180;
            QMatrix4x4 matrix;
            matrix.setToIdentity();
            // move to "from" point
            matrix.translate( s->vertex(0).x(), s->vertex(0).y(), s->vertex(0).z() );
            // rotate from z-axis to from-->to
            //  axis: cross product of z-axis and from-->to
            matrix.rotate( angle, v.z(), 0.0f,-v.x());
            matrix.scale(1,length,1);
            // stock the transformation
            for(int i=0; i<4; i++)
                transfo1_voronoi->push_back((float)matrix.data()[i]);
            for(int i=4; i<8; i++)
                transfo2_voronoi->push_back((float)matrix.data()[i]);
            for(int i=8; i<12; i++)
                transfo3_voronoi->push_back((float)matrix.data()[i]);
            for(int i=12; i<16; i++)
                transfo4_voronoi->push_back((float)matrix.data()[i]);
        } else if (const Ray_3 *r = CGAL::object_cast<Ray_3>(&o)) {
            drawEdge( r->point(0),  // the source of the ray
                      r->point(1),  // another point on the ray, different from the source
                      pos_voronoi );
            Vector_3 v = r->point(1) - r->point(0);
            float length = sqrt( v.squared_length() );
            // normalize
            v = v / length;
            // compute the angle: cos theta = v.z/1.0
            GLfloat angle = acos( v.y() ) / M_PI * 180;
            QMatrix4x4 matrix;
            matrix.setToIdentity();
            // move to "from" point
            matrix.translate(r->point(0).x(), r->point(0).y(), r->point(0).z() );
            // rotate from z-axis to from-->to
            //  axis: cross product of z-axis and from-->to
            matrix.rotate( angle, v.z(), 0.0f,-v.x());
            matrix.scale(1,length,1);
            // stock the transformation
            for(int i=0; i<4; i++)
                transfo1_voronoi->push_back((float)matrix.data()[i]);
            for(int i=4; i<8; i++)
                transfo2_voronoi->push_back((float)matrix.data()[i]);
            for(int i=8; i<12; i++)
                transfo3_voronoi->push_back((float)matrix.data()[i]);
            for(int i=12; i<16; i++)
                transfo4_voronoi->push_back((float)matrix.data()[i]);
        }
    }//end-for-edges
    // Draw facets
    for(facets_iterator fit = m_pScene->m_dt.finite_facets_begin();
        fit != m_pScene->m_dt.finite_facets_end(); ++fit) {
        drawFacet( m_pScene->m_dt.triangle(*fit), pos_facets);
    }//end-for-facets
    // Draw the newly inserted point
    drawVertex( m_newPt, pos_newPoint);
    // Draw conflict region
    for(QList<Facet>::iterator fit = m_boundaryFacets.begin();
        fit < m_boundaryFacets.end(); ++fit) {
        if( m_pScene->m_dt.is_infinite(*fit) )  continue;
        drawFacet( m_pScene->m_dt.triangle(*fit), pos_newFacet); //semi-transparent purple
    }//end-for-facets
    // Highlight the selected vertices
    for(QList<int>::iterator vit=m_vidSeled.begin(); vit<m_vidSeled.end(); ++vit) {
        drawVertex( m_pScene->m_vhArray.at(*vit)->point(), pos_selectedVertex );
    }//end-for-seledpts
    if( m_isMoving ) {
        // Highlight the moving point
        drawVertex( m_pScene->m_vhArray.at( m_vidMoving )->point(), pos_movingPoint );
    }//end-if-v
    // Draw the nearest neighbor
    if( m_nearestNb != NULL ) {
        drawVertex( m_queryPt, pos_queryPoint);
        drawVertex( m_nearestNb->point(), pos_nearest_neighbor);
    }
    if( m_hasEmptyS ) {
        // Draw the query point
        drawVertex( m_queryPt, pos_queryPoint );
        // Draw the cell containing that point
        for(int i=0; i<4; ++i) {
            if( m_pScene->m_dt.is_infinite(m_cellContain, i) )  continue;
            drawFacet( m_pScene->m_dt.triangle( m_cellContain, i ), pos_emptyFacet );
        }//end-for-facets
        // Draw the sphere
        pos_emptySphere->push_back(m_centerPt.x());
        pos_emptySphere->push_back(m_centerPt.y());
        pos_emptySphere->push_back(m_centerPt.z());
        draw_sphere(m_fREmptyS, 35, points_emptySphere, normals_emptySphere);
    }



    // Draw all points during incremental mode
    if( !m_incrementalPts.isEmpty() ) {
        // draw the rest to-be-inserted vertices
        for(QList<Point_3>::iterator pit=m_incrementalPts.begin();
            pit < m_incrementalPts.end(); ++pit) {
            drawVertex( (*pit), incremental_points);
        }

        switch( m_curStep ) {
        case NEWPT:

            // Highlight the next-to-insert point
            drawVertex( m_curIncPt, incremental_next_point );
            break;
        case CELL:  // show the tetrahedron that contains the point
            // Highlight the next-to-insert vertex
            drawVertex( m_curIncPt, incremental_next_point );
            // Draw the cell containing that point
            for(int i=0; i<4; ++i) {
                if( m_pScene->m_dt.is_infinite(m_cellContain, i) )  continue;
                drawFacet( m_pScene->m_dt.triangle( m_cellContain, i ), incremental_facet );
            }//end-for-facets
            break;
        case CONFLICT:  // show the conflict region
            // Highlight the next-to-insert vertex
            drawVertex( m_curIncPt, incremental_next_point );
            // Draw conflict region
            for(QList<Facet>::iterator fit = m_boundaryFacets.begin();
                fit < m_boundaryFacets.end(); ++fit) {
                if( m_pScene->m_dt.is_infinite(*fit) )  continue;
                drawFacet( m_pScene->m_dt.triangle(*fit), incremental_conflict ); //semi-transparent purple
            }//end-for-facets
            break;
        default:
            break;
        }//end-of=switch
    }//end-if-incpts



}

void Viewer::initialize_buffers()
{
    rendering_program.bind();
    {
        //Points
        vao[0].bind();
        buffers[0].bind();
        buffers[0].allocate(pos_points->data(), pos_points->size()*sizeof(float));
        poly_vertexLocation[0] = rendering_program.attributeLocation("vertex");
        rendering_program.enableAttributeArray(poly_vertexLocation[0]);
        rendering_program.setAttributeBuffer(poly_vertexLocation[0],GL_FLOAT,0,3);
        buffers[0].release();
        vao[0].release();
        //Delaunay Edges
        vao[1].bind();
        buffers[1].bind();
        buffers[1].allocate(pos_delaunay->data(), pos_delaunay->size()*sizeof(float));
        poly_vertexLocation[0] = rendering_program.attributeLocation("vertex");
        rendering_program.enableAttributeArray(poly_vertexLocation[0]);
        rendering_program.setAttributeBuffer(poly_vertexLocation[0],GL_FLOAT,0,3);
        buffers[1].release();
        vao[1].release();

        //Voronoi Edges
        vao[2].bind();
        buffers[2].bind();
        buffers[2].allocate(pos_voronoi->data(), pos_voronoi->size()*sizeof(float));
        poly_vertexLocation[0] = rendering_program.attributeLocation("vertex");
        rendering_program.enableAttributeArray(poly_vertexLocation[0]);
        rendering_program.setAttributeBuffer(poly_vertexLocation[0],GL_FLOAT,0,3);
        buffers[2].release();
        vao[2].release();

        //Facets
        vao[3].bind();
        buffers[3].bind();
        buffers[3].allocate(pos_facets->data(), pos_facets->size()*sizeof(float));
        poly_vertexLocation[0] = rendering_program.attributeLocation("vertex");
        rendering_program.enableAttributeArray(poly_vertexLocation[0]);
        rendering_program.setAttributeBuffer(poly_vertexLocation[0],GL_FLOAT,0,3);
        buffers[3].release();
        vao[3].release();

        //New point
        vao[4].bind();
        buffers[4].bind();
        buffers[4].allocate(pos_newPoint->data(), pos_newPoint->size()*sizeof(float));
        poly_vertexLocation[0] = rendering_program.attributeLocation("vertex");
        rendering_program.enableAttributeArray(poly_vertexLocation[0]);
        rendering_program.setAttributeBuffer(poly_vertexLocation[0],GL_FLOAT,0,3);
        buffers[4].release();
        vao[4].release();
        //New Facet
        vao[5].bind();
        buffers[5].bind();
        buffers[5].allocate(pos_newFacet->data(), pos_newFacet->size()*sizeof(float));
        poly_vertexLocation[0] = rendering_program.attributeLocation("vertex");
        rendering_program.enableAttributeArray(poly_vertexLocation[0]);
        rendering_program.setAttributeBuffer(poly_vertexLocation[0],GL_FLOAT,0,3);
        buffers[5].release();
        vao[5].release();

        //Selected Point
        vao[6].bind();
        buffers[6].bind();
        buffers[6].allocate(pos_selectedVertex->data(), pos_selectedVertex->size()*sizeof(float));
        poly_vertexLocation[0] = rendering_program.attributeLocation("vertex");
        rendering_program.enableAttributeArray(poly_vertexLocation[0]);
        rendering_program.setAttributeBuffer(poly_vertexLocation[0],GL_FLOAT,0,3);
        buffers[6].release();
        vao[6].release();

        //Moving Point
        vao[7].bind();
        buffers[7].bind();
        buffers[7].allocate(pos_movingPoint->data(), pos_movingPoint->size()*sizeof(float));
        poly_vertexLocation[0] = rendering_program.attributeLocation("vertex");
        rendering_program.enableAttributeArray(poly_vertexLocation[0]);
        rendering_program.setAttributeBuffer(poly_vertexLocation[0],GL_FLOAT,0,3);
        buffers[7].release();
        vao[7].release();

        //Querry Point
        vao[8].bind();
        buffers[8].bind();
        buffers[8].allocate(pos_queryPoint->data(), pos_queryPoint->size()*sizeof(float));
        poly_vertexLocation[0] = rendering_program.attributeLocation("vertex");
        rendering_program.enableAttributeArray(poly_vertexLocation[0]);
        rendering_program.setAttributeBuffer(poly_vertexLocation[0],GL_FLOAT,0,3);
        buffers[8].release();
        vao[8].release();

        //Nearest Neighbor
        vao[9].bind();
        buffers[9].bind();
        buffers[9].allocate(pos_nearest_neighbor->data(), pos_nearest_neighbor->size()*sizeof(float));
        poly_vertexLocation[0] = rendering_program.attributeLocation("vertex");
        rendering_program.enableAttributeArray(poly_vertexLocation[0]);
        rendering_program.setAttributeBuffer(poly_vertexLocation[0],GL_FLOAT,0,3);
        buffers[9].release();
        vao[9].release();

        //Facet empty Sphere
        vao[10].bind();
        buffers[10].bind();
        buffers[10].allocate(pos_emptyFacet->data(), pos_emptyFacet->size()*sizeof(float));
        poly_vertexLocation[0] = rendering_program.attributeLocation("vertex");
        rendering_program.enableAttributeArray(poly_vertexLocation[0]);
        rendering_program.setAttributeBuffer(poly_vertexLocation[0],GL_FLOAT,0,3);
        buffers[10].release();
        vao[10].release();

        vao[21].bind();
        buffers[28].bind();
        buffers[28].allocate(incremental_next_point->data(), incremental_next_point->size()*sizeof(float));
        poly_vertexLocation[0] = rendering_program.attributeLocation("vertex");
        rendering_program.enableAttributeArray(poly_vertexLocation[0]);
        rendering_program.setAttributeBuffer(poly_vertexLocation[0],GL_FLOAT,0,3);
        buffers[28].release();
        vao[21].release();

        vao[22].bind();
        buffers[29].bind();
        buffers[29].allocate(incremental_facet->data(), incremental_facet->size()*sizeof(float));
        poly_vertexLocation[0] = rendering_program.attributeLocation("vertex");
        rendering_program.enableAttributeArray(poly_vertexLocation[0]);
        rendering_program.setAttributeBuffer(poly_vertexLocation[0],GL_FLOAT,0,3);
        buffers[29].release();
        vao[22].release();

        vao[23].bind();
        buffers[30].bind();
        buffers[30].allocate(incremental_conflict->data(), incremental_conflict->size()*sizeof(float));
        poly_vertexLocation[0] = rendering_program.attributeLocation("vertex");
        rendering_program.enableAttributeArray(poly_vertexLocation[0]);
        rendering_program.setAttributeBuffer(poly_vertexLocation[0],GL_FLOAT,0,3);
        buffers[30].release();
        vao[23].release();

        vao[24].bind();
        buffers[31].bind();
        buffers[31].allocate(incremental_points->data(), incremental_points->size()*sizeof(float));
        poly_vertexLocation[0] = rendering_program.attributeLocation("vertex");
        rendering_program.enableAttributeArray(poly_vertexLocation[0]);
        rendering_program.setAttributeBuffer(poly_vertexLocation[0],GL_FLOAT,0,3);
        buffers[31].release();
        vao[24].release();
    }

    rendering_program.release();
    rendering_program_spheres.bind();
    {
        //Track Ball
        vao[11].bind();
        buffers[11].bind();
        buffers[11].allocate(pos_trackBall->data(), pos_trackBall->size()*sizeof(float));
        centerLocation[0] = rendering_program_spheres.attributeLocation("center");
        rendering_program_spheres.enableAttributeArray(centerLocation[0]);
        rendering_program_spheres.setAttributeBuffer(centerLocation[0],GL_FLOAT,0,3);
        buffers[11].release();

        buffers[12].bind();
        buffers[12].allocate(normals_trackBall->data(), normals_trackBall->size()*sizeof(float));
        normalsLocation[0] = rendering_program_spheres.attributeLocation("normal");
        rendering_program_spheres.enableAttributeArray(normalsLocation[0]);
        rendering_program_spheres.setAttributeBuffer(normalsLocation[0],GL_FLOAT,0,3);
        buffers[12].release();

        buffers[13].bind();
        buffers[13].allocate(points_trackBall->data(), points_trackBall->size()*sizeof(float));
        poly_vertexLocation[1] = rendering_program_spheres.attributeLocation("vertex");
        rendering_program_spheres.enableAttributeArray(poly_vertexLocation[1]);
        rendering_program_spheres.setAttributeBuffer(poly_vertexLocation[1],GL_FLOAT,0,3);
        buffers[13].release();

        if(extension_is_found)
        {
            glVertexAttribDivisor(centerLocation[0],1);
            glVertexAttribDivisor(normalsLocation[0],0);
        }
        vao[11].release();

        //Empty Sphere
        vao[12].bind();
        if(extension_is_found)
        {
            buffers[14].bind();
            buffers[14].allocate(pos_emptySphere->data(), pos_emptySphere->size()*sizeof(float));
            centerLocation[0] = rendering_program_spheres.attributeLocation("center");
            rendering_program_spheres.enableAttributeArray(centerLocation[0]);
            rendering_program_spheres.setAttributeBuffer(centerLocation[0],GL_FLOAT,0,3);
            buffers[14].release();
        }
        buffers[32].bind();
        buffers[32].allocate(normals_emptySphere->data(), normals_emptySphere->size()*sizeof(float));
        normalsLocation[0] = rendering_program_spheres.attributeLocation("normal");
        rendering_program_spheres.enableAttributeArray(normalsLocation[0]);
        rendering_program_spheres.setAttributeBuffer(normalsLocation[0],GL_FLOAT,0,3);
        buffers[32].release();

        buffers[16].bind();
        buffers[16].allocate(points_emptySphere->data(), points_emptySphere->size()*sizeof(float));
        poly_vertexLocation[1] = rendering_program_spheres.attributeLocation("vertex");
        rendering_program_spheres.enableAttributeArray(poly_vertexLocation[1]);
        rendering_program_spheres.setAttributeBuffer(poly_vertexLocation[1],GL_FLOAT,0,3);
        buffers[16].release();

        if(extension_is_found)
        {
            glVertexAttribDivisor(centerLocation[0],1);
            glVertexAttribDivisor(normalsLocation[0],0);
        }
        vao[12].release();

        //Vertex Sphere
        vao[13].bind();
        buffers[0].bind();
        centerLocation[0] = rendering_program_spheres.attributeLocation("center");
        rendering_program_spheres.enableAttributeArray(centerLocation[0]);
        rendering_program_spheres.setAttributeBuffer(centerLocation[0],GL_FLOAT,0,3);
        buffers[0].release();

        buffers[15].bind();
        buffers[15].allocate(normals_sphere->data(), normals_sphere->size()*sizeof(float));
        normalsLocation[0] = rendering_program_spheres.attributeLocation("normal");
        rendering_program_spheres.enableAttributeArray(normalsLocation[0]);
        rendering_program_spheres.setAttributeBuffer(normalsLocation[0],GL_FLOAT,0,3);
        buffers[15].release();

        buffers[17].bind();
        buffers[17].allocate(points_sphere->data(), points_sphere->size()*sizeof(float));
        poly_vertexLocation[1] = rendering_program_spheres.attributeLocation("vertex");
        rendering_program_spheres.enableAttributeArray(poly_vertexLocation[1]);
        rendering_program_spheres.setAttributeBuffer(poly_vertexLocation[1],GL_FLOAT,0,3);
        buffers[17].release();

        if(extension_is_found)
        {
            glVertexAttribDivisor(centerLocation[0],1);
            glVertexAttribDivisor(normalsLocation[0],0);
        }
        vao[13].release();

        //New point Sphere
        vao[14].bind();
        buffers[4].bind();
        centerLocation[0] = rendering_program_spheres.attributeLocation("center");
        rendering_program_spheres.enableAttributeArray(centerLocation[0]);
        rendering_program_spheres.setAttributeBuffer(centerLocation[0],GL_FLOAT,0,3);
        buffers[4].release();

        buffers[15].bind();
        normalsLocation[0] = rendering_program_spheres.attributeLocation("normal");
        rendering_program_spheres.enableAttributeArray(normalsLocation[0]);
        rendering_program_spheres.setAttributeBuffer(normalsLocation[0],GL_FLOAT,0,3);
        buffers[15].release();

        buffers[17].bind();
        poly_vertexLocation[1] = rendering_program_spheres.attributeLocation("vertex");
        rendering_program_spheres.enableAttributeArray(poly_vertexLocation[1]);
        rendering_program_spheres.setAttributeBuffer(poly_vertexLocation[1],GL_FLOAT,0,3);
        buffers[17].release();

        if(extension_is_found)
        {
            glVertexAttribDivisor(centerLocation[0],1);
            glVertexAttribDivisor(normalsLocation[0],0);
        }
        vao[14].release();

        //Selected point Sphere
        vao[15].bind();
        buffers[6].bind();
        centerLocation[0] = rendering_program_spheres.attributeLocation("center");
        rendering_program_spheres.enableAttributeArray(centerLocation[0]);
        rendering_program_spheres.setAttributeBuffer(centerLocation[0],GL_FLOAT,0,3);
        buffers[6].release();

        buffers[15].bind();
        normalsLocation[0] = rendering_program_spheres.attributeLocation("normal");
        rendering_program_spheres.enableAttributeArray(normalsLocation[0]);
        rendering_program_spheres.setAttributeBuffer(normalsLocation[0],GL_FLOAT,0,3);
        buffers[15].release();

        buffers[17].bind();
        poly_vertexLocation[1] = rendering_program_spheres.attributeLocation("vertex");
        rendering_program_spheres.enableAttributeArray(poly_vertexLocation[1]);
        rendering_program_spheres.setAttributeBuffer(poly_vertexLocation[1],GL_FLOAT,0,3);
        buffers[17].release();

        if(extension_is_found)
        {
            glVertexAttribDivisor(centerLocation[0],1);
            glVertexAttribDivisor(normalsLocation[0],0);
        }
        vao[15].release();

        //Moving point Sphere
        vao[16].bind();
        buffers[7].bind();
        centerLocation[0] = rendering_program_spheres.attributeLocation("center");
        rendering_program_spheres.enableAttributeArray(centerLocation[0]);
        rendering_program_spheres.setAttributeBuffer(centerLocation[0],GL_FLOAT,0,3);
        buffers[7].release();

        buffers[15].bind();
        normalsLocation[0] = rendering_program_spheres.attributeLocation("normal");
        rendering_program_spheres.enableAttributeArray(normalsLocation[0]);
        rendering_program_spheres.setAttributeBuffer(normalsLocation[0],GL_FLOAT,0,3);
        buffers[15].release();

        buffers[17].bind();
        poly_vertexLocation[1] = rendering_program_spheres.attributeLocation("vertex");
        rendering_program_spheres.enableAttributeArray(poly_vertexLocation[1]);
        rendering_program_spheres.setAttributeBuffer(poly_vertexLocation[1],GL_FLOAT,0,3);
        buffers[17].release();

        if(extension_is_found)
        {
            glVertexAttribDivisor(centerLocation[0],1);
            glVertexAttribDivisor(normalsLocation[0],0);
        }
        vao[16].release();

        //Querry point Sphere
        vao[17].bind();
        buffers[8].bind();
        centerLocation[0] = rendering_program_spheres.attributeLocation("center");
        rendering_program_spheres.enableAttributeArray(centerLocation[0]);
        rendering_program_spheres.setAttributeBuffer(centerLocation[0],GL_FLOAT,0,3);
        buffers[8].release();

        buffers[15].bind();
        normalsLocation[0] = rendering_program_spheres.attributeLocation("normal");
        rendering_program_spheres.enableAttributeArray(normalsLocation[0]);
        rendering_program_spheres.setAttributeBuffer(normalsLocation[0],GL_FLOAT,0,3);
        buffers[15].release();

        buffers[17].bind();
        poly_vertexLocation[1] = rendering_program_spheres.attributeLocation("vertex");
        rendering_program_spheres.enableAttributeArray(poly_vertexLocation[1]);
        rendering_program_spheres.setAttributeBuffer(poly_vertexLocation[1],GL_FLOAT,0,3);
        buffers[17].release();

        if(extension_is_found)
        {
            glVertexAttribDivisor(centerLocation[0],1);
            glVertexAttribDivisor(normalsLocation[0],0);
        }
        vao[17].release();

        //Nearest neighbor Sphere
        vao[18].bind();
        buffers[9].bind();
        centerLocation[0] = rendering_program_spheres.attributeLocation("center");
        rendering_program_spheres.enableAttributeArray(centerLocation[0]);
        rendering_program_spheres.setAttributeBuffer(centerLocation[0],GL_FLOAT,0,3);
        buffers[9].release();

        buffers[15].bind();
        normalsLocation[0] = rendering_program_spheres.attributeLocation("normal");
        rendering_program_spheres.enableAttributeArray(normalsLocation[0]);
        rendering_program_spheres.setAttributeBuffer(normalsLocation[0],GL_FLOAT,0,3);
        buffers[15].release();

        buffers[17].bind();
        poly_vertexLocation[1] = rendering_program_spheres.attributeLocation("vertex");
        rendering_program_spheres.enableAttributeArray(poly_vertexLocation[1]);
        rendering_program_spheres.setAttributeBuffer(poly_vertexLocation[1],GL_FLOAT,0,3);
        buffers[17].release();

        if(extension_is_found)
        {
            glVertexAttribDivisor(centerLocation[0],1);
            glVertexAttribDivisor(normalsLocation[0],0);
        }
        vao[18].release();

        //incremental list
        vao[25].bind();
        buffers[31].bind();
        rendering_program_spheres.enableAttributeArray(centerLocation[0]);
        rendering_program_spheres.setAttributeBuffer(centerLocation[0],GL_FLOAT,0,3);
        buffers[31].release();

        buffers[15].bind();
        normalsLocation[0] = rendering_program_spheres.attributeLocation("normal");
        rendering_program_spheres.enableAttributeArray(normalsLocation[0]);
        rendering_program_spheres.setAttributeBuffer(normalsLocation[0],GL_FLOAT,0,3);
        buffers[15].release();

        buffers[17].bind();
        poly_vertexLocation[1] = rendering_program_spheres.attributeLocation("vertex");
        rendering_program_spheres.enableAttributeArray(poly_vertexLocation[1]);
        rendering_program_spheres.setAttributeBuffer(poly_vertexLocation[1],GL_FLOAT,0,3);
        buffers[17].release();

        if(extension_is_found)
        {
            glVertexAttribDivisor(centerLocation[0],1);
            glVertexAttribDivisor(normalsLocation[0],0);
        }
        vao[25].release();

        //incremental next point
        vao[26].bind();
        buffers[28].bind();
        rendering_program_spheres.enableAttributeArray(centerLocation[0]);
        rendering_program_spheres.setAttributeBuffer(centerLocation[0],GL_FLOAT,0,3);
        buffers[28].release();

        buffers[15].bind();
        normalsLocation[0] = rendering_program_spheres.attributeLocation("normal");
        rendering_program_spheres.enableAttributeArray(normalsLocation[0]);
        rendering_program_spheres.setAttributeBuffer(normalsLocation[0],GL_FLOAT,0,3);
        buffers[15].release();

        buffers[17].bind();
        poly_vertexLocation[1] = rendering_program_spheres.attributeLocation("vertex");
        rendering_program_spheres.enableAttributeArray(poly_vertexLocation[1]);
        rendering_program_spheres.setAttributeBuffer(poly_vertexLocation[1],GL_FLOAT,0,3);
        buffers[17].release();

        if(extension_is_found)
        {
            glVertexAttribDivisor(centerLocation[0],1);
            glVertexAttribDivisor(normalsLocation[0],0);
        }
        vao[26].release();

    }
    rendering_program_spheres.release();
    rendering_program_cylinders.bind();
    {
        vao[19].bind();
        buffers[18].bind();
        buffers[18].allocate(points_cylinder->data(), points_cylinder->size()*sizeof(float));
        poly_vertexLocation[2] = rendering_program_cylinders.attributeLocation("vertex");
        rendering_program_cylinders.enableAttributeArray(poly_vertexLocation[2]);
        rendering_program_cylinders.setAttributeBuffer(poly_vertexLocation[2],GL_FLOAT,0,3);
        buffers[18].release();

        buffers[19].bind();
        buffers[19].allocate(normals_cylinder->data(), normals_cylinder->size()*sizeof(float));
        normalsLocation[1] = rendering_program_cylinders.attributeLocation("normal");
        rendering_program_cylinders.enableAttributeArray(normalsLocation[1]);
        rendering_program_cylinders.setAttributeBuffer(normalsLocation[1],GL_FLOAT,0,3);
        buffers[19].release();

        buffers[20].bind();
        buffers[20].allocate(transfo1_delaunay->data(), transfo1_delaunay->size()*sizeof(float));
        centerLocation[1] = rendering_program_cylinders.attributeLocation("transfo1");
        rendering_program_cylinders.enableAttributeArray(centerLocation[1]);
        rendering_program_cylinders.setAttributeBuffer(centerLocation[1],GL_FLOAT,0,4);
        buffers[20].release();

        buffers[21].bind();
        buffers[21].allocate(transfo2_delaunay->data(), transfo2_delaunay->size()*sizeof(float));
        centerLocation[2] = rendering_program_cylinders.attributeLocation("transfo2");
        rendering_program_cylinders.enableAttributeArray(centerLocation[2]);
        rendering_program_cylinders.setAttributeBuffer(centerLocation[2],GL_FLOAT,0,4);
        buffers[21].release();


        buffers[22].bind();
        buffers[22].allocate(transfo3_delaunay->data(), transfo3_delaunay->size()*sizeof(float));
        centerLocation[3] = rendering_program_cylinders.attributeLocation("transfo3");
        rendering_program_cylinders.enableAttributeArray(centerLocation[3]);
        rendering_program_cylinders.setAttributeBuffer(centerLocation[3],GL_FLOAT,0,4);
        buffers[22].release();

        buffers[23].bind();
        buffers[23].allocate(transfo4_delaunay->data(), transfo4_delaunay->size()*sizeof(float));
        centerLocation[4] = rendering_program_cylinders.attributeLocation("transfo4");
        rendering_program_cylinders.enableAttributeArray(centerLocation[4]);
        rendering_program_cylinders.setAttributeBuffer(centerLocation[4],GL_FLOAT,0,4);
        buffers[23].release();

        if(extension_is_found)
        {
            glVertexAttribDivisor(centerLocation[1],1);
            glVertexAttribDivisor(centerLocation[2],1);
            glVertexAttribDivisor(centerLocation[3],1);
            glVertexAttribDivisor(centerLocation[4],1);
            glVertexAttribDivisor(normalsLocation[1],0);
        }

        vao[19].release();

        vao[20].bind();
        buffers[18].bind();
        poly_vertexLocation[2] = rendering_program_cylinders.attributeLocation("vertex");
        rendering_program_cylinders.enableAttributeArray(poly_vertexLocation[2]);
        rendering_program_cylinders.setAttributeBuffer(poly_vertexLocation[2],GL_FLOAT,0,3);
        buffers[18].release();

        buffers[19].bind();
        normalsLocation[1] = rendering_program_cylinders.attributeLocation("normal");
        rendering_program_cylinders.enableAttributeArray(normalsLocation[1]);
        rendering_program_cylinders.setAttributeBuffer(normalsLocation[1],GL_FLOAT,0,3);
        buffers[19].release();

        buffers[24].bind();
        buffers[24].allocate(transfo1_voronoi->data(), transfo1_voronoi->size()*sizeof(float));
        centerLocation[1] = rendering_program_cylinders.attributeLocation("transfo1");
        rendering_program_cylinders.enableAttributeArray(centerLocation[1]);
        rendering_program_cylinders.setAttributeBuffer(centerLocation[1],GL_FLOAT,0,4);
        buffers[24].release();

        buffers[25].bind();
        buffers[25].allocate(transfo2_voronoi->data(), transfo2_voronoi->size()*sizeof(float));
        centerLocation[2] = rendering_program_cylinders.attributeLocation("transfo2");
        rendering_program_cylinders.enableAttributeArray(centerLocation[2]);
        rendering_program_cylinders.setAttributeBuffer(centerLocation[2],GL_FLOAT,0,4);
        buffers[25].release();


        buffers[26].bind();
        buffers[26].allocate(transfo3_voronoi->data(), transfo3_voronoi->size()*sizeof(float));
        centerLocation[3] = rendering_program_cylinders.attributeLocation("transfo3");
        rendering_program_cylinders.enableAttributeArray(centerLocation[3]);
        rendering_program_cylinders.setAttributeBuffer(centerLocation[3],GL_FLOAT,0,4);
        buffers[26].release();

        buffers[27].bind();
        buffers[27].allocate(transfo4_voronoi->data(), transfo4_voronoi->size()*sizeof(float));
        centerLocation[4] = rendering_program_cylinders.attributeLocation("transfo4");
        rendering_program_cylinders.enableAttributeArray(centerLocation[4]);
        rendering_program_cylinders.setAttributeBuffer(centerLocation[4],GL_FLOAT,0,4);
        buffers[27].release();
        if(extension_is_found)
        {
            glVertexAttribDivisor(centerLocation[1],1);
            glVertexAttribDivisor(centerLocation[2],1);
            glVertexAttribDivisor(centerLocation[3],1);
            glVertexAttribDivisor(centerLocation[4],1);
            glVertexAttribDivisor(normalsLocation[1],0);
        }

        vao[20].release();
    }
    rendering_program_cylinders.release();
    are_buffers_initialized = true;
}

void Viewer::attrib_buffers(QGLViewer* viewer)
{
    QMatrix4x4 mvpMatrix;
    QMatrix4x4 mvMatrix;
    double mat[16];
    viewer->camera()->getModelViewProjectionMatrix(mat);
    for(int i=0; i < 16; i++)
    {
        mvpMatrix.data()[i] = (float)mat[i];
    }
    viewer->camera()->getModelViewMatrix(mat);
    for(int i=0; i < 16; i++)
    {
        mvMatrix.data()[i] = (float)mat[i];
    }
    QVector4D	position(0.0f,0.0f,1.0f,1.0f );
    // Ambient
    ambient[0] = 0.29225f;
    ambient[1] = 0.29225f;
    ambient[2] = 0.29225f;
    ambient[3] = 1.0f;
    // Diffuse
    diffuse[0] = 0.50754f;
    diffuse[1] = 0.50754f;
    diffuse[2] = 0.50754f;
    diffuse[3] = 1.0f;
    // Specular
    specular[0] = 0.508273f;
    specular[1] = 0.508273f;
    specular[2] = 0.508273f;
    specular[3] = 1.0f;
    // Shininess
    shininess = 51.2f;

    rendering_program.bind();
    mvpLocation[0] = rendering_program.uniformLocation("mvp_matrix");
    colorLocation[0] = rendering_program.uniformLocation("color");
    rendering_program.setUniformValue(mvpLocation[0], mvpMatrix);

    rendering_program.release();


    rendering_program_spheres.bind();
    colorLocation[1] = rendering_program_spheres.uniformLocation("color");
    mvpLocation[1] = rendering_program_spheres.uniformLocation("mvp_matrix");
    mvLocation[0] = rendering_program_spheres.uniformLocation("mv_matrix");
    lightLocation[0] = rendering_program_spheres.uniformLocation("light_pos");
    lightLocation[1] = rendering_program_spheres.uniformLocation("light_diff");
    lightLocation[2] = rendering_program_spheres.uniformLocation("light_spec");
    lightLocation[3] = rendering_program_spheres.uniformLocation("light_amb");
    lightLocation[4] = rendering_program_spheres.uniformLocation("spec_power");

    rendering_program_spheres.setUniformValue(lightLocation[0], position);
    rendering_program_spheres.setUniformValue(mvpLocation[1], mvpMatrix);
    rendering_program_spheres.setUniformValue(mvLocation[0], mvMatrix);
    rendering_program_spheres.setUniformValue(lightLocation[1], diffuse);
    rendering_program_spheres.setUniformValue(lightLocation[2], specular);
    rendering_program_spheres.setUniformValue(lightLocation[3], ambient);
    rendering_program_spheres.setUniformValue(lightLocation[4], shininess);

    rendering_program_spheres.release();

    rendering_program_cylinders.bind();
    colorLocation[2] = rendering_program_cylinders.uniformLocation("color");
    mvpLocation[2] = rendering_program_cylinders.uniformLocation("mvp_matrix");
    mvLocation[1] = rendering_program_cylinders.uniformLocation("mv_matrix");
    lightLocation[5] = rendering_program_cylinders.uniformLocation("light_pos");
    lightLocation[6] = rendering_program_cylinders.uniformLocation("light_diff");
    lightLocation[7] = rendering_program_cylinders.uniformLocation("light_spec");
    lightLocation[8] = rendering_program_cylinders.uniformLocation("light_amb");
    lightLocation[9] = rendering_program_cylinders.uniformLocation("spec_power");

    rendering_program_cylinders.setUniformValue(lightLocation[5], position);
    rendering_program_cylinders.setUniformValue(lightLocation[6], diffuse);
    rendering_program_cylinders.setUniformValue(lightLocation[7], specular);
    rendering_program_cylinders.setUniformValue(lightLocation[8], ambient);
    rendering_program_cylinders.setUniformValue(lightLocation[9], shininess);
    rendering_program_cylinders.setUniformValue(mvpLocation[2], mvpMatrix);
    rendering_program_cylinders.setUniformValue(mvLocation[1], mvMatrix);


    rendering_program_cylinders.release();
}


/*************************************************************/
/*  Draw functions */

void Viewer::draw()
{
    glEnable(GL_DEPTH_TEST);
    if(!are_buffers_initialized)
        initialize_buffers();
    QFont fontPrompt("Arial", 8);
    attrib_buffers(this);
    if(m_isFlat || !extension_is_found)
    {

        if(m_showVertex)
        {
            rendering_program.bind();
            glPointSize(8.0);
            vao[0].bind();
            rendering_program.setUniformValue(colorLocation[0], m_colorVertex);
            glDrawArrays(GL_POINTS, 0, pos_points->size()/3);
            vao[0].release();
            rendering_program.release();
        }
        if(m_showDEdge)
        {
            rendering_program.bind();
            vao[1].bind();
            rendering_program.setUniformValue(colorLocation[0], m_colorDEdge);
            glDrawArrays(GL_LINES, 0, pos_delaunay->size()/3);
            vao[1].release();
            rendering_program.release();
        }
        if(m_showVEdge)
        {
            rendering_program.bind();
            vao[2].bind();
            rendering_program.setUniformValue(colorLocation[0], m_colorVEdge);
            glDrawArrays(GL_LINES, 0, pos_voronoi->size()/3);
            vao[2].release();
            rendering_program.release();
        }
        // Insert point mode
        if( m_curMode == INSERT_PT) {
            // Show prompt messages
            qglColor( ::Qt::black );
            drawText( width()-200, 20, tr("Shift+Left: Insert a point"), fontPrompt );
            drawText( width()-200, 40, tr("Hold Left: Move the point"), fontPrompt );
            drawText( width()-200, 60, tr("Return: Insert to DT"), fontPrompt );
            drawText( width()-200, 80, tr("Escape: Cancel insertion"), fontPrompt );
            drawText( width()-200, 100, tr("Shift+Wheel: Resize trackball"), fontPrompt );
            if( m_hasNewPt ) {
                rendering_program.bind();
                vao[4].bind();
                color.setRgbF(1.0,0.0,0.0);
                rendering_program.setUniformValue(colorLocation[0], color);
                glDrawArrays(GL_POINTS, 0, pos_newPoint->size()/3);
                vao[4].release();

                vao[5].bind();
                color.setRgb(215, 80, 0, 96);
                rendering_program.setUniformValue(colorLocation[0], color);
                glDrawArrays(GL_TRIANGLES, 0, pos_newFacet->size()/3);
                vao[5].release();
                rendering_program.release();
            }
        }
        else if( m_curMode == SELECT) {
            // Show prompt messages
            qglColor( ::Qt::black );
            drawText( width()-200, 20, tr("Shift+Left: Select"), fontPrompt );
            drawText( width()-200, 40, tr("Ctrl+Left: Add selection"),
                      QFont("Arial", 14) );
            drawText( width()-200, 60, tr("Escape: Cancel selection"), fontPrompt );
            drawText( width()-200, 80, tr("DEL: Delete selected"), fontPrompt );
            rendering_program.bind();
            vao[6].bind();
            color.setRgbF(1.0,0.0,0.0);
            rendering_program.setUniformValue(colorLocation[0], color);
            glDrawArrays(GL_POINTS, 0, pos_selectedVertex->size()/3);
            vao[6].release();
            rendering_program.release();
        }
        else if( m_curMode == MOVE ) {
            // Show prompt messages
            qglColor( ::Qt::black );
            drawText( width()-200, 20, tr("Left Click: Select"), fontPrompt );
            if(m_isMoving)
                drawText( width()-200, 40, tr("Shift+Wheel: Resize trackball"), fontPrompt );
            rendering_program.bind();
            vao[7].bind();
            color.setRgbF(1.0,0.0,0.0);
            rendering_program.setUniformValue(colorLocation[0], color);
            glDrawArrays(GL_POINTS, 0, pos_movingPoint->size()/3);
            vao[7].release();
            rendering_program.release();
        }
        else if( m_curMode == FINDNB ) {
            // Show prompt messages
            qglColor( ::Qt::black );
            drawText( width()-200, 20, tr("Shift+Left: Place query point"), fontPrompt );
            drawText( width()-200, 40, tr("Shift+Wheel: Resize trackball"), fontPrompt );
            rendering_program.bind();
            vao[8].bind();
            color.setRgbF(1.0,0.0,0.0);
            rendering_program.setUniformValue(colorLocation[0], color);
            glDrawArrays(GL_POINTS, 0, pos_queryPoint->size()/3);
            vao[8].release();

            vao[9].bind();
            color.setRgbF(1.0,0.0,0.0);
            rendering_program.setUniformValue(colorLocation[0], color);
            glDrawArrays(GL_POINTS, 0, pos_nearest_neighbor->size()/3);
            vao[9].release();
            rendering_program.release();
        }
        else if(m_curMode == EMPTYSPH){
            // Show prompt messages
            qglColor( ::Qt::black );
            drawText( width()-200, 20, tr("Shift+Left: Place query point"), fontPrompt );
            drawText( width()-200, 40, tr("Press S: Show/Hide trackball"), fontPrompt );
            drawText( width()-200, 60, tr("Shift+Wheel: Resize trackball"), fontPrompt );
            rendering_program.bind();
            vao[8].bind();
            color.setRgbF(1.0,0.0,0.0);
            rendering_program.setUniformValue(colorLocation[0], color);
            glDrawArrays(GL_POINTS, 0, pos_queryPoint->size()/3);
            vao[8].release();
            vao[10].bind();
            rendering_program.setUniformValue(colorLocation[0], m_colorFacet);
            glDrawArrays(GL_TRIANGLES, 0, pos_emptyFacet->size()/3);
            vao[10].release();
            rendering_program.release();
        }
        // Draw all points during incremental mode
        if( !m_incrementalPts.isEmpty() ) {
            // draw the rest to-be-inserted vertices
            rendering_program.bind();
            glPointSize(8.0);
            vao[24].bind();
            color.setRgbF(0.7,0.7,0.7);
            rendering_program.setUniformValue(colorLocation[0],color);
            glDrawArrays(GL_POINTS, 0, incremental_points->size()/3);
            vao[24].release();
            rendering_program.release();
            switch( m_curStep ) {
            case NEWPT:
                // Show prompt messages
                qglColor( ::Qt::black );
                drawText( 10, 20, tr("Highlight the next-to-insert point"), fontPrompt );
                // Highlight the next-to-insert point
                rendering_program.bind();
                glPointSize(8.0);
                vao[21].bind();
                color.setRgbF(1.0,0.0,0.0);
                rendering_program.setUniformValue(colorLocation[0], color);
                glDrawArrays(GL_POINTS, 0, incremental_next_point->size()/3);
                vao[21].release();
                rendering_program.release();
                break;
            case CELL:  // show the tetrahedron that contains the point
                // Show prompt messages
                qglColor( ::Qt::black );
                drawText( 10, 20, tr("Show the tetrahedron containing the point"), fontPrompt );
                drawText( 10, 40, tr("(Only finite facets are drawn)"), fontPrompt );
                // Highlight the next-to-insert vertex
                rendering_program.bind();
                glPointSize(8.0);
                vao[21].bind();
                color.setRgbF(1.0,0.0,0.0);
                rendering_program.setUniformValue(colorLocation[0],  color);
                glDrawArrays(GL_POINTS, 0, incremental_next_point->size()/3);
                vao[21].release();
                rendering_program.release();
                // Draw the cell containing that point
                rendering_program.bind();
                vao[22].bind();
                rendering_program.setUniformValue(colorLocation[0], m_colorFacet);
                glDrawArrays(GL_TRIANGLES, 0, incremental_facet->size()/3);
                vao[22].release();
                rendering_program.release();
                break;
            case CONFLICT:  // show the conflict region
                // Show prompt messages
                qglColor( ::Qt::black );
                drawText( 10, 20, tr("Show the conflict region"), fontPrompt );
                // Highlight the next-to-insert vertex
                rendering_program.bind();
                glPointSize(8.0);
                vao[21].bind();
                color.setRgbF(1.0,0.0,0.0);
                rendering_program.setUniformValue(colorLocation[0], color);
                glDrawArrays(GL_POINTS, 0, incremental_next_point->size()/3);
                vao[21].release();
                rendering_program.release();
                // Draw conflict region
                rendering_program.bind();
                vao[23].bind();
                color.setRgb(215, 80, 0, 96);
                rendering_program.setUniformValue(colorLocation[0], color);
                glDrawArrays(GL_TRIANGLES, 0, incremental_facet->size()/3);
                vao[23].release();
                rendering_program.release();

                break;
            default:
                break;
            }//end-of=switch
        }//end-if-incpts

    }
    else
    {
        if(m_showVertex)
        {
            rendering_program_spheres.bind();
            vao[13].bind();
            rendering_program_spheres.setUniformValue(colorLocation[1], m_colorVertex);
            glDrawArraysInstanced(GL_TRIANGLES, 0, points_sphere->size()/3, pos_points->size()/3);
            vao[13].release();
            rendering_program.release();
        }
        if(m_showDEdge)
        {
            rendering_program_cylinders.bind();
            vao[19].bind();
            rendering_program_cylinders.setUniformValue(colorLocation[2], m_colorDEdge);
            glDrawArraysInstanced(GL_TRIANGLES, 0, points_cylinder->size()/3, transfo1_delaunay->size()/4);
            vao[19].release();
            rendering_program_cylinders.release();
        }
        if(m_showVEdge)
        {
            rendering_program_cylinders.bind();
            vao[20].bind();
            rendering_program_cylinders.setUniformValue(colorLocation[2], m_colorVEdge);
            glDrawArraysInstanced(GL_TRIANGLES, 0, points_cylinder->size()/3, transfo1_voronoi->size()/4);
            vao[20].release();
            rendering_program_cylinders.release();
        }

        if( m_curMode == INSERT_PT) {
            // Show prompt messages
            qglColor( ::Qt::black );
            drawText( width()-200, 20, tr("Shift+Left: Insert a point"), fontPrompt );
            drawText( width()-200, 40, tr("Hold Left: Move the point"), fontPrompt );
            drawText( width()-200, 60, tr("Return: Insert to DT"), fontPrompt );
            drawText( width()-200, 80, tr("Escape: Cancel insertion"), fontPrompt );
            drawText( width()-200, 100, tr("Shift+Wheel: Resize trackball"), fontPrompt );
            if( m_hasNewPt ) {
                rendering_program_spheres.bind();
                vao[14].bind();
                color.setRgbF(1.0,0.0,0.0);
                rendering_program_spheres.setUniformValue(colorLocation[1], color);
                glDrawArraysInstanced(GL_TRIANGLES, 0, points_sphere->size()/3, pos_newPoint->size()/3);
                vao[14].release();
                rendering_program_spheres.release();
                rendering_program.bind();
                vao[5].bind();
                color.setRgb(215, 80, 0, 96);
                rendering_program.setUniformValue(colorLocation[0], color);
                glDrawArrays(GL_TRIANGLES, 0, pos_newFacet->size()/3);
                vao[5].release();
                rendering_program.release();
            }
        }
        else if( m_curMode == SELECT) {
            // Show prompt messages
            qglColor( ::Qt::black );
            drawText( width()-200, 20, tr("Shift+Left: Select"), fontPrompt );
            drawText( width()-200, 40, tr("Ctrl+Left: Add selection"),
                      QFont("Arial", 14) );
            drawText( width()-200, 60, tr("Escape: Cancel selection"), fontPrompt );
            drawText( width()-200, 80, tr("DEL: Delete selected"), fontPrompt );
            rendering_program_spheres.bind();
            vao[15].bind();
            color.setRgbF(1.0,0.0,0.0);
            rendering_program_spheres.setUniformValue(colorLocation[1], color);
            glDrawArraysInstanced(GL_TRIANGLES, 0, points_sphere->size()/3, pos_selectedVertex->size()/3);
            vao[15].release();
            rendering_program_spheres.release();
        }
        else if( m_curMode == MOVE ) {
            // Show prompt messages
            qglColor( ::Qt::black );
            drawText( width()-200, 20, tr("Left Click: Select"), fontPrompt );
            if(m_isMoving)
                drawText( width()-200, 40, tr("Shift+Wheel: Resize trackball"), fontPrompt );
            rendering_program_spheres.bind();
            vao[16].bind();
            color.setRgbF(1.0,0.0,0.0);
            rendering_program_spheres.setUniformValue(colorLocation[1], color);
            glDrawArraysInstanced(GL_TRIANGLES, 0, points_sphere->size()/3, pos_movingPoint->size()/3);
            vao[16].release();
            rendering_program_spheres.release();
        }
        else if( m_curMode == FINDNB ) {
            // Show prompt messages
            qglColor( ::Qt::black );
            drawText( width()-200, 20, tr("Shift+Left: Place query point"), fontPrompt );
            drawText( width()-200, 40, tr("Shift+Wheel: Resize trackball"), fontPrompt );
            rendering_program_spheres.bind();
            vao[17].bind();
            color.setRgbF(1.0,0.0,0.0);
            rendering_program_spheres.setUniformValue(colorLocation[1], color);
            glDrawArraysInstanced(GL_TRIANGLES, 0, points_sphere->size()/3, pos_queryPoint->size()/3);
            vao[17].release();
            rendering_program_spheres.release();
            rendering_program_spheres.bind();
            vao[18].bind();
            color.setRgbF(1.0,0.0,0.0);
            rendering_program_spheres.setUniformValue(colorLocation[1], color);
            glDrawArraysInstanced(GL_TRIANGLES, 0, points_sphere->size()/3, pos_nearest_neighbor->size()/3);
            vao[18].release();
            rendering_program_spheres.release();
        }
        else if(m_curMode == EMPTYSPH){
            // Show prompt messages
            qglColor( ::Qt::black );
            drawText( width()-200, 20, tr("Shift+Left: Place query point"), fontPrompt );
            drawText( width()-200, 40, tr("Press S: Show/Hide trackball"), fontPrompt );
            drawText( width()-200, 60, tr("Shift+Wheel: Resize trackball"), fontPrompt );
            rendering_program_spheres.bind();
            vao[17].bind();
            color.setRgbF(1.0,0.0,0.0);
            rendering_program_spheres.setUniformValue(colorLocation[1], color);
            glDrawArraysInstanced(GL_TRIANGLES, 0, points_sphere->size()/3, pos_queryPoint->size()/3);
            vao[17].release();
            rendering_program_spheres.release();

            rendering_program.bind();
            vao[10].bind();
            rendering_program.setUniformValue(colorLocation[0], m_colorFacet);
            glDrawArrays(GL_TRIANGLES, 0, pos_emptyFacet->size()/3);
            vao[10].release();
            rendering_program.release();
        }
        // Draw all points during incremental mode
        if( !m_incrementalPts.isEmpty() ) {
            // draw the rest to-be-inserted vertices
            rendering_program_spheres.bind();
            vao[25].bind();
            color.setRgbF(0.7,0.7,0.7);
            rendering_program_spheres.setUniformValue(colorLocation[1],color);
            glDrawArraysInstanced(GL_TRIANGLES, 0, points_sphere->size()/3, incremental_points->size()/3);
            vao[25].release();
            rendering_program_spheres.release();
            switch( m_curStep ) {
            case NEWPT:
                // Show prompt messages
                qglColor( ::Qt::black );
                drawText( 10, 20, tr("Highlight the next-to-insert point"), fontPrompt );
                // Highlight the next-to-insert point
                rendering_program_spheres.bind();
                vao[26].bind();
                color.setRgbF(1.0,0.0,0.0);
                rendering_program_spheres.setUniformValue(colorLocation[1],color);
                glDrawArraysInstanced(GL_TRIANGLES, 0, points_sphere->size()/3, incremental_next_point->size()/3);
                vao[26].release();
                rendering_program_spheres.release();
                break;
            case CELL:  // show the tetrahedron that contains the point
                // Show prompt messages
                qglColor( ::Qt::black );
                drawText( 10, 20, tr("Show the tetrahedron containing the point"), fontPrompt );
                drawText( 10, 40, tr("(Only finite facets are drawn)"), fontPrompt );
                // Highlight the next-to-insert vertex
                rendering_program_spheres.bind();
                vao[26].bind();
                color.setRgbF(1.0,0.0,0.0);
                rendering_program_spheres.setUniformValue(colorLocation[1],  color);
                glDrawArraysInstanced(GL_TRIANGLES, 0, points_sphere->size()/3, incremental_next_point->size()/3);
                vao[26].release();
                rendering_program_spheres.release();
                // Draw the cell containing that point
                rendering_program.bind();
                vao[22].bind();
                rendering_program.setUniformValue(colorLocation[0], m_colorFacet);
                glDrawArrays(GL_TRIANGLES, 0, incremental_facet->size()/3);
                vao[22].release();
                rendering_program.release();
                break;
            case CONFLICT:  // show the conflict region
                // Show prompt messages
                qglColor( ::Qt::black );
                drawText( 10, 20, tr("Show the conflict region"), fontPrompt );
                // Highlight the next-to-insert vertex
                rendering_program_spheres.bind();
                vao[26].bind();
                color.setRgbF(1.0,0.0,0.0);
                rendering_program_spheres.setUniformValue(colorLocation[1],  color);
                glDrawArraysInstanced(GL_TRIANGLES, 0, points_sphere->size()/3, incremental_next_point->size()/3);
                vao[26].release();
                rendering_program_spheres.release();
                // Draw conflict region
                rendering_program.bind();
                vao[23].bind();
                color.setRgb(215, 80, 0, 96);
                rendering_program.setUniformValue(colorLocation[0], color);
                glDrawArrays(GL_TRIANGLES, 0, incremental_facet->size()/3);
                vao[23].release();
                rendering_program.release();

                break;
            default:
                break;
            }//end-of=switch
        }//end-if-incpts
    }

    if(m_showFacet)
    {
        rendering_program.bind();
        vao[3].bind();
        rendering_program.setUniformValue(colorLocation[0], m_colorFacet);
        glDrawArrays(GL_TRIANGLES, 0, pos_facets->size()/3);
        vao[3].release();
        rendering_program.release();
    }
    if( m_curMode == INSERT_V  ) {
        // Show prompt messages
        qglColor( ::Qt::black );
        drawText( width()-200, 20, tr("Shift+Left: Insert a vertex"), fontPrompt );
        drawText( width()-200, 40, tr("Shift+Wheel: Resize trackball"), fontPrompt );

    }
    if(m_curMode != NONE && m_curMode !=  SELECT && m_showTrackball)
    {

        rendering_program_spheres.bind();
        vao[11].bind();
        rendering_program_spheres.setUniformValue(colorLocation[1], m_colorTrackball);
        if(extension_is_found)
            glDrawArraysInstanced(GL_TRIANGLES, 0, points_trackBall->size()/3, 1);
        else
            glDrawArrays(GL_TRIANGLES, 0, points_trackBall->size()/3);
        vao[11].release();
        rendering_program_spheres.release();
    }
    if(m_curMode ==EMPTYSPH)
    {
        rendering_program_spheres.bind();
        vao[12].bind();
        rendering_program_spheres.setUniformValue(colorLocation[1], m_colorEmptySphere);
        if(extension_is_found)
            glDrawArraysInstanced(GL_TRIANGLES, 0, points_emptySphere->size()/3, pos_emptySphere->size()/3);
        else if(pos_emptySphere->size()>0)
        {
            rendering_program_spheres.setAttributeValue("center", QVector4D(pos_emptySphere->at(0), pos_emptySphere->at(1), pos_emptySphere->at(2), 1.0));
            glDrawArrays(GL_TRIANGLES, 0, points_emptySphere->size()/3);
        }
        vao[12].release();
        rendering_program_spheres.release();
    }

}

void Viewer::drawVertex(const Point_3& p, std::vector<float> *vertices)
{

    vertices->push_back(p.x()); vertices->push_back(p.y()); vertices->push_back(p.z());


}

void Viewer::drawEdge(const Point_3& from, const Point_3& to, std::vector<float> *vertices)
{
    vertices->push_back( from.x()); vertices->push_back(from.y()); vertices->push_back(from.z());
    vertices->push_back( to.x()); vertices->push_back(to.y()); vertices->push_back(to.z());
}

void Viewer::drawFacet(const Triangle_3& t, std::vector<float> *vertices)
{
    Point_3 p0 = t.vertex(0);
    Point_3 p1 = t.vertex(1);
    Point_3 p2 = t.vertex(2);
    vertices->push_back( p0.x()); vertices->push_back(p0.y()); vertices->push_back(p0.z());
    vertices->push_back( p1.x()); vertices->push_back(p1.y()); vertices->push_back(p1.z());
    vertices->push_back( p2.x()); vertices->push_back(p2.y()); vertices->push_back(p2.z());

}


/*************************************************************/
/*  Select functions */

void Viewer::drawWithNames()
{
    for(int i=0; i<m_pScene->m_vhArray.size(); ++i) {
        // push a name for each point onto the name stack
        // note: it can NOT be used between glBegin and glEnd
        ::glPushName( i );

        // draw the point
        ::glBegin(GL_POINTS);
        Point_3& p = m_pScene->m_vhArray.at(i)->point();
        ::glVertex3f(p.x(), p.y(), p.z());
        ::glEnd();

        // pop one name off the top of the name stack
        ::glPopName();
    }//end-for-points

    // push a name for the newly inserted point
    if( m_curMode == INSERT_PT && m_hasNewPt ) {
        ::glPushName( ::GLuint(-1) );
        ::glBegin(GL_POINTS);
        ::glVertex3f(m_newPt.x(), m_newPt.y(), m_newPt.z());
        ::glEnd();
        ::glPopName();
    }//end-if-newPt
}

void Viewer::endSelection(const QPoint& /*point*/)
{
    // flush GL buffers
    ::glFlush();

    // reset GL_RENDER mode (was GL_SELECT) and get the number of selected points
    size_t nSel = ::glRenderMode(GL_RENDER);

    /* No selection */
    if( nSel <= 0 ) {
        if( m_curMode == SELECT )
            m_isPress = false;
    }//end-if-notselected

    // each hit record has 4 data: # of names in name stack, min and max depth of old hits,
    //  name stack contents [see glSelectBuffer man page for more details]
    //  i.e. (selectBuffer())[4*i+3] is the id pushed on the stack

    /* Check whether the new point is clicked on */
    else if( m_curMode == INSERT_PT ) {
        if( m_hasNewPt && (selectBuffer())[3] == ::GLuint(-1) )
            m_isMoving = true;
    }//end-if-inspt

    /* Check whether vertex is clicked on */
    else if( m_curMode == MOVE ) {
        m_isMoving = true;
        m_vidMoving = (selectBuffer())[3];
        // compute the corresponding size of trackball, i.e. selectedV is on the ball
        Point_3 p = m_pScene->m_vhArray.at( m_vidMoving )->point();
        m_fRadius = sqrt( p.x()*p.x() + p.y()*p.y() + p.z()*p.z() );
    }//end-if-move

    /* Store current selections */
    else { // m_curMode == SELECT
        if( m_selMode == NORMAL ) {
            // remove the old selections
            for(QList<int>::iterator vit=m_vidSeled.begin();
                vit < m_vidSeled.end(); ++vit) {
                m_pScene->m_vhArray.at(*vit)->setSeled( false );
            }
            m_vidSeled.clear();

            // record the new selections
            for(std::size_t i=0; i<nSel; ++i) {
                m_vidSeled.push_back( (selectBuffer())[4*i+3] );
                m_pScene->m_vhArray.at( m_vidSeled.back() )->setSeled();
            }
        } else {
            for(std::size_t i=0; i<nSel; ++i) {
                if( !m_vidSeled.contains( (selectBuffer())[4*i+3] ) ) {
                    m_vidSeled.push_back( (selectBuffer())[4*i+3] );
                    m_pScene->m_vhArray.at( (selectBuffer())[4*i+3] )->setSeled();
                }//end-if-contain
            }//end-for
        }//end-if-add
    }//end-if-sel
}

/*************************************************************/
/*  Mouse and Keyboard functions */

void Viewer::mousePressEvent(QMouseEvent *event)
{
    // button() holds the button that caused the event
    //  note: for mouse move event, button() always return Qt::NoButton
    // modifiers() holds the keyboard modifier flags at the time of the event
    // buttons() holds the button state when the event was generated,
    //  i.e. all buttons that are pressed down
    // pos() holds the mouse cursor's position relative to the receiving widget

    // Get event modifiers key
    const Qt::KeyboardModifiers modifiers = event->modifiers();

    if( m_curMode == INSERT_V
            && event->button() == Qt::LeftButton && modifiers == Qt::SHIFT ) {
        m_isPress = true;
    }//end-if-insv

    else if(m_curMode == INSERT_PT && event->button() == Qt::LeftButton ) {
        /* shift+left to insert */
        if( modifiers == Qt::SHIFT ) {
            if( m_pScene->m_dt.is_valid() && m_pScene->m_dt.dimension() == 3 )
                m_isPress = true;
            else
                displayMessage( tr("There exists no triangulation yet.") );
            m_hasNewPt = false;
        } else { /* left button to move */
            m_isMoving = false;
            // define selection window (default was 3)
            setSelectRegionWidth( 10 );
            setSelectRegionHeight( 10 );
            // perform the selection
            select( event->pos() );
            if( m_isMoving )
                // redraw window
            {  changed(); updateGL();}
            else
                // if no point is selected, then regular action (rotation) will be performed
                QGLViewer::mousePressEvent(event);
        }//end-if-shift
    }//end-if-inspt

    else if( m_curMode == SELECT && event->button() == Qt::LeftButton ) {
        // set the selection mode
        switch( modifiers ) {
        case Qt::SHIFT : // select
            m_isPress = true;
            m_selMode = NORMAL;
            // initialize multiple selection window
            m_rectSel = QRect( event->pos(), event->pos() );
            // redraw window
            updateGL();
            break;
        case Qt::CTRL : // add selection
            m_isPress = true;
            m_selMode = ADD;
            // initialize multiple selection window
            m_rectSel = QRect( event->pos(), event->pos() );
            // redraw window
            updateGL();
            break;
        default: // rotate
            QGLViewer::mousePressEvent(event);
            break;
        }
    }//end-if-select

    else if(m_curMode == MOVE && event->button() == Qt::LeftButton ) {
        m_isMoving = false;
        // define selection window (default was 3)
        setSelectRegionWidth( 10 );
        setSelectRegionHeight( 10 );
        // perform the selection
        select( event->pos() );
        if( m_isMoving ) // redraw window
        {  changed(); updateGL();}
        else // if no point is selected, then regular action (rotation) will be performed
            QGLViewer::mousePressEvent(event);
    }//end-if-move

    else if( m_curMode == FINDNB
             && event->button() == Qt::LeftButton && modifiers == Qt::SHIFT ) {
        if( m_pScene->m_dt.is_valid() && m_pScene->m_dt.dimension() == 3 )
            m_isPress = true;
        else
            displayMessage( tr("There exists no triangulation yet.") );
    }//end-if-findnb

    else if( m_curMode == EMPTYSPH
             && event->button() == Qt::LeftButton && modifiers == Qt::SHIFT ) {
        if( m_pScene->m_dt.is_valid() && m_pScene->m_dt.dimension() == 3 )
            m_isPress = true;
        else
            displayMessage( tr("There exists no triangulation yet.") );
        m_hasEmptyS = false;
    }//end-if-emptyS

    else
        QGLViewer::mousePressEvent(event);

}

void Viewer::mouseMoveEvent(QMouseEvent *event)
{
    if( m_curMode == INSERT_PT && m_isMoving ) {
        Vec pt;
        if( computeIntersect( event->pos(), pt ) ) {
            m_newPt = Point_3(pt.x, pt.y, pt.z);
            // compute the conflict hole induced by point p
            computeConflict( m_newPt );
            // redraw
            changed(); updateGL();
        }//end-if-compute
    }//end-if-inspt

    else if( m_curMode == SELECT && m_isPress ) {
        // update multiple selection window
        m_rectSel.setBottomRight( event->pos() );
        // redraw
        //changed(); updateGL();
    }//end-if-sel

    else if( m_curMode == MOVE && m_isMoving ) {
        Vec pt;
        if( computeIntersect( event->pos(), pt ) ) {
            // note: QList::operator[] return a modifiable reference;
            //   while QList::at return a const reference
            // move_if_no_collision moves the point stored in v to pt
            //  if there is not already another vertex placed on pt,
            //  the triangulation is modified s.t. the new position of v is pt;
            //  otherwise, the vertex at point pt is returned.
            Vertex_handle vh = m_pScene->m_dt.move_if_no_collision(
                        m_pScene->m_vhArray.at( m_vidMoving ),
                        Point_3( pt.x, pt.y, pt.z ) );
            int id1 = m_pScene->m_vhArray.indexOf( vh );
            int id2 = m_pScene->m_vhArray.indexOf( vh, m_vidMoving+1 );
            // remove the duplicate in vhArray
            if( id1 != m_vidMoving )
                m_pScene->m_vhArray.removeAt( id1 );
            else if( id2 != -1 )
                m_pScene->m_vhArray.removeAt( id2 );
            m_pScene->m_vhArray[m_vidMoving] = vh;
        }//end-if-compute

        // redraw
        //  changed(); updateGL();
    }//end-if-move

    else
        QGLViewer::mouseMoveEvent(event);
}

void Viewer::mouseReleaseEvent(QMouseEvent *event)
{
    /* INS_V mode - Shift+Left: compute and insert a vertex */
    if( m_curMode == INSERT_V && m_isPress ) {
        m_isPress = false;
        Vec pt;
        if( computeIntersect( event->pos(), pt ) ) {
            m_pScene->m_vhArray.push_back( m_pScene->m_dt.insert( Point_3( pt.x, pt.y, pt.z ) ) );
        }//end-if-compute

        // redraw
        changed(); updateGL();
    }//end-if-ins

    /* INS_PT mode - Shift+Left: compute and insert a point */
    else if( m_curMode == INSERT_PT && m_isPress ) {
        m_isPress = false;
        Vec pt;
        if( computeIntersect( event->pos(), pt ) ) {
            m_hasNewPt = true;
            m_newPt = Point_3(pt.x, pt.y, pt.z);
            // compute the conflict hole induced by point p
            computeConflict( m_newPt );
        }//end-if-compute

        // redraw
        changed(); updateGL();
    }//end-if-inspt

    /* INS_PT mode - Left: compute and insert a point */
    else if( m_curMode == INSERT_PT && m_isMoving ) {
        m_isMoving = false;
        Vec pt;
        if( computeIntersect( event->pos(), pt ) ) {
            m_newPt = Point_3(pt.x, pt.y, pt.z);
            // compute the conflict hole induced by point p
            computeConflict( m_newPt );
        }//end-if-compute

        // redraw
        changed(); updateGL();
    }//end-if-inspt

    /* SEL mode - Left: terminate multiple point selection */
    else if( m_curMode == SELECT && m_isPress ) {
        // might swap left/right and top/bottom to make rectanle valid
        m_rectSel = m_rectSel.normalized();

        if( m_rectSel.width() == 1 && m_rectSel.height() == 1 ) { /* select a point */
            // set a default selection window
            setSelectRegionWidth( 10 );
            setSelectRegionHeight( 10 );
            // compute rectangle center and perform selection
            select( m_rectSel.center() );
            if( m_isPress ) {
                m_isPress = false;
            } else {
                displayMessage( tr("No point is selected.") );
            }
        } else {  /* select multiple points, ie. selection window > 1 */
            // define selection window
            if( m_rectSel.width() < 10 )
                setSelectRegionWidth( 10 );
            else
                setSelectRegionWidth( m_rectSel.width() );
            if( m_rectSel.height() < 10 )
                setSelectRegionHeight( 10 );
            else
                setSelectRegionHeight( m_rectSel.height() );
            // compute rectangle center and perform selection
            select( m_rectSel.center() );
            if( m_isPress ) {
                m_isPress = false;
                displayMessage( QString::number(m_vidSeled.size()) + tr(" points are selected") );
            } else { // empty window will cancel the current selection
                for(QList<int>::iterator iit = m_vidSeled.begin(); iit < m_vidSeled.end(); ++iit)
                    m_pScene->m_vhArray.at(*iit)->setSeled( false );
                m_vidSeled.clear();
            }
        }//end-if-selwindow

        // update display to show
        changed(); updateGL();
    }//end-if-select

    /* MOVE mode - Left: terminate point moving */
    else if( m_curMode == MOVE && m_isMoving ) {
        Vec pt;
        if( computeIntersect( event->pos(), pt ) ) {
            // note: QList::operator[] return a modifiable reference;
            //   while QList::at return a const reference
            // move_if_no_collision moves the point stored in v to pt
            //  if there is not already another vertex placed on pt,
            //  the triangulation is modified s.t. the new position of v is pt;
            //  otherwise, the vertex at point pt is returned.
            Vertex_handle vh = m_pScene->m_dt.move_if_no_collision(
                        m_pScene->m_vhArray.at( m_vidMoving ),
                        Point_3( pt.x, pt.y, pt.z ) );
            int id1 = m_pScene->m_vhArray.indexOf( vh );
            int id2 = m_pScene->m_vhArray.indexOf( vh, m_vidMoving+1 );
            // remove the duplicate in vhArray
            if( id1 != m_vidMoving )
                m_pScene->m_vhArray.removeAt( id1 );
            else if( id2 != -1 )
                m_pScene->m_vhArray.removeAt( id2 );
            m_pScene->m_vhArray[m_vidMoving] = vh;
        }//end-if-compute

        // redraw
        changed(); updateGL();
    }//end-if-move

    /* FindNb mode - Shift+Left: find the nearest neighbor of the point */
    else if( m_curMode == FINDNB && m_isPress ) {
        m_isPress = false;
        Vec pt;
        if( computeIntersect( event->pos(), pt ) ) {
            m_queryPt = Point_3( pt.x, pt.y, pt.z );
            m_nearestNb = m_pScene->m_dt.nearest_vertex( m_queryPt );
        }//end-if-compute

        // redraw
        changed(); updateGL();
    }//end-if-findnb

    /* EmptySphere mode - Shift+Left: show the empty sphere of the cell */
    else if( m_curMode == EMPTYSPH && m_isPress ) {
        m_isPress = false;
        Vec pt;
        m_hasEmptyS = computeIntersect( event->pos(), pt );
        if( m_hasEmptyS ) {
            m_queryPt = Point_3( pt.x, pt.y, pt.z );
            // find the cell that contains point p in its interior
            m_cellContain = m_pScene->m_dt.locate( m_queryPt );
            // show error if point is outside the convex hull
            if( m_pScene->m_dt.is_infinite( m_cellContain ) ) {
                m_hasEmptyS = false;
                displayMessage( tr("Query point is outside the convex hull!") );
            } else { /* compute the empty sphere */
                // find the circumcenter of the four vertices of c
                m_centerPt = m_pScene->m_dt.dual( m_cellContain );
                // compute the radius of the empty sphere
                m_fREmptyS = sqrt( CGAL::squared_distance( m_centerPt,
                                                           m_cellContain->vertex(0)->point() ) );
            }
        }//end-if-compute
        // redraw
        changed(); updateGL();
    }//end-if-emptysphere

    else
        QGLViewer::mouseReleaseEvent(event);

}

void Viewer::wheelEvent(QWheelEvent *event)
{
    // Get event modifiers key
    const Qt::KeyboardModifiers modifiers = event->modifiers();

    if( (m_curMode == INSERT_V || m_curMode == FINDNB || m_curMode == EMPTYSPH )
            && modifiers == Qt::SHIFT ) {
        // delta() returns the distance that the wheel is rotated, in eighths of a degree.
        //  note: most mouse types work in steps of 15 degrees
        //  positive value: rotate forwards away from the user;
        //  negative value: rotate backwards toward the user.
        m_fRadius += (event->delta()*1.f / m_iStep ); // inc-/decrease by 0.1 per step
        if( m_fRadius < 0.1f )
            m_fRadius = 0.1f;

        // redraw
        changed(); updateGL();
    }//end-if-insv

    else if( m_curMode == INSERT_PT && modifiers == Qt::SHIFT ) {
        // delta() returns the distance that the wheel is rotated, in eighths of a degree.
        //  note: most mouse types work in steps of 15 degrees
        //  positive value: rotate forwards away from the user;
        //  negative value: rotate backwards toward the user.
        float origR = m_fRadius;
        m_fRadius += (event->delta()*1.f / m_iStep ); // inc-/decrease by 0.1 per step
        if( m_fRadius < 0.1f )
            m_fRadius = 0.1f;
        // update the new point and its conflict region
        if( m_hasNewPt ) {
            origR = m_fRadius / origR;
            m_newPt = Point_3( m_newPt.x()*origR, m_newPt.y()*origR, m_newPt.z()*origR );
            // compute the conflict hole induced by point p
            computeConflict( m_newPt );
        }//end-if-conflict

        // redraw
        changed(); updateGL();
    }//end-if-inspt

    // resize the trackball when moving a point
    else if( m_curMode == MOVE && modifiers == Qt::SHIFT && m_isMoving ) {
        float origR = m_fRadius;
        m_fRadius += (event->delta()*1.f / m_iStep ); // inc-/decrease by 0.1 per step
        if( m_fRadius < 0.1f )
            m_fRadius = 0.1f;
        origR = m_fRadius / origR;
        Point_3 pt = m_pScene->m_vhArray.at( m_vidMoving )->point();
        // note: QList::operator[] return a modifiable reference;
        //   while QList::at return a const reference
        // move_if_no_collision moves the point stored in v to pt
        //  if there is not already another vertex placed on pt,
        //  the triangulation is modified s.t. the new position of v is pt;
        //  otherwise, the vertex at point pt is returned.
        Vertex_handle vh = m_pScene->m_dt.move_if_no_collision(
                    m_pScene->m_vhArray.at( m_vidMoving ),
                    Point_3( pt.x()*origR, pt.y()*origR, pt.z()*origR ) );
        int id1 = m_pScene->m_vhArray.indexOf( vh );
        int id2 = m_pScene->m_vhArray.indexOf( vh, m_vidMoving+1 );
        // remove the duplicate in vhArray
        if( id1 != m_vidMoving )
            m_pScene->m_vhArray.removeAt( id1 );
        else if( id2 != -1 )
            m_pScene->m_vhArray.removeAt( id2 );
        m_pScene->m_vhArray[m_vidMoving] = vh;

        // redraw
        changed(); updateGL();
    }//end-if-move

    else
        QGLViewer::wheelEvent(event);
}

void Viewer::keyPressEvent(QKeyEvent *event)
{
    // Get event modifiers key
    const Qt::KeyboardModifiers modifiers = event->modifiers();

    /* Insert the newly inserted point as a vertex */
    if( m_curMode == INSERT_PT && m_hasNewPt
            && ( event->key()==Qt::Key_Return || event->key()==Qt::Key_Enter )
            && modifiers==Qt::NoButton ) {
        Facet& f = m_boundaryFacets.first(); // a boundary facet, i.e. a pair (cell_handle, i)
        // insert_in_hole will create a new vertex by starring a hole
        //   i.e. delete all conflict cells, create a new vertex,
        //        and for each boundary facet, create a new cell with the new vertex
        // it takes in an iterator range of conflict cells which specifies a hole
        //   and (begin, i) is a boundary facet that begin is one of the conflict cell
        //   but begin->neighbor(i) is not
        // it returns the handle of the new vertex
        m_pScene->m_vhArray.push_back( m_pScene->m_dt.insert_in_hole( m_newPt, // the point
                                                                      m_conflictCells.begin(), // cell_begin
                                                                      m_conflictCells.end(), // cell_end
                                                                      f.first, // cell_handle begin
                                                                      f.second ) ); // integer i

        m_hasNewPt = false;
        // erase old conflict hole info
        m_boundaryFacets.clear();
        m_conflictCells.clear();

        // redraw
        changed(); updateGL();
    }//end-if-insVertex

    /* Cancel the newly inserted point and its conflict region */
    else if( m_curMode == INSERT_PT && m_hasNewPt
             && event->key()==Qt::Key_Escape && modifiers==Qt::NoButton ) {
        m_hasNewPt = false;
        // erase old conflict hole info
        m_boundaryFacets.clear();
        m_conflictCells.clear();

        // redraw
        changed(); updateGL();
    }//end-if-escapeIns

    /* Delete selected points */
    else if( m_curMode == SELECT
             && event->key()==Qt::Key_Delete && modifiers==Qt::NoButton ) {
        // sort selected id's in descending order
        qSort(m_vidSeled.begin(), m_vidSeled.end(), qGreater<int>());
        for(QList<int>::iterator vit=m_vidSeled.begin(); vit<m_vidSeled.end(); ++vit) {
            // remove the selected point from DT and vertex_handle_array
            // note: QList::takeAt will removes the item at index position i and returns it.
            m_pScene->m_dt.remove( m_pScene->m_vhArray.takeAt( *vit ) );
        }
        // clear the selection buffer
        m_vidSeled.clear();

        // redraw
        changed(); updateGL();
    }//end-if-del

    /* Cancel the selection */
    else if( m_curMode == SELECT
             && event->key()==Qt::Key_Escape && modifiers==Qt::NoButton ) {
        // clear the selection buffer
        for(QList<int>::iterator iit=m_vidSeled.begin(); iit<m_vidSeled.end(); ++iit) {
            m_pScene->m_vhArray.at(*iit)->setSeled( false );
        }
        m_vidSeled.clear();

        // redraw
        changed(); updateGL();
    }//end-if-escapeSel


    /* Show/hide the trackball when drawing the empty sphere */
    else if( m_curMode == EMPTYSPH
             && event->key()==Qt::Key_S && modifiers==Qt::NoButton ) {
        m_showTrackball = !m_showTrackball;
        // redraw
        changed(); updateGL();
    }//end-if-showBall

    else
        QGLViewer::keyPressEvent(event);


    // redraw
    changed(); updateGL();
}

/*************************************************************/
/*  Computation functions */

bool Viewer::computeIntersect( const QPoint & pos, Vec & pt )
{
    Vec eye, dir;
    // Compute eye position and direction to the clicked point,
    //  used to draw a representation of the intersecting line
    camera()->convertClickToLine( pos, eye, dir );
    // Compute the intersection point with the sphere
    //  note that the center of the sphere is at the origin (0, 0, 0)
    //  thus, (1) pt = eye + t*dir and (2) dist( pt, origin ) = radius
    //  i.e. (x_eye + t*x_dir)^2 + (y_eye + t*y_dir)^2 + (z_eye + t*z_dir)^2 = r^2
    //  --> t^2( dir*dir ) + 2t( eye*dir ) + eye*eye - r^2 = 0
    //      where "dir*dir" is the dot product of vector dir
    //  we need to solve t and the smaller t (nearer to eye position) is what we want
    float a = dir*dir;
    float b = eye*dir;
    float c = eye*eye - m_fRadius*m_fRadius;
    float delta = b*b - a*c;
    if( delta < 0 ) {
        displayMessage( tr("Point is not on the sphere!") );
        return false;
    } else {
        float t = ( (-1.)*b - sqrt(delta) ) / a;
        pt = eye + t*dir;
        return true;
    }
}

void Viewer::computeConflict( Point_3 pt )
{
    // find the cell that contains point p in its interior
    m_cellContain = m_pScene->m_dt.locate( pt );
    // erase old conflict hole info
    m_boundaryFacets.clear();
    m_conflictCells.clear();
    // show msg if point is outside the convex hull
    if( m_pScene->m_dt.is_infinite( m_cellContain ) )
        displayMessage( tr("Note: point is outside the convex hull.") );
    // compute the conflict hole induced by point p
    m_pScene->m_dt.find_conflicts( pt, // the point
                                   m_cellContain, // starting cell that must be in conflict
                                   std::back_inserter(m_boundaryFacets), // the facets on the boundary
                                   std::back_inserter(m_conflictCells) ); // the cells in conflict
}

/*************************************************************/
/*  Animation functions */

void Viewer::toggleIncremental(bool on) {
    if( on ) {  // play
        if( m_incrementalPts.isEmpty() ) {
            /* start play */
            if( m_pScene->m_dt.number_of_vertices() == 0 ) {
                CGAL::Random_points_in_cube_3<Point_3> pts_generator(1.0);
                CGAL::cpp11::copy_n( pts_generator, 100, std::back_inserter(m_incrementalPts) );
            } else {
                for(QList<Vertex_handle>::iterator vit = m_pScene->m_vhArray.begin();
                    vit < m_pScene->m_vhArray.end(); ++vit) {
                    m_incrementalPts.push_back( (*vit)->point() );
                }//end-for
                // erase existing vertices
                initClean();
            }//end-if-pts
            // sorts points in a way that improves space locality
            CGAL::spatial_sort( m_incrementalPts.begin(), m_incrementalPts.end() );
            // set the current to "hightlight the new point"
            m_curStep = INIT;
        }/* else resume play */

        // set up the timer
        m_pTimer->start(1000);
    } else { // pause
        m_pTimer->stop();
    }

    // redraw
    changed(); updateGL();
}

void Viewer::stopIncremental() {
    if( !m_incrementalPts.isEmpty() ) {
        // will call toggleIncremental to stop the timer
    Q_EMIT( stopIncAnimation() );

        // insert the rest points
        for(QList<Point_3>::iterator pit=m_incrementalPts.begin();
            pit < m_incrementalPts.end(); ++pit) {
            Vertex_handle hint;
            if( m_pScene->m_vhArray.isEmpty() ) {
                hint = m_pScene->m_dt.insert( *pit );
            } else {
                hint = m_pScene->m_vhArray.last();
                hint = m_pScene->m_dt.insert( *pit, hint );
            }
            m_pScene->m_vhArray.push_back( hint );
        }
        m_incrementalPts.clear();
    }

    // redraw
    changed(); updateGL();
}

void Viewer::incremental_insert() {
    Vertex_handle hint;
    if( !m_incrementalPts.isEmpty() ) {
        switch( m_curStep ) {
        case INIT:  // end of INIT: get the next-to-insert point
            m_curIncPt = m_incrementalPts.at(0);
            m_curStep = NEWPT;
            break;
        case NEWPT:  // end of NEWPT: locate the cell containing the point
            if( m_pScene->m_dt.is_valid() && m_pScene->m_dt.dimension() == 3 ) {
                computeConflict( m_curIncPt );
                m_curStep = CELL;
            }
            else {
                // popup the first point and insert it
                m_curIncPt = m_incrementalPts.takeFirst();
                if( m_pScene->m_vhArray.isEmpty() ) {
                    hint = m_pScene->m_dt.insert( m_curIncPt );
                }
                else {
                    hint = m_pScene->m_vhArray.last();
                    hint = m_pScene->m_dt.insert( m_curIncPt, hint );
                }
                m_pScene->m_vhArray.push_back( hint );
                m_curStep = INIT;
            }
            break;
        case CELL:  // end of CELL: compute the conflict region
            m_curStep = CONFLICT;
            break;
        case CONFLICT:  // end of CONFLICT: do the insertion and go back to INIT
            // popup the first point and insert it
            m_curIncPt = m_incrementalPts.takeFirst();
            if( m_pScene->m_vhArray.isEmpty() ) {
                hint = m_pScene->m_dt.insert( m_curIncPt );
            } else {
                hint = m_pScene->m_vhArray.last();
                hint = m_pScene->m_dt.insert( m_curIncPt, hint );
            }
            m_pScene->m_vhArray.push_back( hint );
            m_curStep = INIT;
            break;
        }//end-of-switch
    } else {
        /* if finished, then start over */
        for(QList<Vertex_handle>::iterator vit = m_pScene->m_vhArray.begin();
            vit < m_pScene->m_vhArray.end(); ++vit) {
            m_incrementalPts.push_back( (*vit)->point() );
        }//end-for
        // erase existing vertices
        initClean();
        // set the current to "hightlight the new point"
        m_curStep = INIT;
    }

    // redraw
    changed(); updateGL();
}


void Viewer::draw_cylinder(float R, int prec, std::vector<float> *vertices, std::vector<float> *normals)
{
    vertices->resize(0);
    normals->resize(0);
    // int rings=360/prec, sectors=360/prec;
    // float T, P;
    // float x[4],y[4],z[4];
    //Closing nicely the tubes will cause z-fighting and the spherical parts will get all messy
    /*
//top of the cylinder
    for(int t=0; t<360; t+=sectors)
    {

        vertices->push_back(0);
        vertices->push_back(R+1);
        vertices->push_back(0);


        normals->push_back(0);
        normals->push_back(1);
        normals->push_back(0);



        P = rings*M_PI/180.0;
        T = t*M_PI/180.0;
        x[1] = sin(P) * cos(T) ;
        z[1] = sin(P) * sin(T) ;
        y[1] = cos(P);
        vertices->push_back(R * x[1]);
        vertices->push_back(R * y[1]+1);
        vertices->push_back(R * z[1]);


        normals->push_back(x[1]);
        normals->push_back(y[1]);
        normals->push_back(z[1]);

        //
        P = rings*M_PI/180.0;
        T = (t+sectors)*M_PI/180.0;
        x[2] = sin(P) * cos(T) ;
        z[2] = sin(P) * sin(T) ;
        y[2] = cos(P);
        vertices->push_back(R * x[2]);
        vertices->push_back(R * y[2]+1);
        vertices->push_back(R * z[2]);

        normals->push_back(x[2]);
        normals->push_back(y[2]);
        normals->push_back(z[2]);

    }
    //Body of the sphere
    for (int p=rings; p<90; p+=rings)
        for(int t=0; t<360; t+=sectors)
        {
            //A
            P = p*M_PI/180.0;
            T = t*M_PI/180.0;
            x[0] = sin(P) * cos(T) ;
            z[0] = sin(P) * sin(T) ;
            y[0] = cos(P);

            vertices->push_back(R * x[0]);
            vertices->push_back(R * y[0]+1);
            vertices->push_back(R * z[0]);


            normals->push_back(x[0]);
            normals->push_back(y[0]);
            normals->push_back(z[0]);

            //B
            P = (p+rings)*M_PI/180.0;
            T = t*M_PI/180.0;
            x[1] = sin(P) * cos(T) ;
            z[1] = sin(P) * sin(T) ;
            y[1] = cos(P);
            vertices->push_back(R * x[1]);
            vertices->push_back(R * y[1]+1);
            vertices->push_back(R * z[1]);


            normals->push_back(x[1]);
            normals->push_back(y[1]);
            normals->push_back(z[1]);

            //C
            P = p*M_PI/180.0;
            T = (t+sectors)*M_PI/180.0;
            x[2] = sin(P) * cos(T) ;
            z[2] = sin(P) * sin(T) ;
            y[2] = cos(P);
            vertices->push_back(R * x[2]);
            vertices->push_back(R * y[2]+1);
            vertices->push_back(R * z[2]);


            normals->push_back(x[2]);
            normals->push_back(y[2]);
            normals->push_back(z[2]);
            //D
            P = (p+rings)*M_PI/180.0;
            T = (t+sectors)*M_PI/180.0;
            x[3] = sin(P) * cos(T) ;
            z[3] = sin(P) * sin(T) ;
            y[3] = cos(P);
            vertices->push_back(R * x[3]);
            vertices->push_back(R * y[3]+1);
            vertices->push_back(R * z[3]);


            normals->push_back(x[3]);
            normals->push_back(y[3]);
            normals->push_back(z[3]);



            vertices->push_back(R * x[1]);
            vertices->push_back(R * y[1]+1);
            vertices->push_back(R * z[1]);


            normals->push_back(x[1]);
            normals->push_back(y[1]);
            normals->push_back(z[1]);

            vertices->push_back(R * x[2]);
            vertices->push_back(R * y[2]+1);
            vertices->push_back(R * z[2]);


            normals->push_back(x[2]);
            normals->push_back(y[2]);
            normals->push_back(z[2]);

        }

*/
    //body of the cylinder
    for(int d = 0; d<360; d+= 360/prec)
    {

        //point A1
        float D = d*M_PI/180.0;
        vertices->push_back(R * sin(D));
        vertices->push_back(0);
        vertices->push_back(R * cos(D));

        normals->push_back(sin(D));
        normals->push_back(0);
        normals->push_back(cos(D));

        //point B1
        vertices->push_back(R * sin(D));
        vertices->push_back(1);
        vertices->push_back(R * cos(D));

        normals->push_back(sin(D));
        normals->push_back(0);
        normals->push_back(cos(D));

        //point C1
        D = (d+360/prec)*M_PI/180.0;
        vertices->push_back(R * sin(D));
        vertices->push_back(1);
        vertices->push_back(R * cos(D));

        normals->push_back(sin(D));
        normals->push_back(0);
        normals->push_back(cos(D));

        //point A2
        D = (d+360/prec)*M_PI/180.0;
        vertices->push_back(R * sin(D));
        vertices->push_back(1);
        vertices->push_back(R * cos(D));

        normals->push_back(sin(D));
        normals->push_back(0);
        normals->push_back(cos(D));

        //point B2
        vertices->push_back(R * sin(D));
        vertices->push_back(0);
        vertices->push_back(R * cos(D));

        normals->push_back(sin(D));
        normals->push_back(0);
        normals->push_back(cos(D));

        //point C2
        D = d*M_PI/180.0;
        vertices->push_back(R * sin(D));
        vertices->push_back(0);
        vertices->push_back(R * cos(D));

        normals->push_back(sin(D));
        normals->push_back(0);
        normals->push_back(cos(D));

    }
    /*
   //bottom of the cylinder
       for(int t=0; t<360; t+=sectors)
       {

           vertices->push_back(0);
           vertices->push_back(-R);
           vertices->push_back(0);


           normals->push_back(0);
           normals->push_back(-1);
           normals->push_back(0);



           P = rings*M_PI/180.0;
           T = t*M_PI/180.0;
           x[1] = sin(P) * cos(T) ;
           z[1] = sin(P) * sin(T) ;
           y[1] = cos(P);
           vertices->push_back(R * x[1]);
           vertices->push_back(R * y[1]);
           vertices->push_back(R * z[1]);


           normals->push_back(x[1]);
           normals->push_back(y[1]);
           normals->push_back(z[1]);

           //
           P = rings*M_PI/180.0;
           T = (t+sectors)*M_PI/180.0;
           x[2] = sin(P) * cos(T) ;
           z[2] = sin(P) * sin(T) ;
           y[2] = cos(P);
           vertices->push_back(R * x[2]);
           vertices->push_back(R * y[2]);
           vertices->push_back(R * z[2]);

           normals->push_back(x[2]);
           normals->push_back(y[2]);
           normals->push_back(z[2]);

       }
       //Body of the sphere
       for (int p=90; p<180; p+=rings)
           for(int t=0; t<360; t+=sectors)
           {
               //A
               P = p*M_PI/180.0;
               T = t*M_PI/180.0;
               x[0] = sin(P) * cos(T) ;
               z[0] = sin(P) * sin(T) ;
               y[0] = cos(P);

               vertices->push_back(R * x[0]);
               vertices->push_back(R * y[0]);
               vertices->push_back(R * z[0]);


               normals->push_back(x[0]);
               normals->push_back(y[0]);
               normals->push_back(z[0]);

               //B
               P = (p+rings)*M_PI/180.0;
               T = t*M_PI/180.0;
               x[1] = sin(P) * cos(T) ;
               z[1] = sin(P) * sin(T) ;
               y[1] = cos(P);
               vertices->push_back(R * x[1]);
               vertices->push_back(R * y[1]);
               vertices->push_back(R * z[1]);


               normals->push_back(x[1]);
               normals->push_back(y[1]);
               normals->push_back(z[1]);

               //C
               P = p*M_PI/180.0;
               T = (t+sectors)*M_PI/180.0;
               x[2] = sin(P) * cos(T) ;
               z[2] = sin(P) * sin(T) ;
               y[2] = cos(P);
               vertices->push_back(R * x[2]);
               vertices->push_back(R * y[2]);
               vertices->push_back(R * z[2]);


               normals->push_back(x[2]);
               normals->push_back(y[2]);
               normals->push_back(z[2]);
               //D
               P = (p+rings)*M_PI/180.0;
               T = (t+sectors)*M_PI/180.0;
               x[3] = sin(P) * cos(T) ;
               z[3] = sin(P) * sin(T) ;
               y[3] = cos(P);
               vertices->push_back(R * x[3]);
               vertices->push_back(R * y[3]);
               vertices->push_back(R * z[3]);


               normals->push_back(x[3]);
               normals->push_back(y[3]);
               normals->push_back(z[3]);



               vertices->push_back(R * x[1]);
               vertices->push_back(R * y[1]);
               vertices->push_back(R * z[1]);


               normals->push_back(x[1]);
               normals->push_back(y[1]);
               normals->push_back(z[1]);

               vertices->push_back(R * x[2]);
               vertices->push_back(R * y[2]);
               vertices->push_back(R * z[2]);


               normals->push_back(x[2]);
               normals->push_back(y[2]);
               normals->push_back(z[2]);

           }*/


}

void Viewer::draw_sphere(float R, int prec, std::vector<float> *vertices, std::vector<float> *normals)
{
    vertices->resize(0);
    normals->resize(0);
    int rings=180/prec, sectors=360/prec;
    float T, P;
    float x[4],y[4],z[4];


    //Top of the sphere
    for(int t=0; t<360; t+=sectors)
    {

        vertices->push_back(0);
        vertices->push_back(0);
        vertices->push_back(R);


        normals->push_back(0);
        normals->push_back(0);
        normals->push_back(1);



        P = rings*M_PI/180.0;
        T = t*M_PI/180.0;
        x[1] = sin(P) * cos(T) ;
        y[1] = sin(P) * sin(T) ;
        z[1] = cos(P);
        vertices->push_back(R * x[1]);
        vertices->push_back(R * y[1]);
        vertices->push_back(R * z[1]);


        normals->push_back(x[1]);
        normals->push_back(y[1]);
        normals->push_back(z[1]);

        //
        P = rings*M_PI/180.0;
        T = (t+sectors)*M_PI/180.0;
        x[2] = sin(P) * cos(T) ;
        y[2] = sin(P) * sin(T) ;
        z[2] = cos(P);
        vertices->push_back(R * x[2]);
        vertices->push_back(R * y[2]);
        vertices->push_back(R * z[2]);

        normals->push_back(x[2]);
        normals->push_back(y[2]);
        normals->push_back(z[2]);

    }

    //Body of the sphere
    for (int p=rings; p<180-rings; p+=rings)
        for(int t=0; t<360; t+=sectors)
        {
            //A
            P = p*M_PI/180.0;
            T = t*M_PI/180.0;
            x[0] = sin(P) * cos(T) ;
            y[0] = sin(P) * sin(T) ;
            z[0] = cos(P);

            vertices->push_back(R * x[0]);
            vertices->push_back(R * y[0]);
            vertices->push_back(R * z[0]);


            normals->push_back(x[0]);
            normals->push_back(y[0]);
            normals->push_back(z[0]);

            //B
            P = (p+rings)*M_PI/180.0;
            T = t*M_PI/180.0;
            x[1] = sin(P) * cos(T) ;
            y[1] = sin(P) * sin(T) ;
            z[1] = cos(P);
            vertices->push_back(R * x[1]);
            vertices->push_back(R * y[1]);
            vertices->push_back(R * z[1]);


            normals->push_back(x[1]);
            normals->push_back(y[1]);
            normals->push_back(z[1]);

            //C
            P = p*M_PI/180.0;
            T = (t+sectors)*M_PI/180.0;
            x[2] = sin(P) * cos(T) ;
            y[2] = sin(P) * sin(T) ;
            z[2] = cos(P);
            vertices->push_back(R * x[2]);
            vertices->push_back(R * y[2]);
            vertices->push_back(R * z[2]);


            normals->push_back(x[2]);
            normals->push_back(y[2]);
            normals->push_back(z[2]);
            //D
            P = (p+rings)*M_PI/180.0;
            T = (t+sectors)*M_PI/180.0;
            x[3] = sin(P) * cos(T) ;
            y[3] = sin(P) * sin(T) ;
            z[3] = cos(P);
            vertices->push_back(R * x[3]);
            vertices->push_back(R * y[3]);
            vertices->push_back(R * z[3]);


            normals->push_back(x[3]);
            normals->push_back(y[3]);
            normals->push_back(z[3]);



            vertices->push_back(R * x[1]);
            vertices->push_back(R * y[1]);
            vertices->push_back(R * z[1]);


            normals->push_back(x[1]);
            normals->push_back(y[1]);
            normals->push_back(z[1]);

            vertices->push_back(R * x[2]);
            vertices->push_back(R * y[2]);
            vertices->push_back(R * z[2]);


            normals->push_back(x[2]);
            normals->push_back(y[2]);
            normals->push_back(z[2]);

        }
    //Bottom of the sphere
    for(int t=0; t<360; t+=sectors)
    {


        vertices->push_back(0);
        vertices->push_back(0);
        vertices->push_back(-R);


        normals->push_back(0);
        normals->push_back(0);
        normals->push_back(-1);


        P = (180-rings)*M_PI/180.0;
        T = t*M_PI/180.0;
        x[1] = sin(P) * cos(T) ;
        y[1] = sin(P) * sin(T) ;
        z[1] = cos(P);
        vertices->push_back(R * x[1]);
        vertices->push_back(R * y[1]);
        vertices->push_back(R * z[1]);


        normals->push_back(x[1]);
        normals->push_back(y[1]);
        normals->push_back(z[1]);


        P = (180-rings)*M_PI/180.0;
        T = (t+sectors)*M_PI/180.0;
        x[2] = sin(P) * cos(T) ;
        y[2] = sin(P) * sin(T) ;
        z[2] = cos(P);
        vertices->push_back(R * x[2]);
        vertices->push_back(R * y[2]);
        vertices->push_back(R * z[2]);


        normals->push_back(x[2]);
        normals->push_back(y[2]);
        normals->push_back(z[2]);

    }

}
