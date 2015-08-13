#include "Scene.h"


void Scene::compile_shaders()
{

    if(! buffers[0].create() || !buffers[1].create() || !buffers[2].create() || !buffers[3].create()
            || !buffers[4].create() || !buffers[5].create() || !buffers[6].create()
            || !buffers[7].create() || !buffers[8].create() || !buffers[9].create()
            || !buffers[10].create() || !buffers[11].create() || !buffers[12].create()
            || !buffers[13].create() || !buffers[14].create() || !buffers[15].create()
            || !buffers[16].create() || !buffers[17].create() || !buffers[18].create()
            || !buffers[19].create() || !buffers[20].create()
            || !buffers[21].create() || !buffers[22].create() || !buffers[23].create())
    {
        std::cerr<<"VBO Creation FAILED"<<std::endl;
    }

    if(!vao[0].create() || !vao[1].create() || !vao[2].create() || !vao[3].create() ||
            !vao[4].create() || !vao[5].create() ||!vao[6].create() || !vao[7].create()
            || !vao[8].create() || !vao[9].create() || !vao[10].create() || !vao[11].create())
    {
        std::cerr<<"VAO Creation FAILED"<<std::endl;
    }
    //prepare the sphere model that will be instanced rendered in draw.
    draw_sphere(0.02f,15);
    //prepare the tube model that will be instanced rendered when drawing the cube.
    draw_cylinder(0.01f,25, points_cube, normals_cylinder);

    //prepare the tube model that will be instanced rendered when drawing the triangulation.
    draw_cylinder(0.005f,25, points_cylinder, normals_cylinder);


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
        "   vec4 diffuse = abs(dot(N,L)) * light_diff; \n"
        "   vec4 specular = pow(max(dot(R,V), 0.0), spec_power) * light_spec; \n"

        "gl_FragColor = light_amb + diffuse + specular; \n"
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
        "   gl_Position =  mvp_matrix * transfo * vertex; \n"
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

void Scene::compute_elements()
{
    pos_points.resize(0);
    transfo1_cylinder.resize(0);
    transfo2_cylinder.resize(0);
    transfo3_cylinder.resize(0);
    transfo4_cylinder.resize(0);
    transfo1_cube.resize(0);
    transfo2_cube.resize(0);
    transfo3_cube.resize(0);
    transfo4_cube.resize(0);
    transfo1_square.resize(0);
    transfo2_square.resize(0);
    transfo3_square.resize(0);
    transfo4_square.resize(0);

    // Draw vertices
    {
        for (Point_iterator pit = p3dt.periodic_points_begin(it_type) ;
             pit != p3dt.periodic_points_end(it_type) ; pit++)
        {
            if (pit->first.z()!=0 && in_plane)
            {
                continue;
            }
            pos_points.push_back(pit->first.x());pos_points.push_back(pit->first.y());pos_points.push_back(pit->first.z());
        }

    }

    // Draw segments
    {
        Segment_set segments_to_draw;
        primitives_from_geom_it(segments_to_draw);
        if (cube_clipping && !two_color_clipping) segment_clipping(segments_to_draw);
        pos_tube.resize(0);

        for (Segment_set::iterator it = segments_to_draw.begin() ;
             it != segments_to_draw.end(); it++)
        {
            Point p1 = it->source(), p2 = it->target();
            if (in_plane && (p1.z()!=0. || p2.z()!=0.)) continue;
            pos_tube.push_back(p1.x()); pos_tube.push_back(p1.y()); pos_tube.push_back(p1.z());
            pos_tube.push_back(p2.x()); pos_tube.push_back(p2.y()); pos_tube.push_back(p2.z());

            Point p = p1;
            Vector v = p2-p1;

            FT len = (FT)std::sqrt(CGAL_NTS to_double(v*v));

            // normalize
            v = v / len;
            double angle = 0.0;
            if(std::sqrt(CGAL_NTS to_double(v.x()*v.x()+v.y()*v.y())) > 1)
                angle = 90.0f;
            else
                angle =acos(v.y()/std::sqrt(v.x()*v.x()+v.y()*v.y()+v.z()*v.z()))*180.0/M_PI;//asin(std::sqrt(CGAL_NTS to_double(v.x()*v.x()+v.y()*v.y())))/M_PI*180.0;

            Vector axis;
            axis = Vector(v.z(), 0, -v.x());


            QMatrix4x4 matrix;
            matrix.setToIdentity();

            matrix.translate(CGAL_NTS to_double(p.x()),
                             CGAL_NTS to_double(p.y()),
                             CGAL_NTS to_double(p.z()));

            matrix.rotate(angle,CGAL_NTS to_double(axis.x()),
                          CGAL_NTS to_double(axis.y()),
                          CGAL_NTS to_double(axis.z()));
            matrix.scale(1,CGAL_NTS to_double(len),1);

            for(int i=0; i<4; i++)
                transfo1_cylinder.push_back((float)matrix.data()[i]);
            for(int i=4; i<8; i++)
                transfo2_cylinder.push_back((float)matrix.data()[i]);
            for(int i=8; i<12; i++)
                transfo3_cylinder.push_back((float)matrix.data()[i]);
            for(int i=12; i<16; i++)
                transfo4_cylinder.push_back((float)matrix.data()[i]);

        }
    }
    // Draw cube
    {

        QMatrix4x4 matrix;
        matrix.setToIdentity();

        for (float x=0.0; x<2.0; x+=1.0) {
            for (float y=0.0; y<2.0; y+=1.0) {
                matrix.translate(x,0,y);
                //   matrix = matrix.transposed();
                for(int i=0; i<4; i++)
                    transfo1_cube.push_back((float)matrix.data()[i]);
                for(int i=4; i<8; i++)
                    transfo2_cube.push_back((float)matrix.data()[i]);
                for(int i=8; i<12; i++)
                    transfo3_cube.push_back((float)matrix.data()[i]);
                for(int i=12; i<16; i++)
                    transfo4_cube.push_back((float)matrix.data()[i]);
                matrix.translate(-x,0,-y);
            }
        }

        for (float x=0.0; x<2.0; x+=1.0) {
            for (float y=0.0; y<2.0; y+=1.0) {
                matrix.translate(x,y,0);
                matrix.rotate(90,1.0,0.0,0.0);
                for(int i=0; i<4; i++)
                    transfo1_cube.push_back((float)matrix.data()[i]);
                for(int i=4; i<8; i++)
                    transfo2_cube.push_back((float)matrix.data()[i]);
                for(int i=8; i<12; i++)
                    transfo3_cube.push_back((float)matrix.data()[i]);
                for(int i=12; i<16; i++)
                    transfo4_cube.push_back((float)matrix.data()[i]);
                matrix.rotate(90,-1.0,0.0,0.0);
                matrix.translate(-x,-y,0);
            }
        }

        for (float x=0.0; x<2.0; x+=1.0) {
            for (float y=0.0; y<2.0; y+=1.0) {
                matrix.translate(0.0,x,y);
                matrix.rotate(90,0.0,0.0,-1.0);
                for(int i=0; i<4; i++)
                    transfo1_cube.push_back((float)matrix.data()[i]);
                for(int i=4; i<8; i++)
                    transfo2_cube.push_back((float)matrix.data()[i]);
                for(int i=8; i<12; i++)
                    transfo3_cube.push_back((float)matrix.data()[i]);
                for(int i=12; i<16; i++)
                    transfo4_cube.push_back((float)matrix.data()[i]);
                matrix.rotate(90,0.0,0.0,1.0);
                matrix.translate(0.0,-x,-y);
            }
        }


        pos_cube.resize(24*3);

        pos_cube[0]=0.0;  pos_cube[3]=1.0;  pos_cube[6]=0.0; pos_cube[9]= 1.0;
        pos_cube[1]=0.0;  pos_cube[4]=0.0;  pos_cube[7]=1.0; pos_cube[10]=1.0;
        pos_cube[2]=0.0;  pos_cube[5]=0.0;  pos_cube[8]=0.0; pos_cube[11]=0.0;

        pos_cube[12]=0.0;  pos_cube[15]=1.0;  pos_cube[18]=0.0; pos_cube[21]=1.0;
        pos_cube[13]=0.0;  pos_cube[16]=0.0;  pos_cube[19]=1.0; pos_cube[22]=1.0;
        pos_cube[14]=1.0;  pos_cube[17]=1.0;  pos_cube[20]=1.0; pos_cube[23]=1.0;

        pos_cube[24]=0.0;  pos_cube[27]=0.0;  pos_cube[30]=1.0; pos_cube[33]=1.0;
        pos_cube[25]=0.0;  pos_cube[28]=1.0;  pos_cube[31]=0.0; pos_cube[34]=1.0;
        pos_cube[26]=0.0;  pos_cube[29]=0.0;  pos_cube[32]=0.0; pos_cube[35]=0.0;

        pos_cube[36]=0.0;  pos_cube[39]=0.0;  pos_cube[42]=1.0; pos_cube[45]=1.0;
        pos_cube[37]=0.0;  pos_cube[40]=1.0;  pos_cube[43]=0.0; pos_cube[46]=1.0;
        pos_cube[38]=1.0;  pos_cube[41]=1.0;  pos_cube[44]=1.0; pos_cube[47]=1.0;

        pos_cube[48]=0.0;  pos_cube[51]=0.0;  pos_cube[54]=1.0; pos_cube[57]=1.0;
        pos_cube[49]=0.0;  pos_cube[52]=0.0;  pos_cube[55]=0.0; pos_cube[58]=0.0;
        pos_cube[50]=0.0;  pos_cube[53]=1.0;  pos_cube[56]=0.0; pos_cube[59]=1.0;

        pos_cube[60]=0.0;  pos_cube[63]=0.0;  pos_cube[66]=1.0; pos_cube[69]=1.0;
        pos_cube[61]=1.0;  pos_cube[64]=1.0;  pos_cube[67]=1.0; pos_cube[70]=1.0;
        pos_cube[62]=0.0;  pos_cube[65]=1.0;  pos_cube[68]=0.0; pos_cube[71]=1.0;
    }
    //Draw square
    {
        pos_square.resize(24);
        pos_square[0]=0.0; pos_square[3]=1.0; pos_square[6]=0.0; pos_square[9]=1.0;
        pos_square[1]=0.0; pos_square[4]=0.0; pos_square[7]=1.0; pos_square[10]=1.0;
        pos_square[2]=0.0; pos_square[5]=0.0; pos_square[8]=0.0; pos_square[11]=0.0;

        pos_square[12]=0.0; pos_square[15]=0.0; pos_square[18]=1.0; pos_square[21]=1.0;
        pos_square[13]=0.0; pos_square[16]=1.0; pos_square[19]=0.0; pos_square[22]=1.0;
        pos_square[14]=0.0; pos_square[17]=0.0; pos_square[20]=0.0; pos_square[23]=0.0;

        QMatrix4x4 matrix;
        matrix.setToIdentity();

        for (float x=0.0; x<2.0; x+=1.0) {

            matrix.translate(x,0,0);
            for(int i=0; i<4; i++)
                transfo1_square.push_back((float)matrix.data()[i]);
            for(int i=4; i<8; i++)
                transfo2_square.push_back((float)matrix.data()[i]);
            for(int i=8; i<12; i++)
                transfo3_square.push_back((float)matrix.data()[i]);
            for(int i=12; i<16; i++)
                transfo4_square.push_back((float)matrix.data()[i]);
            matrix.translate(-x,0,0);
        }

        for (float y=0.0; y<2.0; y+=1.0) {

            matrix.translate(0,y,0);
            matrix.rotate(90,0,0,-1);
            for(int i=0; i<4; i++)
                transfo1_square.push_back((float)matrix.data()[i]);
            for(int i=4; i<8; i++)
                transfo2_square.push_back((float)matrix.data()[i]);
            for(int i=8; i<12; i++)
                transfo3_square.push_back((float)matrix.data()[i]);
            for(int i=12; i<16; i++)
                transfo4_square.push_back((float)matrix.data()[i]);
            matrix.rotate(90,0,0,1);
            matrix.translate(0,-y,0);
        }



    }
}

void Scene::initialize_buffers()
{
    rendering_program.bind();

    vao[0].bind();
    buffers[0].bind();
    buffers[0].allocate(pos_points.data(), static_cast<int>(pos_points.size()*sizeof(float)));
    poly_vertexLocation[0] = rendering_program.attributeLocation("vertex");
    rendering_program.enableAttributeArray(poly_vertexLocation[0]);
    rendering_program.setAttributeBuffer(poly_vertexLocation[0],GL_FLOAT,0,3);
    buffers[0].release();

    vao[0].release();

    vao[1].bind();
    buffers[1].bind();
    buffers[1].allocate(pos_tube.data(), static_cast<int>(pos_tube.size()*sizeof(float)));
    poly_vertexLocation[0] = rendering_program.attributeLocation("vertex");
    rendering_program.enableAttributeArray(poly_vertexLocation[0]);
    rendering_program.setAttributeBuffer(poly_vertexLocation[0],GL_FLOAT,0,3);
    buffers[1].release();
    vao[1].release();

    vao[11].bind();
    buffers[23].bind();
    buffers[23].allocate(pos_square.data(), static_cast<int>(pos_square.size()*sizeof(float)));
    poly_vertexLocation[0] = rendering_program.attributeLocation("vertex");
    rendering_program.enableAttributeArray(poly_vertexLocation[0]);
    rendering_program.setAttributeBuffer(poly_vertexLocation[0],GL_FLOAT,0,3);
    buffers[23].release();
    vao[11].release();

    vao[2].bind();
    buffers[2].bind();
    buffers[2].allocate(pos_cube.data(), static_cast<int>(pos_cube.size()*sizeof(float)));
    poly_vertexLocation[0] = rendering_program.attributeLocation("vertex");
    rendering_program.enableAttributeArray(poly_vertexLocation[0]);
    rendering_program.setAttributeBuffer(poly_vertexLocation[0],GL_FLOAT,0,3);
    buffers[2].release();
    vao[2].release();

    rendering_program.release();

    rendering_program_spheres.bind();

    vao[3].bind();
    buffers[3].bind();
    buffers[3].allocate(points_spheres.data(), static_cast<int>(points_spheres.size()*sizeof(float)));
    poly_vertexLocation[1] = rendering_program_spheres.attributeLocation("vertex");
    rendering_program_spheres.enableAttributeArray(poly_vertexLocation[1]);
    rendering_program_spheres.setAttributeBuffer(poly_vertexLocation[1],GL_FLOAT,0,3);
    buffers[3].release();

    buffers[4].bind();
    buffers[4].allocate(normals_spheres.data(), static_cast<int>(normals_spheres.size()*sizeof(float)));
    rendering_program_spheres.bind();
    normalsLocation[0] = rendering_program_spheres.attributeLocation("normal");
    rendering_program_spheres.enableAttributeArray(normalsLocation[0]);
    rendering_program_spheres.setAttributeBuffer(normalsLocation[0],GL_FLOAT,0,3);
    buffers[4].release();


    buffers[0].bind();
    rendering_program_spheres.bind();
    centerLocation[0] = rendering_program_spheres.attributeLocation("center");
    rendering_program_spheres.enableAttributeArray(centerLocation[0]);
    rendering_program_spheres.setAttributeBuffer(centerLocation[0],GL_FLOAT,0,3);
    buffers[0].release();
    if(extension_is_found)
    {
        glVertexAttribDivisor(centerLocation[0],1);
        glVertexAttribDivisor(normalsLocation[0],0);
    }
    vao[3].release();
    rendering_program_spheres.release();

    rendering_program_cylinders.bind();
    vao[8].bind();
    buffers[8].bind();
    buffers[8].allocate(points_cube->data(), static_cast<int>(points_cube->size()*sizeof(float)));
    poly_vertexLocation[2] = rendering_program_cylinders.attributeLocation("vertex");
    rendering_program_cylinders.enableAttributeArray(poly_vertexLocation[2]);
    rendering_program_cylinders.setAttributeBuffer(poly_vertexLocation[2],GL_FLOAT,0,3);
    buffers[8].release();

    buffers[9].bind();
    buffers[9].allocate(normals_cylinder->data(), static_cast<int>(normals_cylinder->size()*sizeof(float)));
    normalsLocation[1] = rendering_program_cylinders.attributeLocation("normal");
    rendering_program_cylinders.enableAttributeArray(normalsLocation[1]);
    rendering_program_cylinders.setAttributeBuffer(normalsLocation[1],GL_FLOAT,0,3);
    buffers[9].release();

    buffers[10].bind();
    buffers[10].allocate(transfo1_cube.data(), static_cast<int>(transfo1_cube.size()*sizeof(float)));
    centerLocation[1] = rendering_program_cylinders.attributeLocation("transfo1");
    rendering_program_cylinders.enableAttributeArray(centerLocation[1]);
    rendering_program_cylinders.setAttributeBuffer(centerLocation[1],GL_FLOAT,0,4);
    buffers[10].release();

    buffers[11].bind();
    buffers[11].allocate(transfo2_cube.data(), static_cast<int>(transfo2_cube.size()*sizeof(float)));
    centerLocation[2] = rendering_program_cylinders.attributeLocation("transfo2");
    rendering_program_cylinders.enableAttributeArray(centerLocation[2]);
    rendering_program_cylinders.setAttributeBuffer(centerLocation[2],GL_FLOAT,0,4);
    buffers[11].release();


    buffers[12].bind();
    buffers[12].allocate(transfo3_cube.data(), static_cast<int>(transfo3_cube.size()*sizeof(float)));
    centerLocation[3] = rendering_program_cylinders.attributeLocation("transfo3");
    rendering_program_cylinders.enableAttributeArray(centerLocation[3]);
    rendering_program_cylinders.setAttributeBuffer(centerLocation[3],GL_FLOAT,0,4);
    buffers[12].release();

    buffers[13].bind();
    buffers[13].allocate(transfo4_cube.data(), static_cast<int>(transfo4_cube.size()*sizeof(float)));
    centerLocation[4] = rendering_program_cylinders.attributeLocation("transfo4");
    rendering_program_cylinders.enableAttributeArray(centerLocation[4]);
    rendering_program_cylinders.setAttributeBuffer(centerLocation[4],GL_FLOAT,0,4);
    buffers[13].release();

    if(extension_is_found)
    {
        glVertexAttribDivisor(centerLocation[1],1);
        glVertexAttribDivisor(centerLocation[2],1);
        glVertexAttribDivisor(centerLocation[3],1);
        glVertexAttribDivisor(centerLocation[4],1);

        glVertexAttribDivisor(normalsLocation[1],0);
    }

    vao[8].release();
    vao[9].bind();
    buffers[14].bind();
    buffers[14].allocate(points_cylinder->data(), static_cast<int>(points_cylinder->size()*sizeof(float)));
    poly_vertexLocation[2] = rendering_program_cylinders.attributeLocation("vertex");
    rendering_program_cylinders.enableAttributeArray(poly_vertexLocation[2]);
    rendering_program_cylinders.setAttributeBuffer(poly_vertexLocation[2],GL_FLOAT,0,3);
    buffers[14].release();

    buffers[9].bind();
    normalsLocation[1] = rendering_program_cylinders.attributeLocation("normal");
    rendering_program_cylinders.enableAttributeArray(normalsLocation[1]);
    rendering_program_cylinders.setAttributeBuffer(normalsLocation[1],GL_FLOAT,0,3);
    buffers[9].release();

    buffers[15].bind();
    buffers[15].allocate(transfo1_cylinder.data(), static_cast<int>(transfo1_cylinder.size()*sizeof(float)));
    centerLocation[1] = rendering_program_cylinders.attributeLocation("transfo1");
    rendering_program_cylinders.enableAttributeArray(centerLocation[1]);
    rendering_program_cylinders.setAttributeBuffer(centerLocation[1],GL_FLOAT,0,4);
    buffers[15].release();

    buffers[16].bind();
    buffers[16].allocate(transfo2_cylinder.data(), static_cast<int>(transfo2_cylinder.size()*sizeof(float)));
    centerLocation[2] = rendering_program_cylinders.attributeLocation("transfo2");
    rendering_program_cylinders.enableAttributeArray(centerLocation[2]);
    rendering_program_cylinders.setAttributeBuffer(centerLocation[2],GL_FLOAT,0,4);
    buffers[16].release();


    buffers[17].bind();
    buffers[17].allocate(transfo3_cylinder.data(), static_cast<int>(transfo3_cylinder.size()*sizeof(float)));
    centerLocation[3] = rendering_program_cylinders.attributeLocation("transfo3");
    rendering_program_cylinders.enableAttributeArray(centerLocation[3]);
    rendering_program_cylinders.setAttributeBuffer(centerLocation[3],GL_FLOAT,0,4);
    buffers[17].release();

    buffers[18].bind();
    buffers[18].allocate(transfo4_cylinder.data(), static_cast<int>(transfo4_cylinder.size()*sizeof(float)));
    centerLocation[4] = rendering_program_cylinders.attributeLocation("transfo4");
    rendering_program_cylinders.enableAttributeArray(centerLocation[4]);
    rendering_program_cylinders.setAttributeBuffer(centerLocation[4],GL_FLOAT,0,4);
    buffers[18].release();

    if(extension_is_found)
    {
        glVertexAttribDivisor(centerLocation[1],1);
        glVertexAttribDivisor(centerLocation[2],1);
        glVertexAttribDivisor(centerLocation[3],1);
        glVertexAttribDivisor(centerLocation[4],1);

        glVertexAttribDivisor(normalsLocation[1],0);
    }

    vao[9].release();

    vao[10].bind();
    buffers[8].bind();
    poly_vertexLocation[2] = rendering_program_cylinders.attributeLocation("vertex");
    rendering_program_cylinders.enableAttributeArray(poly_vertexLocation[2]);
    rendering_program_cylinders.setAttributeBuffer(poly_vertexLocation[2],GL_FLOAT,0,3);
    buffers[8].release();

    buffers[9].bind();
    normalsLocation[1] = rendering_program_cylinders.attributeLocation("normal");
    rendering_program_cylinders.enableAttributeArray(normalsLocation[1]);
    rendering_program_cylinders.setAttributeBuffer(normalsLocation[1],GL_FLOAT,0,3);
    buffers[9].release();

    buffers[19].bind();
    buffers[19].allocate(transfo1_square.data(), static_cast<int>(transfo1_square.size()*sizeof(float)));
    centerLocation[1] = rendering_program_cylinders.attributeLocation("transfo1");
    rendering_program_cylinders.enableAttributeArray(centerLocation[1]);
    rendering_program_cylinders.setAttributeBuffer(centerLocation[1],GL_FLOAT,0,4);
    buffers[19].release();

    buffers[20].bind();
    buffers[20].allocate(transfo2_square.data(), static_cast<int>(transfo2_square.size()*sizeof(float)));
    centerLocation[2] = rendering_program_cylinders.attributeLocation("transfo2");
    rendering_program_cylinders.enableAttributeArray(centerLocation[2]);
    rendering_program_cylinders.setAttributeBuffer(centerLocation[2],GL_FLOAT,0,4);
    buffers[20].release();


    buffers[21].bind();
    buffers[21].allocate(transfo3_square.data(), static_cast<int>(transfo3_square.size()*sizeof(float)));
    centerLocation[3] = rendering_program_cylinders.attributeLocation("transfo3");
    rendering_program_cylinders.enableAttributeArray(centerLocation[3]);
    rendering_program_cylinders.setAttributeBuffer(centerLocation[3],GL_FLOAT,0,4);
    buffers[21].release();

    buffers[22].bind();
    buffers[22].allocate(transfo4_square.data(), static_cast<int>(transfo4_square.size()*sizeof(float)));
    centerLocation[4] = rendering_program_cylinders.attributeLocation("transfo4");
    rendering_program_cylinders.enableAttributeArray(centerLocation[4]);
    rendering_program_cylinders.setAttributeBuffer(centerLocation[4],GL_FLOAT,0,4);
    buffers[22].release();

    if(extension_is_found)
    {
        glVertexAttribDivisor(centerLocation[1],1);
        glVertexAttribDivisor(centerLocation[2],1);
        glVertexAttribDivisor(centerLocation[3],1);
        glVertexAttribDivisor(centerLocation[4],1);

        glVertexAttribDivisor(normalsLocation[1],0);
    }
    vao[10].release();

    rendering_program_cylinders.release();
    are_buffers_initialized = true;

}

void Scene::attrib_buffers(QGLViewer* viewer)
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


    rendering_program.bind();
    mvpLocation[0] = rendering_program.uniformLocation("mvp_matrix");
    colorLocation[0] = rendering_program.uniformLocation("color");
    rendering_program.setUniformValue(mvpLocation[0], mvpMatrix);

    rendering_program.release();


    rendering_program_spheres.bind();
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

void Scene::init() {
    // undo from QGLViewer internal initializeGL function
    // glDisable(GL_COLOR_MATERIAL);
    initializeOpenGLFunctions();
    glDrawArraysInstanced = (PFNGLDRAWARRAYSINSTANCEDARBPROC)ui->viewer->context()->getProcAddress("glDrawArraysInstancedARB");
    if(!glDrawArraysInstanced)
    {
        qDebug()<<"glDrawArraysInstancedARB : extension not found. Wireframe mode will be forced.";
        extension_is_found = false;
    }
    else
        extension_is_found = true;

    glVertexAttribDivisor = (PFNGLVERTEXATTRIBDIVISORARBPROC)ui->viewer->context()->getProcAddress("glVertexAttribDivisorARB");
    if(!glDrawArraysInstanced)
    {
        qDebug()<<"glVertexAttribDivisorARB : extension not found. Wireframe mode will be forced.";
        extension_is_found = false;
    }
    else
        extension_is_found = true;

    wireframe = !extension_is_found;
    // camera
    // only 2.7 gets an 'f' as VC++ warns if we don't
    ui->viewer->camera()->setPosition(Vec(0.5,0.5,2.7f));
    ui->viewer->camera()->lookAt(Vec(0.5,0.5,0.5));

    // scene inits
    ui->viewer->setSceneCenter(qglviewer::Vec(0.5,0.5,0.5));
    ui->viewer->setSceneRadius(2.0);
    ui->viewer->setBackgroundColor(Qt::white);
    ui->viewer->setForegroundColor(Qt::red);

    // OpenGL inits
    glPointSize(10.0);
    glLineWidth(1.0);
    glEnable(GL_POINT_SMOOTH);
    glEnable(GL_LINE_SMOOTH);

    // Scene OpenGL state
    compile_shaders();
    init_scene(EMPTY);


}

// Draws the triangulation
void Scene::draw() {
    glEnable(GL_DEPTH_TEST);
    if(!are_buffers_initialized)
        initialize_buffers();
    gl_draw_location();

    gl_draw_conflict();

    //// Draw the triangulation itself that is stored in the list.

    if(wireframe)
    {
        //draw the points
        vao[0].bind();
        change_material(materials[VERTEX_COLOR]);
        attrib_buffers(ui->viewer);
        rendering_program.bind();

        glPointSize(5);
        ::glEnable(GL_POINT_SMOOTH);

        rendering_program.setUniformValue(colorLocation[0], color);
        glDrawArrays(GL_POINTS, 0, static_cast<GLsizei>(pos_points.size()/3));
        rendering_program.release();
        vao[0].release();
        //draw the moving ball
        if (flying_ball) {
            vao[4].bind();
            change_material(materials[FLYING_BALL_COLOR]);
            attrib_buffers(ui->viewer);
            rendering_program.bind();

            glPointSize(5);
            ::glEnable(GL_POINT_SMOOTH);
            rendering_program.setUniformValue(colorLocation[0], color);
            glDrawArrays(GL_POINTS, 0, 1);
            rendering_program.release();
            vao[4].release();
        }
        //draw the lines
        vao[1].bind();
        change_material(materials[EDGE_COLOR]);
        attrib_buffers(ui->viewer);
        rendering_program.bind();

        rendering_program.setUniformValue(colorLocation[0], color);
        glDrawArrays(GL_LINES, 0,  static_cast<GLsizei>(pos_tube.size()/3));
        rendering_program.release();
        vao[1].release();

        //draw the cube
        if(!in_plane)
        {
            vao[2].bind();
            change_material(materials[DOMAIN_COLOR]);
            attrib_buffers(ui->viewer);
            rendering_program.bind();
            rendering_program.setUniformValue(colorLocation[0], color);
            glDrawArrays(GL_LINES, 0,  static_cast<GLsizei>(pos_cube.size()/3));
            rendering_program.release();
            vao[2].release();
        }
        else
        {
            vao[11].bind();
            change_material(materials[DOMAIN_COLOR]);
            attrib_buffers(ui->viewer);
            rendering_program.bind();
            rendering_program.setUniformValue(colorLocation[0], color);
            glDrawArrays(GL_LINES, 0,  static_cast<GLsizei>(pos_square.size()/3));
            rendering_program.release();
            vao[11].release();
        }
    }
    else
    {
        if(!in_plane)
        {
            //cube
            vao[8].bind();
            change_material(materials[DOMAIN_COLOR]);
            attrib_buffers(ui->viewer);
            rendering_program_cylinders.bind();
            glDrawArraysInstanced(GL_TRIANGLES, 0,  static_cast<GLsizei>(points_cube->size()/3),  static_cast<GLsizei>(transfo1_cube.size()/4));
            rendering_program_cylinders.release();
            vao[8].release();
        }
        else
        {
            //square
            vao[10].bind();
            change_material(materials[DOMAIN_COLOR]);
            attrib_buffers(ui->viewer);
            rendering_program_cylinders.bind();
            glDrawArraysInstanced(GL_TRIANGLES, 0,  static_cast<GLsizei>(points_cube->size()/3),  static_cast<GLsizei>(transfo1_square.size()/4));
            rendering_program_cylinders.release();
            vao[10].release();
        }
        //draw the spheres
        vao[3].bind();
        change_material(materials[VERTEX_COLOR]);
        attrib_buffers(ui->viewer);
        rendering_program_spheres.bind();
        glDrawArraysInstanced(GL_TRIANGLES, 0,  static_cast<GLsizei>(points_spheres.size()/3),  static_cast<GLsizei>(pos_points.size()/3));
        rendering_program_spheres.release();
        vao[3].release();

        //draw the moving ball
        if (flying_ball) {
            vao[7].bind();
            change_material(materials[FLYING_BALL_COLOR]);
            attrib_buffers(ui->viewer);
            rendering_program_spheres.bind();
            glDrawArraysInstanced(GL_TRIANGLES, 0,  static_cast<GLsizei>(points_spheres.size()/3),1);
            rendering_program_spheres.release();
            vao[7].release();
        }
        //draw the triangulation
        vao[9].bind();
        change_material(materials[EDGE_COLOR]);
        attrib_buffers(ui->viewer);
        rendering_program_cylinders.bind();
        glDrawArraysInstanced(GL_TRIANGLES, 0,  static_cast<GLsizei>(points_cylinder->size()/3),  static_cast<GLsizei>(transfo1_cylinder.size()/4));
        rendering_program_cylinders.release();
        vao[9].release();
    }
    //draw the triangles
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
    if(dlocate)
    {
        vao[5].bind();
        rendering_program.bind();
        color.setRgbF(0.f, 0.f , 0.5f, 0.5f);
        rendering_program.setUniformValue(colorLocation[0], color);
        glDrawArrays(GL_TRIANGLES, 0,  static_cast<GLsizei>(pos_location.size()/3));
        rendering_program.release();
        vao[5].release();
    }
    if(dconflict)
    {
        vao[6].bind();
        rendering_program.bind();
        color.setRgbF(0.69f, 0.18f , 0.26f, 0.6f);
        rendering_program.setUniformValue(colorLocation[0], color);
        glDrawArrays(GL_TRIANGLES, 0,  static_cast<GLsizei>(pos_conflict.size()/3));
        rendering_program.release();
        vao[6].release();
    }
    glDisable(GL_BLEND);


}

void Scene::load_points(const QString& fileName) {
    p3dt.clear();
    std::vector<Point> points;
    std::ifstream ifs(fileName.toLatin1().data() );
    std::copy(std::istream_iterator<Point>(ifs),
              std::istream_iterator<Point>(),
              std::back_inserter(points));
    std::random_shuffle(points.begin(), points.end());
    p3dt.insert(points.begin(), points.end());

    QString snv;
    int nv = static_cast<int>(p3dt.number_of_vertices());
    snv.setNum(nv);
    changed();
    Q_EMIT message(QString("|V| = ") + snv, 0);


    draw();
}

// update the position of the moving point
void Scene::update_position()
{
    double x = moving_point.x() +0.01023;
    double y = moving_point.y() +0.003123;
    double z = (in_plane ? 0.0 : moving_point.z() +0.02567);
    if(x>1.)x-=1.;
    if(y>1.)y-=1.;
    if(z>1.)z-=1.;
    moving_point = Point(x,y,z);
    float moving_ball[] = {float(moving_point.x()),
                           float(moving_point.y()),
                           float(moving_point.z())};
    vao[4].bind();
    buffers[5].bind();
    buffers[5].allocate(moving_ball, 3*sizeof(float));
    rendering_program.bind();
    poly_vertexLocation[0] = rendering_program.attributeLocation("vertex");
    rendering_program.enableAttributeArray(poly_vertexLocation[0]);
    rendering_program.setAttributeBuffer(poly_vertexLocation[0],GL_FLOAT,0,3);
    buffers[5].release();

    rendering_program.release();
    vao[4].release();

    vao[7].bind();
    buffers[3].bind();
    poly_vertexLocation[1] = rendering_program_spheres.attributeLocation("vertex");
    rendering_program_spheres.enableAttributeArray(poly_vertexLocation[1]);
    rendering_program_spheres.setAttributeBuffer(poly_vertexLocation[1],GL_FLOAT,0,3);
    buffers[3].release();

    buffers[4].bind();
    rendering_program_spheres.bind();
    normalsLocation[1] = rendering_program_spheres.attributeLocation("normal");
    rendering_program_spheres.enableAttributeArray(normalsLocation[1]);
    rendering_program_spheres.setAttributeBuffer(normalsLocation[1],GL_FLOAT,0,3);
    buffers[4].release();


    buffers[5].bind();
    rendering_program_spheres.bind();
    centerLocation[0] = rendering_program_spheres.attributeLocation("center");
    rendering_program_spheres.enableAttributeArray(centerLocation[0]);
    rendering_program_spheres.setAttributeBuffer(centerLocation[0],GL_FLOAT,0,3);
    buffers[5].release();
    if(extension_is_found)
    {
        glVertexAttribDivisor(centerLocation[0],1);
        glVertexAttribDivisor(normalsLocation[1],0);
    }
    vao[7].release();

    ui->viewer->update();
}

// some initialization templates
void Scene::init_scene(Init ID) {
    bool temp_flags[] = {dlocate, dconflict};
    dlocate = false;
    dconflict = false;
    p3dt.clear();
    RandPts rp(0.5);
    Point pt2;
    switch (ID) {
    case GRID:
        p3dt.insert_dummy_points();
        break;
    case SINGLE:
        p3dt.insert(Point(0.3,0.4,0.5));
        break;
    case PLANE:
        for (int i=0 ; i<10 ; i++) {
            pt2 = *rp+Vector(0.5,0.5,0.5);
            rp++;
            p3dt.insert(Point(pt2.x(),pt2.y(),0.0));
        }
        break;
    case RANDOM:
        do {
            p3dt.insert(*rp+Vector(0.5,0.5,0.5));
            rp++;
        }
        while (p3dt.number_of_vertices()<30);
    default:
        break;
    }
    dlocate = temp_flags[0];
    dconflict = temp_flags[1];
    changed();
}

// get the offset that is common to all points of a triangle in the
// triangulation
inline void Scene::get_tri_offsets(const Cell_handle ch, int i,
                                   Offset &off0, Offset &off1, Offset &off2) const {
    off0 = p3dt.get_offset(ch,(i+1)&3);
    off1 = p3dt.get_offset(ch,(i+2)&3);
    off2 = p3dt.get_offset(ch,(i+3)&3);
    if (it_type == P3DT::UNIQUE || it_type == P3DT::UNIQUE_COVER_DOMAIN) {
        int diff_offx = (std::min)((std::min)(off0.x(),off1.x()),off2.x());
        int diff_offy = (std::min)((std::min)(off0.y(),off1.y()),off2.y());
        int diff_offz = (std::min)((std::min)(off0.z(),off1.z()),off2.z());
        Offset diff_off(diff_offx, diff_offy, diff_offz);
        off0 -= diff_off;
        off1 -= diff_off;
        off2 -= diff_off;
    }
}

// get the offset that is common to all points of a tetrahedron in the
// triangulation
inline void Scene::get_tet_offsets(const Cell_handle ch,
                                   Offset &off0, Offset &off1, Offset &off2, Offset &off3) const {
    off0 = p3dt.get_offset(ch,0);
    off1 = p3dt.get_offset(ch,1);
    off2 = p3dt.get_offset(ch,2);
    off3 = p3dt.get_offset(ch,3);
    if (it_type == P3DT::UNIQUE || it_type == P3DT::UNIQUE_COVER_DOMAIN) {
        int diff_offx = (std::min)((std::min)(off0.x(),off1.x()),
                                   (std::min)(off2.x(),off3.x()));
        int diff_offy = (std::min)((std::min)(off0.y(),off1.y()),
                                   (std::min)(off2.y(),off3.y()));
        int diff_offz = (std::min)((std::min)(off0.z(),off1.z()),
                                   (std::min)(off2.z(),off3.z()));
        Offset diff_off(diff_offx, diff_offy, diff_offz);
        off0 -= diff_off;
        off1 -= diff_off;
        off2 -= diff_off;
        off3 -= diff_off;
    }
}

// return an integer that encodes the translations which have to be
// applied to the triangle to draw
inline int Scene::get_tri_drawing_offsets(const Cell_handle ch, int i) const {
    Offset off0, off1, off2;
    // if drawing boundary cells multiply is not activated then there is
    // nothing to do.
    switch( it_type ) {
    case P3DT::UNIQUE_COVER_DOMAIN:
        get_tri_offsets(ch,i,off0,off1,off2);
        break;
    case P3DT::STORED_COVER_DOMAIN:
        off0 = p3dt.int_to_off(ch->offset((i+1)&3));
        off1 = p3dt.int_to_off(ch->offset((i+2)&3));
        off2 = p3dt.int_to_off(ch->offset((i+3)&3));
        break;
    default:
        return 0;
    }

    CGAL_assertion(off0.x() == 0 || off0.x() == 1);
    CGAL_assertion(off0.y() == 0 || off0.y() == 1);
    CGAL_assertion(off0.z() == 0 || off0.z() == 1);
    CGAL_assertion(off1.x() == 0 || off1.x() == 1);
    CGAL_assertion(off1.y() == 0 || off1.y() == 1);
    CGAL_assertion(off1.z() == 0 || off1.z() == 1);
    CGAL_assertion(off2.x() == 0 || off2.x() == 1);
    CGAL_assertion(off2.y() == 0 || off2.y() == 1);
    CGAL_assertion(off2.z() == 0 || off2.z() == 1);

    int offx = ( ((off0.x() == 0 && off1.x() == 0 && off2.x() == 0)
                  || (off0.x() == 1 && off1.x() == 1 && off2.x() == 1)) ? 0 : 1);
    int offy = ( ((off0.y() == 0 && off1.y() == 0 && off2.y() == 0)
                  || (off0.y() == 1 && off1.y() == 1 && off2.y() == 1)) ? 0 : 1);
    int offz = ( ((off0.z() == 0 && off1.z() == 0 && off2.z() == 0)
                  || (off0.z() == 1 && off1.z() == 1 && off2.z() == 1)) ? 0 : 1);

    return( 4*offx + 2*offy + offz );
}

// return an integer that encodes the translations which have to be
// applied to the tetrahedron to draw
inline int Scene::get_tet_drawing_offsets(const Cell_handle ch) const {
    Offset off0, off1, off2, off3;
    // if drawing boundary cells multiply is not activated then there is
    // nothing to do.
    switch( it_type ) {
    case P3DT::UNIQUE_COVER_DOMAIN:
        get_tet_offsets(ch,off0,off1,off2,off3);
        break;
    case P3DT::STORED_COVER_DOMAIN:
        off0 = p3dt.int_to_off(ch->offset(0));
        off1 = p3dt.int_to_off(ch->offset(1));
        off2 = p3dt.int_to_off(ch->offset(2));
        off3 = p3dt.int_to_off(ch->offset(3));
        break;
    default:
        return 0;
    }

    CGAL_assertion(off0.x() == 0 || off0.x() == 1);
    CGAL_assertion(off0.y() == 0 || off0.y() == 1);
    CGAL_assertion(off0.z() == 0 || off0.z() == 1);
    CGAL_assertion(off1.x() == 0 || off1.x() == 1);
    CGAL_assertion(off1.y() == 0 || off1.y() == 1);
    CGAL_assertion(off1.z() == 0 || off1.z() == 1);
    CGAL_assertion(off2.x() == 0 || off2.x() == 1);
    CGAL_assertion(off2.y() == 0 || off2.y() == 1);
    CGAL_assertion(off2.z() == 0 || off2.z() == 1);
    CGAL_assertion(off3.x() == 0 || off3.x() == 1);
    CGAL_assertion(off3.y() == 0 || off3.y() == 1);
    CGAL_assertion(off3.z() == 0 || off3.z() == 1);

    int offx = ( ((off0.x() == 0 && off1.x() == 0
                   && off2.x() == 0 && off3.x() == 0)
                  || (off0.x() == 1 && off1.x() == 1
                      && off2.x() == 1 && off3.x() == 1)) ? 0 : 1);
    int offy = ( ((off0.y() == 0 && off1.y() == 0
                   && off2.y() == 0 && off3.y() == 0)
                  || (off0.y() == 1 && off1.y() == 1
                      && off2.y() == 1 && off3.y() == 1)) ? 0 : 1);
    int offz = ( ((off0.z() == 0 && off1.z() == 0
                   && off2.z() == 0 && off3.z() == 0)
                  || (off0.z() == 1 && off1.z() == 1
                      && off2.z() == 1 && off3.z() == 1)) ? 0 : 1);

    return( 4*offx + 2*offy + offz );
}

// construct a triangle from a given facet, given vertex offsets and a
// common offset
inline Triangle Scene::construct_triangle(const Cell_handle ch, int i,
                                          const Offset& off0, const Offset& off1, const Offset& off2, int off) const {
    if (it_type == P3DT::STORED || it_type == P3DT::UNIQUE) {
        CGAL_assertion( off == 0 );
        return p3dt.construct_triangle(
                    ch->vertex((i+1)&3)->point(), ch->vertex((i+2)&3)->point(),
                    ch->vertex((i+3)&3)->point(), off0, off1, off2);
    }
    Offset diff_off((off>>2)&1,(off>>1)&1,off&1);
    switch (it_type) {
    case P3DT::STORED_COVER_DOMAIN:
        return p3dt.construct_triangle(
                    ch->vertex((i+1)&3)->point(), ch->vertex((i+2)&3)->point(),
                    ch->vertex((i+3)&3)->point(),
                    p3dt.combine_offsets(off0,-diff_off),
                    p3dt.combine_offsets(off1,-diff_off),
                    p3dt.combine_offsets(off2,-diff_off));
        break;
    case P3DT::UNIQUE_COVER_DOMAIN:
        return p3dt.construct_triangle(
                    ch->vertex((i+1)&3)->point(), ch->vertex((i+2)&3)->point(),
                    ch->vertex((i+3)&3)->point(),
                    off0-diff_off, off1-diff_off, off2-diff_off);
        break;
    default:
        CGAL_assertion(false);
        return Triangle();
    }
}

// construct a triangle from a given cell, given vertex offsets and a
// common offset
inline Tetrahedron Scene::construct_tetrahedron(const Cell_handle ch,
                                                const Offset& off0, const Offset& off1, const Offset& off2,
                                                const Offset& off3, int off) const {
    if (it_type == P3DT::STORED || it_type == P3DT::UNIQUE) {
        CGAL_assertion( off == 0 );
        return p3dt.construct_tetrahedron(
                    ch->vertex(0)->point(), ch->vertex(1)->point(),
                    ch->vertex(2)->point(), ch->vertex(3)->point(),
                    off0, off1, off2, off3);
    }
    Offset diff_off((off>>2)&1,(off>>1)&1,off&1);
    switch (it_type) {
    case P3DT::STORED_COVER_DOMAIN:
        return p3dt.construct_tetrahedron(
                    ch->vertex(0)->point(), ch->vertex(1)->point(),
                    ch->vertex(2)->point(), ch->vertex(3)->point(),
                    p3dt.combine_offsets(off0,-diff_off),
                    p3dt.combine_offsets(off1,-diff_off),
                    p3dt.combine_offsets(off2,-diff_off),
                    p3dt.combine_offsets(off3,-diff_off));
        break;
    case P3DT::UNIQUE_COVER_DOMAIN:
        return p3dt.construct_tetrahedron(
                    ch->vertex(0)->point(), ch->vertex(1)->point(),
                    ch->vertex(2)->point(), ch->vertex(3)->point(),
                    off0-diff_off, off1-diff_off, off2-diff_off, off3-diff_off);
        break;
    default:
        CGAL_assertion(false);
        return Tetrahedron();
    }
}


// collect primitives (segments, triangles, tetrahedra) from the
// triangulation using the geometric iterators and store them in the
// given segment set
inline void Scene::primitives_from_geom_it(Segment_set& sset) {
    Point p0,p1,p2,p3;
    switch(draw_type) {
    case SEGMENT:
        for ( Segment_iterator sit = p3dt.periodic_segments_begin(it_type) ;
              sit != p3dt.periodic_segments_end(it_type) ; ++sit ) {
            sset.insert(p3dt.segment(*sit));
        }
        break;
    case TRIANGLE:
        for ( Triangle_iterator tit = p3dt.periodic_triangles_begin(it_type) ;
              tit != p3dt.periodic_triangles_end(it_type) ; ++tit ) {
            p0 = p3dt.point(tit->at(0));
            p1 = p3dt.point(tit->at(1));
            p2 = p3dt.point(tit->at(2));
            sset.insert(p0 < p1 ? Segment(p0,p1) : Segment(p1,p0));
            sset.insert(p0 < p2 ? Segment(p0,p2) : Segment(p2,p0));
            sset.insert(p1 < p2 ? Segment(p1,p2) : Segment(p2,p1));
        }
        break;
    case TETRAHEDRON:
        for ( Tetrahedron_iterator tit = p3dt.periodic_tetrahedra_begin(it_type) ;
              tit != p3dt.periodic_tetrahedra_end(it_type) ; ++tit ) {
            p0 = p3dt.point(tit->at(0));
            p1 = p3dt.point(tit->at(1));
            p2 = p3dt.point(tit->at(2));
            p3 = p3dt.point(tit->at(3));
            sset.insert((p0 < p1) ? Segment(p0,p1) : Segment(p1,p0));
            sset.insert((p0 < p2) ? Segment(p0,p2) : Segment(p2,p0));
            sset.insert((p0 < p3) ? Segment(p0,p3) : Segment(p3,p0));
            sset.insert((p1 < p2) ? Segment(p1,p2) : Segment(p2,p1));
            sset.insert((p1 < p3) ? Segment(p1,p3) : Segment(p3,p1));
            sset.insert((p2 < p3) ? Segment(p2,p3) : Segment(p3,p2));
        }
        break;
    }
}

// clip segments from the given segment set that are partially outside
// of the unit cube/square. Eliminate those who are completely outside
inline void Scene::segment_clipping(Segment_set& sset) {
    Segment_clipper clipper;
    Segment_set sset_tmp;
    for (Segment_set::iterator it = sset.begin() ; it != sset.end() ; ++it) {
        Point s = it->source();
        Point t = it->target();
        if (clipper(s,t)) sset_tmp.insert((s<t?Segment(s,t):Segment(t,s)));
    }
    std::swap(sset, sset_tmp);
}

// clip segments from the given segment set that are partially outside
// of the unit cube/square. Draw their outside part in a different
// color.
// TODO: don't eliminate segments that are completely outside but draw
// them in the different color as well
inline void Scene::segment_2color_clipping (Segment_set& sset) {
    Segment_clipper clipper;
    Segment_set sset_tmp, sset_out;
    for (Segment_set::iterator it = sset.begin() ; it != sset.end() ; ++it) {
        Point s = it->source();
        Point t = it->target();
        if (clipper(s,t)) {
            sset_tmp.insert((s<t?Segment(s,t):Segment(t,s)));
            Point p = it->source();
            Point q = it->target();
            if (Segment(p,s).squared_length() > Segment(p,t).squared_length())
                std::swap(s,t);
            if (p!=s) sset_out.insert((p<s?Segment(p,s):Segment(s,p)));
            if (q!=t) sset_out.insert((q<t?Segment(q,t):Segment(t,q)));
        }
    }


    std::swap(sset, sset_tmp);
}

// Draw the faces of the tetrahedron in which the moving point is currently
// located transparently. It depends on it_type which periodic copies
// of the respective cell will be drawn. In general it will be all
// cells that occur in the draw list
void Scene::gl_draw_location() {
    pos_location.resize(0);
    if (p3dt.number_of_vertices() == 0) return;
    // Do the point location
    Cell_handle ch = p3dt.locate(moving_point);
    std::vector<Projected_triangle> cf;

    // Transparency

    if (in_plane) {
        int i=0;
        int count = 0;
        // Figure out whether there is a facet that is completly contained
        // in the z=0 plane
        for (int j=0 ; j<4 ; j++) {
            if (ch->vertex(j)->point().z() != 0.0 ||
                    p3dt.get_offset(ch,j).z() != 0) {
                i=j;
                count++;
            }
        }
        // If so, compute its triangle(s) and insert it in cf
        if (count==1) {
            Offset off0, off1, off2;
            get_tri_offsets(ch, i, off0, off1, off2);
            int diff_off = get_tri_drawing_offsets(ch, i);
            for (int offs=0 ; offs<=diff_off ; offs++) {
                if ((((~offs)|diff_off)&7)!=7) continue;
                Triangle tri_to_draw = construct_triangle(ch,i,off0,off1,off2,offs);
                Point p = tri_to_draw.vertex(0);
                Point q = tri_to_draw.vertex(1);
                Point r = tri_to_draw.vertex(2);
                cf.push_back(Projected_triangle(.0,Triangle(p,q,r)));
            }
        }
    } else {
        double modelMatrix[16];
        double projMatrix[16];
        int viewport[4];
        glGetDoublev(GL_MODELVIEW_MATRIX, modelMatrix);
        glGetDoublev(GL_PROJECTION_MATRIX, projMatrix);
        glGetIntegerv(GL_VIEWPORT, viewport);

        // Compute the triangles that are the facets of the cell and
        // insert them in cf
        Offset off0,off1,off2,off3;
        get_tet_offsets(ch, off0, off1, off2, off3);
        int diff_off = get_tet_drawing_offsets(ch);

        for (int offs=0 ; offs<=diff_off ; offs++) {
            if ((((~offs)|diff_off)&7)!=7) continue;
            Tetrahedron tet_to_draw = construct_tetrahedron(
                        ch, off0, off1, off2, off3, offs);

            for(int i=0; i < 4; i++){
                Point p = tet_to_draw.vertex((i+1)&3);
                Point q = tet_to_draw.vertex((i+2)&3);
                Point r = tet_to_draw.vertex((i+3)&3);
                Vector c= (Vector(Point(),p)+Vector(Point(),q)+Vector(Point(),r))/3.;
                Point cp = Point(c.x(),c.y(),c.z());
                // project facet center
                double px,py,pz;
                gluProject(cp.x(),cp.y(),cp.z(),
                           modelMatrix, projMatrix, viewport,
                           &px,&py,&pz);
                cf.push_back(Projected_triangle(pz,Triangle(p,q,r)));
            }
        }

        // Sort cf according to their z coordinates to enable transparency
        std::sort(cf.begin(), cf.end(), Projected_triangle::closer);
    }

    // Draw all triangles from cf
    for (std::vector<Projected_triangle >::iterator cfit = cf.begin() ; cfit != cf.end() ; cfit++) {
        Point p = cfit->t().vertex(0);
        Point q = cfit->t().vertex(1);
        Point r = cfit->t().vertex(2);
        pos_location.push_back(p.x()); pos_location.push_back(p.y()); pos_location.push_back(p.z());
        pos_location.push_back(q.x()); pos_location.push_back(q.y()); pos_location.push_back(q.z());
        pos_location.push_back(r.x()); pos_location.push_back(r.y()); pos_location.push_back(r.z());

    }
    vao[5].bind();
    buffers[6].bind();
    buffers[6].allocate(pos_location.data(), static_cast<int>(pos_location.size()*sizeof(float)));
    rendering_program.bind();
    poly_vertexLocation[0] = rendering_program.attributeLocation("vertex");
    rendering_program.enableAttributeArray(poly_vertexLocation[0]);
    rendering_program.setAttributeBuffer(poly_vertexLocation[0],GL_FLOAT,0,3);
    buffers[6].release();

    rendering_program.release();
    vao[5].release();
}

// Draw the boundary faces of the current conflict region of the
// moving point transparently. It depends on it_type which periodic
// copies of the respective cell will be drawn. In general it will be
// all cells that occur in the draw list.
void Scene::gl_draw_conflict() {
    pos_conflict.resize(0);
    if (p3dt.number_of_vertices() == 0) return;
    Cell_handle ch;
    std::vector<Cell_handle> cic;
    std::vector<Facet> boundary_facets;
    // Find the conflict region
    Cell_handle c = p3dt.locate(moving_point);
    p3dt.find_conflicts(moving_point,c,std::back_inserter(boundary_facets),std::back_inserter(cic),CGAL::Emptyset_iterator());

    std::vector<Projected_triangle> bfm;

    // Transparency
    //glEnable(GL_BLEND);
    //glColor4f(.69f, 0.18f , 0.26f, 0.6f);
    //glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);

    if (in_plane) {
        for (unsigned int k=0 ; k<cic.size(); k++) {
            ch = cic[k];

            int i = 0;
            int count = 0;
            // Figure out whether there is a facet that is completely
            // contained in the z=0 plane
            for (int j=0 ; j<4 ; j++) {
                if (ch->vertex(j)->point().z() != 0.0 ||
                        p3dt.get_offset(ch,j).z() != 0) {
                    i=j;
                    count++;
                }
            }
            // If so, compute its triangle(s) and insert it in bfm
            if (count==1) {
                Offset off0, off1, off2;
                get_tri_offsets(ch, i, off0, off1, off2);
                int diff_off = get_tri_drawing_offsets(ch,i);
                for (int offs = 0 ; offs<=diff_off ; offs++) {
                    if ((((~offs)|diff_off)&7)!=7) continue;
                    Triangle tri_to_draw = construct_triangle(ch,i,off0,off1,off2,offs);
                    Point p = tri_to_draw.vertex(0);
                    Point q = tri_to_draw.vertex(1);
                    Point r = tri_to_draw.vertex(2);
                    bfm.push_back(Projected_triangle(.0,Triangle(p,q,r)));
                }
            }
        }
    } else {
        double modelMatrix[16];
        double projMatrix[16];
        int viewport[4];
        glGetDoublev(GL_MODELVIEW_MATRIX, modelMatrix);
        glGetDoublev(GL_PROJECTION_MATRIX, projMatrix);
        glGetIntegerv(GL_VIEWPORT, viewport);

        for (unsigned int i=0 ; i<boundary_facets.size(); i++) {
            ch = boundary_facets[i].first;
            int j=boundary_facets[i].second;

            // Compute the triangle(s) of the facet and insert them in bfm
            Offset off0, off1, off2;
            get_tri_offsets(ch, j, off0, off1, off2);
            int diff_off = get_tri_drawing_offsets(ch, j);
            for (int offs=0 ; offs<=diff_off ; offs++) {
                if ((((~offs)|diff_off)&7)!=7) continue;
                Triangle tri_to_draw = construct_triangle(ch,j,off0,off1,off2,offs);
                Point p = tri_to_draw.vertex(0);
                Point q = tri_to_draw.vertex(1);
                Point r = tri_to_draw.vertex(2);
                Vector c= (Vector(Point(),p)+Vector(Point(),q)+Vector(Point(),r))/3.;
                Point cp = Point(c.x(),c.y(),c.z());
                // project facet center
                double px,py,pz;
                gluProject(cp.x(),cp.y(),cp.z(),
                           modelMatrix, projMatrix, viewport,
                           &px,&py,&pz);
                bfm.push_back(Projected_triangle(pz,Triangle(p,q,r)));
            }
        }

        // Sort bfm according to their z coordinates to enable transparency
        std::sort(bfm.begin(), bfm.end(), Projected_triangle::closer);
    }

    // Draw all triangles from bfm
    // glBegin(GL_TRIANGLES);
    for (std::vector<Projected_triangle >::iterator bfmit = bfm.begin() ;
         bfmit != bfm.end() ; bfmit++) {
        Point p = bfmit->t().vertex(0);
        Point q = bfmit->t().vertex(1);
        Point r = bfmit->t().vertex(2);

        pos_conflict.push_back(p.x()); pos_conflict.push_back(p.y()); pos_conflict.push_back(p.z());
        pos_conflict.push_back(q.x()); pos_conflict.push_back(q.y()); pos_conflict.push_back(q.z());
        pos_conflict.push_back(r.x()); pos_conflict.push_back(r.y()); pos_conflict.push_back(r.z());

    }
    // glEnd();
    vao[6].bind();
    buffers[7].bind();
    buffers[7].allocate(pos_conflict.data(), static_cast<int>(pos_conflict.size()*sizeof(float)));
    rendering_program.bind();
    poly_vertexLocation[0] = rendering_program.attributeLocation("vertex");
    rendering_program.enableAttributeArray(poly_vertexLocation[0]);
    rendering_program.setAttributeBuffer(poly_vertexLocation[0],GL_FLOAT,0,3);
    buffers[7].release();

    rendering_program.release();
    vao[6].release();
    // glEnd();
    //
    // glDisable(GL_BLEND);
}

// provide some color constants for the GLU primitive appearance.
void Scene::change_material(const QString &string) {


    // Change
    if(string == "Silver")
    {
        // Ambient
        ambient[0] = 0.19225f;
        ambient[1] = 0.19225f;
        ambient[2] = 0.19225f;
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
    }

    else if(string == "Gold")
        {
            // Ambient
            ambient[0] = 0.24725f;
            ambient[1] = 0.1995f;
            ambient[2] = 0.0745f;
            ambient[3] = 1.0f;
            // Diffuse
            diffuse[0] = 0.75164f;
            diffuse[1] = 0.60648f;
            diffuse[2] = 0.22648f;
            diffuse[3] = 1.0f;
            // Specular
            specular[0] = 0.928281f;
            specular[1] = 0.855802f;
            specular[2] = 0.666065f;
            specular[3] = 1.0f;
            // Shininess
            shininess = 51.2f;
        }


    else if(string == "Red")
            {
                // Ambient
                ambient[0] = 0.75f;
                ambient[1] = 0.0f;
                ambient[2] = 0.0f;
                ambient[3] = 1.0f;
                // Diffuse
                diffuse[0] = 1.0f;
                diffuse[1] = 0.0f;
                diffuse[2] = 0.0f;
                diffuse[3] = 1.0f;
                // Specular
                specular[0] = 0.75f;
                specular[1] = 0.75f;
                specular[2] = 0.75f;
                specular[3] = 1.0f;
                // Shininess
                shininess = 64.0f;
            }

    else if(string == "Green")
                {
                    // Ambient
                    ambient[0] = 0.0225f;
                    ambient[1] = 0.19125f;
                    ambient[2] = 0.0735f;
                    ambient[3] = 1.0f;
                    // Diffuse
                    diffuse[0] = 0.0828f;
                    diffuse[1] = 0.3038f;
                    diffuse[2] = 0.14048f;
                    diffuse[3] = 1.0f;
                    // Specular
                    specular[0] = 0.086014f;
                    specular[1] = 0.306777f;
                    specular[2] = 0.117622f;
                    specular[3] = 1.0f;
                    // Shininess
                    shininess = 12.8f;
                }

    else if(string == "Black plastic")
                    {
                        // Ambient
                        ambient[0] = 0.0f;
                        ambient[1] = 0.0f;
                        ambient[2] = 0.0f;
                        ambient[3] = 1.0f;
                        // Diffuse
                        diffuse[0] = 0.01f;
                        diffuse[1] = 0.01f;
                        diffuse[2] = 0.01f;
                        diffuse[3] = 1.0f;
                        // Specular
                        specular[0] = 0.5f;
                        specular[1] = 0.5f;
                        specular[2] = 0.5f;
                        specular[3] = 1.0f;
                        // Shininess
                        shininess = 32.0f;
                    }

    color.setRgbF(diffuse[0],diffuse[1],diffuse[2],diffuse[3]);

}

void Scene::changed()
{
    compute_elements();
    are_buffers_initialized = false;
}

void Scene::draw_sphere(float R, int prec)
{

    points_spheres.resize(0);
    int rings=prec, sectors=prec;
    float T, P;
    float x[4],y[4],z[4];


    //Top of the sphere
    for(int t=0; t<360; t+=sectors)
    {

        points_spheres.push_back(0);
        points_spheres.push_back(0);
        points_spheres.push_back(R);


        normals_spheres.push_back(0);
        normals_spheres.push_back(0);
        normals_spheres.push_back(1);



        P = rings*M_PI/180.0;
        T = t*M_PI/180.0;
        x[1] = sin(P) * cos(T) ;
        y[1] = sin(P) * sin(T) ;
        z[1] = cos(P);
        points_spheres.push_back(R * x[1]);
        points_spheres.push_back(R * y[1]);
        points_spheres.push_back(R * z[1]);


        normals_spheres.push_back(x[1]);
        normals_spheres.push_back(y[1]);
        normals_spheres.push_back(z[1]);

        //
        P = rings*M_PI/180.0;
        T = (t+sectors)*M_PI/180.0;
        x[2] = sin(P) * cos(T) ;
        y[2] = sin(P) * sin(T) ;
        z[2] = cos(P);
        points_spheres.push_back(R * x[2]);
        points_spheres.push_back(R * y[2]);
        points_spheres.push_back(R * z[2]);

        normals_spheres.push_back(x[2]);
        normals_spheres.push_back(y[2]);
        normals_spheres.push_back(z[2]);

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

            points_spheres.push_back(R * x[0]);
            points_spheres.push_back(R * y[0]);
            points_spheres.push_back(R * z[0]);


            normals_spheres.push_back(x[0]);
            normals_spheres.push_back(y[0]);
            normals_spheres.push_back(z[0]);

            //B
            P = (p+rings)*M_PI/180.0;
            T = t*M_PI/180.0;
            x[1] = sin(P) * cos(T) ;
            y[1] = sin(P) * sin(T) ;
            z[1] = cos(P);
            points_spheres.push_back(R * x[1]);
            points_spheres.push_back(R * y[1]);
            points_spheres.push_back(R * z[1]);


            normals_spheres.push_back(x[1]);
            normals_spheres.push_back(y[1]);
            normals_spheres.push_back(z[1]);

            //C
            P = p*M_PI/180.0;
            T = (t+sectors)*M_PI/180.0;
            x[2] = sin(P) * cos(T) ;
            y[2] = sin(P) * sin(T) ;
            z[2] = cos(P);
            points_spheres.push_back(R * x[2]);
            points_spheres.push_back(R * y[2]);
            points_spheres.push_back(R * z[2]);


            normals_spheres.push_back(x[2]);
            normals_spheres.push_back(y[2]);
            normals_spheres.push_back(z[2]);
            //D
            P = (p+rings)*M_PI/180.0;
            T = (t+sectors)*M_PI/180.0;
            x[3] = sin(P) * cos(T) ;
            y[3] = sin(P) * sin(T) ;
            z[3] = cos(P);
            points_spheres.push_back(R * x[3]);
            points_spheres.push_back(R * y[3]);
            points_spheres.push_back(R * z[3]);


            normals_spheres.push_back(x[3]);
            normals_spheres.push_back(y[3]);
            normals_spheres.push_back(z[3]);



            points_spheres.push_back(R * x[1]);
            points_spheres.push_back(R * y[1]);
            points_spheres.push_back(R * z[1]);


            normals_spheres.push_back(x[1]);
            normals_spheres.push_back(y[1]);
            normals_spheres.push_back(z[1]);

            points_spheres.push_back(R * x[2]);
            points_spheres.push_back(R * y[2]);
            points_spheres.push_back(R * z[2]);


            normals_spheres.push_back(x[2]);
            normals_spheres.push_back(y[2]);
            normals_spheres.push_back(z[2]);

        }
    //Bottom of the sphere
    for(int t=0; t<360; t+=sectors)
    {


        points_spheres.push_back(0);
        points_spheres.push_back(0);
        points_spheres.push_back(-R);


        normals_spheres.push_back(0);
        normals_spheres.push_back(0);
        normals_spheres.push_back(-1);


        P = (180-rings)*M_PI/180.0;
        T = t*M_PI/180.0;
        x[1] = sin(P) * cos(T) ;
        y[1] = sin(P) * sin(T) ;
        z[1] = cos(P);
        points_spheres.push_back(R * x[1]);
        points_spheres.push_back(R * y[1]);
        points_spheres.push_back(R * z[1]);


        normals_spheres.push_back(x[1]);
        normals_spheres.push_back(y[1]);
        normals_spheres.push_back(z[1]);


        P = (180-rings)*M_PI/180.0;
        T = (t+sectors)*M_PI/180.0;
        x[2] = sin(P) * cos(T) ;
        y[2] = sin(P) * sin(T) ;
        z[2] = cos(P);
        points_spheres.push_back(R * x[2]);
        points_spheres.push_back(R * y[2]);
        points_spheres.push_back(R * z[2]);


        normals_spheres.push_back(x[2]);
        normals_spheres.push_back(y[2]);
        normals_spheres.push_back(z[2]);

    }

}
void Scene::draw_cylinder(float R, int prec, std::vector<float> *vertices, std::vector<float> *normals)
{
    vertices->resize(0);
    int rings=360/prec, sectors=360/prec;
    float T, P;
    float x[4],y[4],z[4];
    //Closing nicely the tubes will cause z-fighting and the spherical parts will get all messy

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

        }


}
