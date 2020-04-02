#include "Viewer.h"
#include <CGAL/point_generators_3.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/Exact_spherical_kernel_3.h>
#include <vector>
#include <CGAL/Qt/CreateOpenGLContext.h>


Viewer::Viewer(QWidget* parent )
  : CGAL::QGLViewer(parent)
{
    extension_is_found = false;
}

void Viewer::compile_shaders()
{
    initializeOpenGLFunctions();
    if(! buffers[0].create() || !buffers[1].create() || !buffers[2].create() || !buffers[3].create()
            || !buffers[4].create() || !buffers[4].create() || !buffers[5].create() || !buffers[6].create()
            || !buffers[7].create() || !buffers[8].create())
    {
        std::cerr<<"VBO Creation FAILED"<<std::endl;
    }

    if(!vao[0].create() || !vao[1].create() || !vao[2].create())
    {
        std::cerr<<"VAO Creation FAILED"<<std::endl;
    }

    //The sphere

    //Vertex source code
    const char vertex_source[] =
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
        "   gl_Position = mvp_matrix * (vertex+vec4(center.xyz, 0.0));\n"
        "}"
    };
    //Vertex source code
    const char fragment_source[] =
    {
        "#version 120 \n"
        "varying highp vec4 fP; \n"
        "varying highp vec3 fN; \n"
        "uniform highp vec4 color; \n"
        "uniform highp vec4 light_pos;  \n"
        "uniform highp vec4 light_diff; \n"
        "uniform highp vec4 light_spec; \n"
        "uniform highp vec4 light_amb;  \n"
        "uniform float spec_power ; \n"

        "void main(void) { \n"

        "   highp vec3 L = light_pos.xyz - fP.xyz; \n"
        "   highp vec3 V = -fP.xyz; \n"

        "   highp vec3 N = normalize(fN); \n"
        "   L = normalize(L); \n"
        "   V = normalize(V); \n"

        "   highp vec3 R = reflect(-L, N); \n"
        "   highp vec4 diffuse = max(dot(N,L), 0.0) * light_diff * color; \n"
        "   highp vec4 specular = pow(max(dot(R,V), 0.0), spec_power) * light_spec; \n"

        "gl_FragColor = light_amb*color + diffuse  ; \n"
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
    if(!extension_is_found)
    {


             //Vertex source code
             const char vertex_source_no_ext[] =
             {
                 "#version 120 \n"
                 "attribute highp vec4 vertex;\n"
                 "uniform highp mat4 mvp_matrix;\n"
                 "void main(void)\n"
                 "{\n"
                 "   gl_PointSize = 4.0;\n"
                 "   gl_Position = mvp_matrix * vertex;\n"
                 "}"
             };
             //Vertex source code
             const char fragment_source_no_ext[] =
             {
                 "#version 120 \n"
                 "uniform highp vec4 color; \n"
                 "void main(void) { \n"
                 "gl_FragColor = color; \n"
                 "} \n"
                 "\n"
             };
             vertex_shader = new QOpenGLShader(QOpenGLShader::Vertex);
             if(!vertex_shader->compileSourceCode(vertex_source_no_ext))
             {
                 std::cerr<<"Compiling vertex source FAILED"<<std::endl;
             }

             fragment_shader= new QOpenGLShader(QOpenGLShader::Fragment);
             if(!fragment_shader->compileSourceCode(fragment_source_no_ext))
             {
                 std::cerr<<"Compiling fragmentsource FAILED"<<std::endl;
             }

             if(!rendering_program_no_ext.addShader(vertex_shader))
             {
                 std::cerr<<"adding vertex shader FAILED"<<std::endl;
             }
             if(!rendering_program_no_ext.addShader(fragment_shader))
             {
                 std::cerr<<"adding fragment shader FAILED"<<std::endl;
             }
             if(!rendering_program_no_ext.link())
             {
                 std::cerr<<"linking Program FAILED"<<std::endl;
             }
             rendering_program_no_ext.bind();
    }



}

void Viewer::initialize_buffers()
{
    //The big white sphere
    vao[0].bind();
    //points of the sphere
    buffers[0].bind();
    buffers[0].allocate(pos_sphere.data(),
                        static_cast<int>(pos_sphere.size()*sizeof(float)));
    vertexLocation[0] = rendering_program.attributeLocation("vertex");
    rendering_program.bind();
    rendering_program.enableAttributeArray(vertexLocation[0]);
    rendering_program.setAttributeBuffer(vertexLocation[0],GL_FLOAT,0,3);
    buffers[0].release();
    //normals of the sphere
    buffers[1].bind();
    buffers[1].allocate(normals.data(),
                        static_cast<int>(normals.size()*sizeof(float)));
    normalsLocation[0] = rendering_program.attributeLocation("normal");
    rendering_program.bind();
    rendering_program.enableAttributeArray(normalsLocation[0]);
    rendering_program.setAttributeBuffer(normalsLocation[0],GL_FLOAT,0,3);
    buffers[1].release();
    //center of the sphere
    buffers[2].bind();
    buffers[2].allocate(trivial_center.data(),
                        static_cast<int>(trivial_center.size()*sizeof(float)));
    trivialCenterLocation = rendering_program.attributeLocation("center");
    rendering_program.bind();
    rendering_program.enableAttributeArray(trivialCenterLocation);
    rendering_program.setAttributeBuffer(trivialCenterLocation,GL_FLOAT,0,3);
    buffers[2].release();
    if(extension_is_found)
    {

        glVertexAttribDivisor(trivialCenterLocation, 1);
        glVertexAttribDivisor(normalsLocation[0], 0);
    }
    vao[0].release();

    //The circles
    vao[1].bind();
    buffers[3].bind();
    buffers[3].allocate(pos_lines.data(),
                        static_cast<int>(pos_lines.size()*sizeof(float)));
    vertexLocation[2] = rendering_program.attributeLocation("vertex");
    rendering_program.bind();
    rendering_program.enableAttributeArray(vertexLocation[2]);
    rendering_program.setAttributeBuffer(vertexLocation[2],GL_FLOAT,0,3);
    buffers[3].release();

    //normals
    buffers[4].bind();
    buffers[4].allocate(normals_lines.data(),
                        static_cast<int>(normals_lines.size()*sizeof(float)));
    normalsLocation[1] = rendering_program.attributeLocation("normal");
    rendering_program.bind();
    rendering_program.enableAttributeArray(normalsLocation[1]);
    rendering_program.setAttributeBuffer(normalsLocation[1],GL_FLOAT,0,3);
    buffers[4].release();
    //center
    buffers[5].bind();
    buffers[5].allocate(trivial_center.data(),
                        static_cast<int>(trivial_center.size()*sizeof(float)));
    trivialCenterLocation = rendering_program.attributeLocation("center");
    rendering_program.bind();
    rendering_program.enableAttributeArray(trivialCenterLocation);
    rendering_program.setAttributeBuffer(trivialCenterLocation,GL_FLOAT,0,3);
    buffers[5].release();
    if(extension_is_found)
    {
        glVertexAttribDivisor(trivialCenterLocation, 1);
        glVertexAttribDivisor(normalsLocation[0], 0);
    }
    rendering_program.release();

    vao[1].release();

    //The little green spheres
    vao[2].bind();
    if(extension_is_found)
    {
        //points of the spheres
        buffers[6].bind();
        buffers[6].allocate(pos_sphere_inter.data(),
                            static_cast<int>(pos_sphere_inter.size()*sizeof(float)));
        vertexLocation[2] = rendering_program.attributeLocation("vertex");
        rendering_program.bind();
        rendering_program.enableAttributeArray(vertexLocation[2]);
        rendering_program.setAttributeBuffer(vertexLocation[2],GL_FLOAT,0,3);
        buffers[6].release();
        //normals of the sphere
        buffers[7].bind();
        buffers[7].allocate(normals_inter.data(),
                            static_cast<int>(normals_inter.size()*sizeof(float)));
        normalsLocation[2] = rendering_program.attributeLocation("normal");
        rendering_program.bind();
        rendering_program.enableAttributeArray(normalsLocation[2]);
        rendering_program.setAttributeBuffer(normalsLocation[2],GL_FLOAT,0,3);
        buffers[7].release();
        //center of the sphere
        buffers[8].bind();
        buffers[8].allocate(pos_points.data(),
                            static_cast<int>(pos_points.size()*sizeof(float)));
        centerLocation = rendering_program.attributeLocation("center");
        rendering_program.bind();
        rendering_program.enableAttributeArray(centerLocation);
        rendering_program.setAttributeBuffer(centerLocation,GL_FLOAT,0,3);
        buffers[8].release();

        glVertexAttribDivisor(centerLocation, 1);
        glVertexAttribDivisor(normalsLocation[1], 0);
    }
    else
    {
        //points of the sphere
        buffers[6].bind();
        buffers[6].allocate(pos_points.data(),
                            static_cast<int>(pos_points.size()*sizeof(float)));
        rendering_program_no_ext.bind();
        rendering_program_no_ext.enableAttributeArray("vertex");
        rendering_program_no_ext.setAttributeBuffer("vertex",GL_FLOAT,0,3);
        buffers[6].release();
    }
    vao[2].release();



}

void Viewer::compute_elements()
{

    //The Central Sphere
    {
        pos_sphere.resize(0);
        trivial_center.resize(0);
        int rings=3,sectors=6;
        float T, P, R = 0.999f;
        float x[4],y[4],z[4];


        //Top of the sphere
        for(int t=0; t<360; t+=sectors)
        {

            pos_sphere.push_back(0);
            pos_sphere.push_back(0);
            pos_sphere.push_back(R);


            normals.push_back(0);
            normals.push_back(0);
            normals.push_back(1);



            P = rings*CGAL_PI/180.0;
            T = t*CGAL_PI/180.0;
            x[1] = sin(P) * cos(T) ;
            y[1] = sin(P) * sin(T) ;
            z[1] = cos(P);
            pos_sphere.push_back(R * x[1]);
            pos_sphere.push_back(R * y[1]);
            pos_sphere.push_back(R * z[1]);


            normals.push_back(x[1]);
            normals.push_back(y[1]);
            normals.push_back(z[1]);

            //
            P = rings*CGAL_PI/180.0;
            T = (t+sectors)*CGAL_PI/180.0;
            x[2] = sin(P) * cos(T) ;
            y[2] = sin(P) * sin(T) ;
            z[2] = cos(P);
            pos_sphere.push_back(R * x[2]);
            pos_sphere.push_back(R * y[2]);
            pos_sphere.push_back(R * z[2]);

            normals.push_back(x[2]);
            normals.push_back(y[2]);
            normals.push_back(z[2]);

        }

        //Body of the sphere
        for (int p=rings; p<180-rings; p+=rings)
            for(int t=0; t<360; t+=sectors)
            {
                //A
                P = p*CGAL_PI/180.0;
                T = t*CGAL_PI/180.0;
                x[0] = sin(P) * cos(T) ;
                y[0] = sin(P) * sin(T) ;
                z[0] = cos(P);

                pos_sphere.push_back(R * x[0]);
                pos_sphere.push_back(R * y[0]);
                pos_sphere.push_back(R * z[0]);


                normals.push_back(x[0]);
                normals.push_back(y[0]);
                normals.push_back(z[0]);

                //B
                P = (p+rings)*CGAL_PI/180.0;
                T = t*CGAL_PI/180.0;
                x[1] = sin(P) * cos(T) ;
                y[1] = sin(P) * sin(T) ;
                z[1] = cos(P);
                pos_sphere.push_back(R * x[1]);
                pos_sphere.push_back(R * y[1]);
                pos_sphere.push_back(R * z[1]);


                normals.push_back(x[1]);
                normals.push_back(y[1]);
                normals.push_back(z[1]);

                //C
                P = p*CGAL_PI/180.0;
                T = (t+sectors)*CGAL_PI/180.0;
                x[2] = sin(P) * cos(T) ;
                y[2] = sin(P) * sin(T) ;
                z[2] = cos(P);
                pos_sphere.push_back(R * x[2]);
                pos_sphere.push_back(R * y[2]);
                pos_sphere.push_back(R * z[2]);


                normals.push_back(x[2]);
                normals.push_back(y[2]);
                normals.push_back(z[2]);
                //D
                P = (p+rings)*CGAL_PI/180.0;
                T = (t+sectors)*CGAL_PI/180.0;
                x[3] = sin(P) * cos(T) ;
                y[3] = sin(P) * sin(T) ;
                z[3] = cos(P);
                pos_sphere.push_back(R * x[3]);
                pos_sphere.push_back(R * y[3]);
                pos_sphere.push_back(R * z[3]);


                normals.push_back(x[3]);
                normals.push_back(y[3]);
                normals.push_back(z[3]);



                pos_sphere.push_back(R * x[1]);
                pos_sphere.push_back(R * y[1]);
                pos_sphere.push_back(R * z[1]);


                normals.push_back(x[1]);
                normals.push_back(y[1]);
                normals.push_back(z[1]);

                pos_sphere.push_back(R * x[2]);
                pos_sphere.push_back(R * y[2]);
                pos_sphere.push_back(R * z[2]);


                normals.push_back(x[2]);
                normals.push_back(y[2]);
                normals.push_back(z[2]);

            }
        //Bottom of the sphere
        for(int t=0; t<360; t+=sectors)
        {


            pos_sphere.push_back(0);
            pos_sphere.push_back(0);
            pos_sphere.push_back(-R);


            normals.push_back(0);
            normals.push_back(0);
            normals.push_back(-1);


            P = (180-rings)*CGAL_PI/180.0;
            T = t*CGAL_PI/180.0;
            x[1] = sin(P) * cos(T) ;
            y[1] = sin(P) * sin(T) ;
            z[1] = cos(P);
            pos_sphere.push_back(R * x[1]);
            pos_sphere.push_back(R * y[1]);
            pos_sphere.push_back(R * z[1]);


            normals.push_back(x[1]);
            normals.push_back(y[1]);
            normals.push_back(z[1]);


            P = (180-rings)*CGAL_PI/180.0;
            T = (t+sectors)*CGAL_PI/180.0;
            x[2] = sin(P) * cos(T) ;
            y[2] = sin(P) * sin(T) ;
            z[2] = cos(P);
            pos_sphere.push_back(R * x[2]);
            pos_sphere.push_back(R * y[2]);
            pos_sphere.push_back(R * z[2]);


            normals.push_back(x[2]);
            normals.push_back(y[2]);
            normals.push_back(z[2]);

        }
        trivial_center.push_back(0.0);trivial_center.push_back(0.0);trivial_center.push_back(0.0);
    }
    //The intersection spheres
    {
        pos_sphere_inter.resize(0);
        int rings=3,sectors=3;
        float T, P, R = 0.005f;
        float x[4],y[4],z[4];


        //Top of the sphere
        for(int t=0; t<360; t+=sectors)
        {

            pos_sphere_inter.push_back(0);
            pos_sphere_inter.push_back(0);
            pos_sphere_inter.push_back(R);


            normals_inter.push_back(0);
            normals_inter.push_back(0);
            normals_inter.push_back(1);



            P = rings*CGAL_PI/180.0;
            T = t*CGAL_PI/180.0;
            x[1] = sin(P) * cos(T) ;
            y[1] = sin(P) * sin(T) ;
            z[1] = cos(P);
            pos_sphere_inter.push_back(R * x[1]);
            pos_sphere_inter.push_back(R * y[1]);
            pos_sphere_inter.push_back(R * z[1]);


            normals_inter.push_back(x[1]);
            normals_inter.push_back(y[1]);
            normals_inter.push_back(z[1]);

            //
            P = rings*CGAL_PI/180.0;
            T = (t+sectors)*CGAL_PI/180.0;
            x[2] = sin(P) * cos(T) ;
            y[2] = sin(P) * sin(T) ;
            z[2] = cos(P);
            pos_sphere_inter.push_back(R * x[2]);
            pos_sphere_inter.push_back(R * y[2]);
            pos_sphere_inter.push_back(R * z[2]);

            normals_inter.push_back(x[2]);
            normals_inter.push_back(y[2]);
            normals_inter.push_back(z[2]);

        }

        //Body of the sphere
        for (int p=rings; p<180-rings; p+=rings)
            for(int t=0; t<360; t+=sectors)
            {
                //A
                P = p*CGAL_PI/180.0;
                T = t*CGAL_PI/180.0;
                x[0] = sin(P) * cos(T) ;
                y[0] = sin(P) * sin(T) ;
                z[0] = cos(P);

                pos_sphere_inter.push_back(R * x[0]);
                pos_sphere_inter.push_back(R * y[0]);
                pos_sphere_inter.push_back(R * z[0]);


                normals_inter.push_back(x[0]);
                normals_inter.push_back(y[0]);
                normals_inter.push_back(z[0]);

                //B
                P = (p+rings)*CGAL_PI/180.0;
                T = t*CGAL_PI/180.0;
                x[1] = sin(P) * cos(T) ;
                y[1] = sin(P) * sin(T) ;
                z[1] = cos(P);
                pos_sphere_inter.push_back(R * x[1]);
                pos_sphere_inter.push_back(R * y[1]);
                pos_sphere_inter.push_back(R * z[1]);


                normals_inter.push_back(x[1]);
                normals_inter.push_back(y[1]);
                normals_inter.push_back(z[1]);

                //C
                P = p*CGAL_PI/180.0;
                T = (t+sectors)*CGAL_PI/180.0;
                x[2] = sin(P) * cos(T) ;
                y[2] = sin(P) * sin(T) ;
                z[2] = cos(P);
                pos_sphere_inter.push_back(R * x[2]);
                pos_sphere_inter.push_back(R * y[2]);
                pos_sphere_inter.push_back(R * z[2]);


                normals_inter.push_back(x[2]);
                normals_inter.push_back(y[2]);
                normals_inter.push_back(z[2]);
                //D
                P = (p+rings)*CGAL_PI/180.0;
                T = (t+sectors)*CGAL_PI/180.0;
                x[3] = sin(P) * cos(T) ;
                y[3] = sin(P) * sin(T) ;
                z[3] = cos(P);
                pos_sphere_inter.push_back(R * x[3]);
                pos_sphere_inter.push_back(R * y[3]);
                pos_sphere_inter.push_back(R * z[3]);


                normals_inter.push_back(x[3]);
                normals_inter.push_back(y[3]);
                normals_inter.push_back(z[3]);



                pos_sphere_inter.push_back(R * x[1]);
                pos_sphere_inter.push_back(R * y[1]);
                pos_sphere_inter.push_back(R * z[1]);


                normals_inter.push_back(x[1]);
                normals_inter.push_back(y[1]);
                normals_inter.push_back(z[1]);

                pos_sphere_inter.push_back(R * x[2]);
                pos_sphere_inter.push_back(R * y[2]);
                pos_sphere_inter.push_back(R * z[2]);


                normals_inter.push_back(x[2]);
                normals_inter.push_back(y[2]);
                normals_inter.push_back(z[2]);

            }
        //Bottom of the sphere
        for(int t=0; t<360; t+=sectors)
        {


            pos_sphere_inter.push_back(0);
            pos_sphere_inter.push_back(0);
            pos_sphere_inter.push_back(-R);


            normals_inter.push_back(0);
            normals_inter.push_back(0);
            normals_inter.push_back(-1);


            P = (180-rings)*CGAL_PI/180.0;
            T = t*CGAL_PI/180.0;
            x[1] = sin(P) * cos(T) ;
            y[1] = sin(P) * sin(T) ;
            z[1] = cos(P);
            pos_sphere_inter.push_back(R * x[1]);
            pos_sphere_inter.push_back(R * y[1]);
            pos_sphere_inter.push_back(R * z[1]);


            normals_inter.push_back(x[1]);
            normals_inter.push_back(y[1]);
            normals_inter.push_back(z[1]);


            P = (180-rings)*CGAL_PI/180.0;
            T = (t+sectors)*CGAL_PI/180.0;
            x[2] = sin(P) * cos(T) ;
            y[2] = sin(P) * sin(T) ;
            z[2] = cos(P);
            pos_sphere_inter.push_back(R * x[2]);
            pos_sphere_inter.push_back(R * y[2]);
            pos_sphere_inter.push_back(R * z[2]);


            normals_inter.push_back(x[2]);
            normals_inter.push_back(y[2]);
            normals_inter.push_back(z[2]);

        }
    }


    //init
    {
        pos_points.resize(0);
        pos_lines.resize(0);
        // Restore previous viewer state.
        restoreStateFromFile();

        //random generator of points within a sphere
        typedef CGAL::Creator_uniform_3<EPIC::FT,EPIC::Point_3>   Creator;
        CGAL::Random_points_in_sphere_3<EPIC::Point_3, Creator>   gen;

        const unsigned nb_circles=20;

        //vector to store input points
        std::vector<EPIC::Point_3> points;
        points.reserve(nb_circles);



        for (unsigned i=0;i<nb_circles;++i){
            EPIC::Point_3 p=*++gen;
            //prevent great circles
            while (p.x()==0 && p.y()==0 && p.z()==0) {  p=*++gen; }

            const EPIC::Point_3 origin(0,0,0);
            const EPIC::Plane_3 plane(p, p-origin);
            EPIC::Vector_3 base1=plane.base1();
            EPIC::Vector_3 base2=plane.base2();
            base1=base1/CGAL::sqrt(base1.squared_length());
            base2=base2/CGAL::sqrt(base2.squared_length());
            const double radius=CGAL::sqrt( CGAL::to_double( 1 - CGAL::squared_distance(origin,p) ) );
            const double nb_pt_per_circle=100;
            const double step=2 * CGAL_PI / nb_pt_per_circle;

            for (double theta = 0; theta < 2 * CGAL_PI-step ; theta += step) {
                const EPIC::Point_3 a=p + ( radius*cos(theta)*base1 + radius*sin(theta)*base2 );
                const EPIC::Point_3 b=p + ( radius*cos(theta+step)*base1 + radius*sin(theta+step)*base2 );
                pos_lines.push_back(a.x());pos_lines.push_back(a.y());pos_lines.push_back(a.z());
                normals_lines.push_back(a.x());normals_lines.push_back(a.y());normals_lines.push_back(a.z());
                pos_lines.push_back(b.x());pos_lines.push_back(b.y());pos_lines.push_back(b.z());
                normals_lines.push_back(b.x());normals_lines.push_back(b.y());normals_lines.push_back(b.z());

            }
            points.push_back(p);

        }

        std::vector<EPIC::Point_3> intersections;
        naive_compute_intersection_points(points,std::back_inserter(intersections));

        //draw points as small spheres
        for (std::vector<EPIC::Point_3>::const_iterator it=intersections.begin();it!=intersections.end();++it){
            pos_points.push_back(it->x()); pos_points.push_back(it->y()); pos_points.push_back(it->z());
        }

    }
}

void Viewer::attrib_buffers(CGAL::QGLViewer* viewer)
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
    // define material
    QVector4D        ambient(0.1f, 0.1f, 0.1f, 1.0f);
    QVector4D        diffuse( 0.9f,
                         0.9f,
                         0.9f,
                         0.0f );

    QVector4D        specular(  0.0f,
                           0.0f,
                           0.0f,
                           0.0f );

    QVector4D        position( -1.2f, 1.2f, .9797958971f, 1.0f  );
    GLfloat shininess =  1.0f;


    rendering_program.bind();
    mvpLocation = rendering_program.uniformLocation("mvp_matrix");
    mvLocation = rendering_program.uniformLocation("mv_matrix");
    colorLocation = rendering_program.uniformLocation("color");
    lightLocation[0] = rendering_program.uniformLocation("light_pos");
    lightLocation[1] = rendering_program.uniformLocation("light_diff");
    lightLocation[2] = rendering_program.uniformLocation("light_spec");
    lightLocation[3] = rendering_program.uniformLocation("light_amb");
    lightLocation[4] = rendering_program.uniformLocation("spec_power");

    rendering_program.setUniformValue(lightLocation[0], position);
    rendering_program.setUniformValue(lightLocation[1], diffuse);
    rendering_program.setUniformValue(lightLocation[2], specular);
    rendering_program.setUniformValue(lightLocation[3], ambient);
    rendering_program.setUniformValue(lightLocation[4], shininess);
    rendering_program.setUniformValue(mvpLocation, mvpMatrix);
    rendering_program.setUniformValue(mvLocation, mvMatrix);

    rendering_program.release();
    if(!extension_is_found)
    {
        rendering_program_no_ext.bind();
        rendering_program_no_ext.setUniformValue("mvp_matrix", mvpMatrix);
        rendering_program_no_ext.release();
    }


}

void Viewer::draw()
{
    glEnable(GL_DEPTH_TEST);
    QColor color;

    //sphere
    vao[0].bind();
    attrib_buffers(this);
    rendering_program.bind();
    color.setRgbF(1.0f, 1.0f, 1.0f);
    rendering_program.setUniformValue(colorLocation, color);
    if(extension_is_found)
        glDrawArraysInstanced(GL_TRIANGLES, 0, static_cast<GLsizei>(pos_sphere.size()/3), 1);
    else
        glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(pos_sphere.size()/3));
    rendering_program.release();
    vao[0].release();

    //intersection
    vao[2].bind();
    attrib_buffers(this);
    color.setRgbF(0.0f, 1.0f, 0.0f);
    if(extension_is_found)
    {
        rendering_program.bind();
        rendering_program.setUniformValue(colorLocation, color);
        glDrawArraysInstanced(GL_TRIANGLES, 0, static_cast<GLsizei>(pos_sphere_inter.size()/3), static_cast<GLsizei>(pos_points.size()/3));
        rendering_program.release();
    }
    else
    {
        rendering_program_no_ext.bind();
        rendering_program_no_ext.setUniformValue(colorLocation, color);
        glDrawArrays(GL_POINTS, 0, static_cast<GLsizei>(pos_points.size()/3));
        rendering_program_no_ext.release();
    }

    vao[2].release();

    //circles

    vao[1].bind();
    attrib_buffers(this);
    rendering_program.bind();
    color.setRgbF(1.0f, 0.0f, 0.0f);
    rendering_program.setUniformValue(colorLocation, color);
    glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(pos_lines.size()/3));
    rendering_program.release();
    vao[1].release();


}

void Viewer::init()
{

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
    compile_shaders();
    compute_elements();
    initialize_buffers();
    glEnable(GL_BLEND);

    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_POINT_SMOOTH);
    glEnable(GL_LINE_SMOOTH);


}

template<class Output_iterator>
void Viewer::naive_compute_intersection_points(const std::vector<EPIC::Point_3>& points,Output_iterator out) const {
    typedef CGAL::Exact_spherical_kernel_3 SK;
    SK::Sphere_3 sphere(SK::Point_3(0,0,0),1);

    //converter point to exact SK point type
    CGAL::Cartesian_converter<EPIC,SK> to_exact;
    std::vector<SK::Circle_3> circles;

    //create circles from points: need to use a converter to change floating points coordinates into an exact NT.
    for (std::vector<EPIC::Point_3>::const_iterator it=points.begin();it!=points.end();++it){
        const SK::Point_3 center=to_exact(*it);
        circles.push_back( SK::Circle_3(sphere,SK::Plane_3(center,center-CGAL::ORIGIN) ) );
    }


    //Look for intersection points among pair of circles: use a naive and quadratic way
    for (std::vector<SK::Circle_3>::const_iterator it_f=circles.begin();it_f!=--circles.end();++it_f){
        std::vector<SK::Circle_3>::const_iterator it_s=it_f;
        ++it_s;
        for (;it_s!=circles.end();++it_s){
            std::vector <CGAL::Object> intersections;
            //ensure_circles are different
            CGAL_precondition(*it_s!=*it_f);
            CGAL::intersection(*it_f,*it_s,std::back_inserter(intersections));
            if (!intersections.empty()){
                for (std::vector <CGAL::Object>::const_iterator it_pt=intersections.begin();it_pt!=intersections.end();++it_pt){
                    const std::pair<SK::Circular_arc_point_3,unsigned>* pt=
                            CGAL::object_cast< std::pair<SK::Circular_arc_point_3,unsigned> > (&(*it_pt));
                    assert(pt!=NULL);
                    *out++=EPIC::Point_3( CGAL::to_double(pt->first.x()),
                                          CGAL::to_double(pt->first.y()),
                                          CGAL::to_double(pt->first.z())
                                          );
                }
            }
        }
    }
}
