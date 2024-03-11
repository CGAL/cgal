/**
 * @file   ui/gl/OpenGLWindow.cpp
 * @author Gernot Walzl
 * @date   2011-12-19
 */

#include "ui/gl/OpenGLWindow.h"

#include "debug.h"
#include "ui/gl/Camera.h"
#include "ui/vecmath.h"
#include <cmath>
#include <stdexcept>

namespace ui { namespace gl {

std::map<int, OpenGLWindowWPtr> OpenGLWindow::windows_;
bool OpenGLWindow::glut_started_ = false;
RecursiveMutex OpenGLWindow::mutex_;


OpenGLWindow::OpenGLWindow(int argc, const char* argv[],
        unsigned int width, unsigned int height,
        const char* title) {
    this->argc_ = argc;
    this->argv_ = argv;
    this->title_ = title;
    this->width_ = width;
    this->height_ = height;
    this->fovy_ = 60.0;
    this->camera_ = Camera::create();
    this->window_ = -1;
    for (unsigned int i = 0; i < 4; i++) {
        this->color_[i] = 0.0f;
    }
}


OpenGLWindow::~OpenGLWindow() {
    Lock l(mutex_);
    this->window_ = -1;
    this->camera_.reset();
}


OpenGLWindowSPtr OpenGLWindow::create(int argc, const char* argv[],
        unsigned int width, unsigned int height, const char* title) {
    return OpenGLWindowSPtr(new OpenGLWindow(
            argc, argv, width, height, title));
}


void OpenGLWindow::createWindow() {
    Lock l(mutex_);
    glutInit(&argc_, const_cast<char**>(argv_));
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
    glutInitWindowSize(width_, height_);
    glutInitWindowPosition(0, 0);
    this->window_ = glutCreateWindow(title_);
    windows_[window_] = shared_from_this();
    glutDisplayFunc(OpenGLWindow::displayFunc);
    glutIdleFunc(OpenGLWindow::displayFunc);
    glutReshapeFunc(OpenGLWindow::reshapeFunc);
    glutKeyboardFunc(OpenGLWindow::keyFunc);
    glutKeyboardUpFunc(OpenGLWindow::keyUpFunc);
    glutSpecialFunc(OpenGLWindow::specialFunc);
    glutSpecialUpFunc(OpenGLWindow::specialUpFunc);
    glutMouseFunc(OpenGLWindow::mouseFunc);
    glutMotionFunc(OpenGLWindow::motionFunc);
    this->init();
}


void OpenGLWindow::initLighting() {
    Lock l(mutex_);
    GLfloat light_color_am[4] = {0.2f, 0.2f, 0.2f, 1.0f};
    GLfloat light_color_diff[4] = {0.8f, 0.8f, 0.8f, 1.0f};
    GLfloat light_color_spec[4] = {0.8f, 0.8f, 0.8f, 1.0f};
    glShadeModel(GL_SMOOTH);
    glLightfv(GL_LIGHT0, GL_AMBIENT, light_color_am);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, light_color_diff);
    glLightfv(GL_LIGHT0, GL_SPECULAR, light_color_spec);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);

    glEnable(GL_COLOR_MATERIAL);
    GLfloat mat_spec[4] = {0.2f, 0.2f, 0.2f, 1.0f};
    GLint mat_shini = 8;
    glMaterialfv(GL_FRONT, GL_SPECULAR, mat_spec);
    glMateriali(GL_FRONT, GL_SHININESS, mat_shini);

    glEnable(GL_NORMALIZE);
}

void OpenGLWindow::initBlending() {
    Lock l(mutex_);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
}

void OpenGLWindow::init() {
    Lock l(mutex_);
    glutSetWindow(window_);
    glClearColor(0.0, 0.0, 0.0, 0.0);
    glClearDepth(1.0);
    glDepthFunc(GL_LESS);
    glEnable(GL_DEPTH_TEST);
    this->initLighting();
    this->initBlending();
    glEnable(GL_CULL_FACE);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(fovy_, (double)width_/(double)height_, 0.001, 1000.0);
    glMatrixMode(GL_MODELVIEW);
}


void OpenGLWindow::setColor(const vec4f rgba) {
    glColor4f(rgba[0], rgba[1], rgba[2], rgba[3]);
    for (unsigned int i = 0; i < 4; i++) {
        this->color_[i] = rgba[i];
    }
}

void OpenGLWindow::getColor(vec4f& out) const {
    for (unsigned int i = 0; i < 4; i++) {
        out[i] = this->color_[i];
    }
}

void OpenGLWindow::drawSphere(const vec3f position, float radius) {
    glPushMatrix();
    glTranslatef(position[0], position[1], position[2]);
    GLUquadric* quad = gluNewQuadric();
    gluQuadricDrawStyle(quad, GLU_FILL);
    gluQuadricNormals(quad, GLU_SMOOTH);
    if (radius > 0.5f) {
        gluSphere(quad, radius, 64, 64);
    } else {
        gluSphere(quad, radius, 8, 8);
    }
    gluDeleteQuadric(quad);
    glPopMatrix();
}

void OpenGLWindow::drawCylinder(const vec3f src, const vec3f dst, float radius) {
    float dx = dst[0] - src[0];
    float dy = dst[1] - src[1];
    float dz = dst[2] - src[2];
    float length = sqrtf(dx*dx + dy*dy + dz*dz);

    glPushMatrix();
    glTranslatef(src[0], src[1], src[2]);

    float angle_xy = atan2f(dy, dx);
    float factor_x = -sinf(angle_xy);
    float factor_y = cosf(angle_xy);
    float angle_z = acosf(dz/length);
    glRotatef(angle_z * (180.0/M_PI), factor_x, factor_y, 0.0f);

    GLUquadric* quad = gluNewQuadric();
    gluQuadricDrawStyle(quad, GLU_FILL);
    gluQuadricNormals(quad, GLU_SMOOTH);
    gluCylinder(quad, radius, radius, length, 8, 1);
    gluDeleteQuadric(quad);
    glPopMatrix();
}

void OpenGLWindow::drawCircularCylinder(const vec3f center, const vec3f axis,
        const vec3f src, const vec3f dst, float radius) {
    vec3f dir_src;
    vec3f dir_dst;
    for (unsigned int i = 0; i < 3; i++) {
        dir_src[i] = src[i] - center[i];
        dir_dst[i] = dst[i] - center[i];
    }
    float radius_sphere = 0.0f;
    for (unsigned int i = 0; i < 3; i++) {
        radius_sphere += dir_src[i] * dir_src[i];
    }
    radius_sphere = sqrtf(radius_sphere);

    vec3f n_axis;
    copy(axis, n_axis);
    normalize(n_axis);
    float angle_axis = angle(axis, dir_src);
    float radius_circle = radius_sphere * cosf(M_PI/2.0f - angle_axis);
    float dist_center_circle = radius_sphere * sinf(M_PI/2.0f - angle_axis);
    vec3f center_circle;
    copy(center, center_circle);
    for (unsigned int i = 0; i < 3; i++) {
        center_circle[i] += dist_center_circle * n_axis[i];
    }

    vec3f dir_axis_src;
    vec3f dir_axis_dst;
    vec3f dir_tmp;
    cross(dir_src, axis, dir_tmp);
    cross(axis, dir_tmp, dir_axis_src);
    normalize(dir_axis_src);
    cross(dir_dst, axis, dir_tmp);
    cross(axis, dir_tmp, dir_axis_dst);
    normalize(dir_axis_dst);

    unsigned int num_cylinders = 32;  // has to be a number with base 2
    vec3f dir[num_cylinders + 1];
    copy(dir_axis_src, dir[0]);
    copy(dir_axis_dst, dir[num_cylinders]);
    unsigned int inc = num_cylinders/2;
    while (inc > 0) {
        for (unsigned int i = inc; i < num_cylinders; i += 2*inc) {
            for (unsigned int j = 0; j < 3; j++) {
                dir[i][j] = dir[i-inc][j] + dir[i+inc][j];
            }
            if (i == num_cylinders/2) {
                if (length(dir[i]) < 0.0005) {
                    // angle between src and dst = M_PI
                    cross(axis, dir[0], dir[i]);
                } else {
                    // fix orientation when angle between src and dst > M_PI
                    vec3f dir_normal;
                    cross(dir[0], dir[i], dir_normal);
                    if (angle(dir_normal, axis) > M_PI/2.0f) {
                        scale(dir[i], -1.0f);
                    }
                }
            }
            normalize(dir[i]);
        }
        inc /= 2;
    }
    for (unsigned int i = 0; i <= num_cylinders; i++) {
        scale(dir[i], radius_circle);
    }

    glPushMatrix();
    glTranslatef(center_circle[0], center_circle[1], center_circle[2]);
    for (unsigned int i = 0; i < num_cylinders; i++) {
        drawCylinder(dir[i], dir[i+1], radius);
    }
    glPopMatrix();
}

void OpenGLWindow::drawTriangle(const vec3f a, const vec3f b, const vec3f c) {
    glBegin(GL_TRIANGLES);
    glVertex3f(a[0], a[1], a[2]);
    glVertex3f(b[0], b[1], b[2]);
    glVertex3f(c[0], c[1], c[2]);
    glEnd();
}

void OpenGLWindow::drawArrow(const vec3f position, const vec3f direction) {
    float length = sqrtf(direction[0]*direction[0] +
                         direction[1]*direction[1] +
                         direction[2]*direction[2]);
    glPushMatrix();
    glTranslatef(position[0], position[1], position[2]);

    float angle_xy = atan2f(direction[1], direction[0]);
    float factor_x = -sinf(angle_xy);
    float factor_y = cosf(angle_xy);
    float angle_z = acosf(direction[2]/length);
    glRotatef(angle_z * (180.0/M_PI), factor_x, factor_y, 0.0f);

    GLUquadric* quad = gluNewQuadric();
    gluQuadricDrawStyle(quad, GLU_FILL);
    gluQuadricNormals(quad, GLU_SMOOTH);
    gluCylinder(quad, length/12.0, length/12.0, length*(2.0/3.0), 8, 1);
    gluDeleteQuadric(quad);

    glTranslatef(0.0, 0.0, length*(2.0/3.0));

    quad = gluNewQuadric();
    gluQuadricDrawStyle(quad, GLU_FILL);
    gluQuadricNormals(quad, GLU_SMOOTH);
    gluCylinder(quad, length/6.0, 0.0, length/3.0, 8, 1);
    gluDeleteQuadric(quad);

    glPopMatrix();
}

void OpenGLWindow::drawText(const std::string& text) {
    glDisable(GL_LIGHTING);
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    gluOrtho2D(0, width_, 0, height_);
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    glTranslatef(10.0, 10.0, 0.0);
    glScalef(0.5, 0.5, 1.0);
    glLineWidth(2.0);
    glColor3f(1.0, 1.0, 1.0);
    const char* c_str = text.c_str();
    for (unsigned int i = 0; i < text.length(); i++) {
        glutStrokeCharacter(GLUT_STROKE_ROMAN, c_str[i]);
    }
    glPopMatrix();
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glEnable(GL_LIGHTING);
}

void OpenGLWindow::drawCrosshair(float size) {
    glDisable(GL_LIGHTING);
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    gluOrtho2D(0, width_, 0, height_);
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    glTranslatef(width_/2.0f, height_/2.0f, 0.0f);
    glLineWidth(2.0);
    glColor3f(1.0, 1.0, 1.0);
    glBegin(GL_LINES);
    glVertex2f(-size, 0.0f);
    glVertex2f(+size, 0.0f);
    glEnd();
    glBegin(GL_LINES);
    glVertex2f(0.0f, -size);
    glVertex2f(0.0f, +size);
    glEnd();
    glPopMatrix();
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glEnable(GL_LIGHTING);
}

void OpenGLWindow::drawContent() {
    throw std::runtime_error("drawContent()"
            " has to be implemented by inherited class.");
}

void OpenGLWindow::handleKeyPressed(int key, int x, int y) {
    throw std::runtime_error("handleKeyPressed(...)"
            " has to be implemented by inherited class.");
}

void OpenGLWindow::handleKeyReleased(int key, int x, int y) {
    throw std::runtime_error("handleKeyReleased(...)"
            " has to be implemented by inherited class.");
}

void OpenGLWindow::handleMousePressed(int button, int x, int y) {
    throw std::runtime_error("handleMousePressed(...)"
            " has to be implemented by inherited class.");
}

void OpenGLWindow::handleMouseReleased(int button, int x, int y) {
    throw std::runtime_error("handleMouseReleased(...)"
            " has to be implemented by inherited class.");
}

void OpenGLWindow::handleMotion(int x, int y) {
    throw std::runtime_error("handleMotion(...)"
            " has to be implemented by inherited class.");
}


void OpenGLWindow::displayFunc() {
    {
        Lock l(mutex_);  // scoped lock
        OpenGLWindowSPtr window = windows_[glutGetWindow()].lock();
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glLoadIdentity();
        CameraSPtr cam = window->camera_;
        gluLookAt(cam->eye(0), cam->eye(1), cam->eye(2),
                  cam->center(0), cam->center(1), cam->center(2),
                  cam->up(0), cam->up(1), cam->up(2));
        //GLfloat light_pos[4] = {10.0f, 10.0f, 10.0f, 1.0f};
        //glLightfv(GL_LIGHT0, GL_POSITION, light_pos);
        window->drawContent();
        glutSwapBuffers();
        // TODO: for each window: glutSetWindow(); glutPostRedisplay();
    }
    thread_sleep(10);
}


void OpenGLWindow::reshapeFunc(int width, int height) {
    Lock l(mutex_);
    if (height <= 0) {
        height = 1;
    }
    OpenGLWindowSPtr window = windows_[glutGetWindow()].lock();
    window->width_ = width;
    window->height_ = height;
    glViewport(0, 0, width, height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(window->fovy_, (double)width/(double)height, 0.01, 100.0);
    glMatrixMode(GL_MODELVIEW);
}


void OpenGLWindow::keyFunc(unsigned char key, int x, int y) {
    Lock l(mutex_);
    OpenGLWindowSPtr window = windows_[glutGetWindow()].lock();
    window->handleKeyPressed(key, x, y);
}

void OpenGLWindow::keyUpFunc(unsigned char key, int x, int y) {
    Lock l(mutex_);
    OpenGLWindowSPtr window = windows_[glutGetWindow()].lock();
    window->handleKeyReleased(key, x, y);
}

void OpenGLWindow::specialFunc(int key, int x, int y) {
    Lock l(mutex_);
    OpenGLWindowSPtr window = windows_[glutGetWindow()].lock();
    // define negative values for special keys
    window->handleKeyPressed(-key, x, y);
}

void OpenGLWindow::specialUpFunc(int key, int x, int y) {
    Lock l(mutex_);
    OpenGLWindowSPtr window = windows_[glutGetWindow()].lock();
    window->handleKeyReleased(-key, x, y);
}

void OpenGLWindow::mouseFunc(int button, int state, int x, int y) {
    Lock l(mutex_);
    OpenGLWindowSPtr window = windows_[glutGetWindow()].lock();
    if (state == GLUT_DOWN) {
        window->handleMousePressed(button, x, y);
    } else if (state == GLUT_UP) {
        window->handleMouseReleased(button, x, y);
    }
}

void OpenGLWindow::motionFunc(int x, int y) {
    Lock l(mutex_);
    OpenGLWindowSPtr window = windows_[glutGetWindow()].lock();
    window->handleMotion(x, y);
}


void OpenGLWindow::mainLoop() {
    if (!glut_started_) {
        glut_started_ = true;
        glutMainLoop();
    }
}


bool OpenGLWindow::writeBMP(const char* filename, unsigned int width, unsigned int height, unsigned char* image) {
    // http://de.wikipedia.org/wiki/Windows_Bitmap
    if (width % 4) {
        DEBUG_VAL("Warning: width should be a multiple of 4.");
        DEBUG_VAR(width);
    }
    bool result = false;
    FILE* fptr = fopen(filename, "w");
    if (fptr) {
        unsigned char bmfh[54] = {
                0x42, 0x4d,   //  [0] bfType (="BM")
                54, 0, 0, 0,  //  [2] bfSize
                0, 0, 0, 0,   //  [6] bfReserved (=0)
                54, 0, 0, 0,  // [10] bfOffBytes (=54)
                40, 0, 0, 0,  // [14] biSize (=40)
                0, 0, 0, 0,   // [18] biWidth
                0, 0, 0, 0,   // [22] biHeight
                1, 0,         // [26] biPlanes (=1)
                24, 0,        // [28] biBitCount (=24)
                0, 0, 0, 0,   // [30] biCompression (=0)
                0, 0, 0, 0,   // [34] biSizeImage
                0, 0, 0, 0,   // [38] biXPelsPerMeter (=0)
                0, 0, 0, 0,   // [42] biYPelsPerMeter (=0)
                0, 0, 0, 0,   // [46] biClrUsed (=0)
                0, 0, 0, 0    // [50] biClrImportant (=0)
        };
        unsigned int size_image = 3 * width * height;
        unsigned int size_file = size_image + 54;
        //bmfh.bfSize = size_file;
        bmfh[2] = size_file & 0xFF;
        bmfh[3] = (size_file >> 8) & 0xFF;
        bmfh[4] = (size_file >> 16) & 0xFF;
        bmfh[5] = (size_file >> 24) & 0xFF;
        //bmfh.biSizeImage = size_image
        bmfh[34] = size_image & 0xFF;
        bmfh[35] = (size_image >> 8) & 0xFF;
        bmfh[36] = (size_image >> 16) & 0xFF;
        bmfh[37] = (size_image >> 24) & 0xFF;
        //bmfh.biWidth = width;
        bmfh[18] = width & 0xFF;
        bmfh[19] = (width >> 8) & 0xFF;
        bmfh[20] = (width >> 16) & 0xFF;
        bmfh[21] = (width >> 24) & 0xFF;
        //bmfh.biHeight = height;
        bmfh[22] = height & 0xFF;
        bmfh[23] = (height >> 8) & 0xFF;
        bmfh[24] = (height >> 16) & 0xFF;
        bmfh[25] = (height >> 24) & 0xFF;

        fwrite((void*)bmfh, 54, 1, fptr);

        unsigned char* frow = new unsigned char[3 * width];
        for (unsigned int h = 0; h < height; h++){
            for (unsigned int w = 0; w < width; w++) {
                for (unsigned int i = 0; i < 3; i++) {
                    frow[3*w + i] = image[(3*width*h + 3*w) + (2-i)];
                }
            }
            fwrite((void*)frow, 3*width, 1, fptr);
        }
        delete[] frow;

        fclose(fptr);
        result = true;
    }
    return result;
}

bool OpenGLWindow::dumpWindow(const char* filename) {
    bool result = false;
    Lock l(mutex_);
    glutSetWindow(window_);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();
    gluLookAt(camera_->eye(0), camera_->eye(1), camera_->eye(2),
              camera_->center(0), camera_->center(1), camera_->center(2),
              camera_->up(0), camera_->up(1), camera_->up(2));
    drawContent();
    glutSwapBuffers();
    unsigned char* image = new unsigned char[3 * width_ * height_];
    glReadPixels(0, 0, width_, height_, GL_RGB, GL_UNSIGNED_BYTE, image);
    result = writeBMP(filename, width_, height_, image);
    delete[] image;
    return result;
}

} }
