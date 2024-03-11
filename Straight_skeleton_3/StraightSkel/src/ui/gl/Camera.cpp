/**
 * @file   ui/gl/Camera.cpp
 * @author Gernot Walzl
 * @date   2011-12-19
 */

#include "ui/gl/Camera.h"

#include "util/StringFactory.h"
#include "util/Configuration.h"
#include <cmath>
#include <stdexcept>

namespace ui { namespace gl {

Camera::Camera() {
    reset();
}

Camera::~Camera() {
    // intentionally does nothing
}

CameraSPtr Camera::create() {
    CameraSPtr result = CameraSPtr(new Camera());
    return result;
}


void Camera::loadConfig(util::ConfigurationSPtr config) {
    if (!config->isLoaded()) {
        return;
    }
    std::string section("ui_gl_Camera");
    eye_[0] = config->getDouble(section, "eye_x");
    eye_[1] = config->getDouble(section, "eye_y");
    eye_[2] = config->getDouble(section, "eye_z");
    center_[0] = config->getDouble(section, "center_x");
    center_[1] = config->getDouble(section, "center_y");
    center_[2] = config->getDouble(section, "center_z");
}


void Camera::reset() {
    this->eye_[0] = 10.0;
    this->eye_[1] = 10.0;
    this->eye_[2] = 10.0;
    this->center_[0] = 0.0;
    this->center_[1] = 0.0;
    this->center_[2] = 0.0;
    this->up_[0] = 0.0;
    this->up_[1] = 0.0;
    this->up_[2] = 1.0;
    this->loadConfig(util::Configuration::getInstance());
}

void Camera::topdown() {
    this->eye_[0] = 0.0;
    this->eye_[1] = -0.001;  // fixes problems with up vector
    this->eye_[2] = 15.0;
    this->center_[0] = 0.0;
    this->center_[1] = 0.0;
    this->center_[2] = 0.0;
}


std::string Camera::toString() {
    std::string result;
    result += "eye=" + util::StringFactory::fromDoubleArr(3, eye_) + "\n";
    result += "center=" + util::StringFactory::fromDoubleArr(3, center_) + "\n";
    result += "up=" + util::StringFactory::fromDoubleArr(3, up_) + "\n";
    return result;
}


GLdouble Camera::eye(unsigned int i) {
    if (i > 2) {
        throw std::out_of_range("");
    }
    return eye_[i];
}

GLdouble Camera::center(unsigned int i) {
    if (i > 2) {
        throw std::out_of_range("");
    }
    return center_[i];
}

GLdouble Camera::up(unsigned int i) {
    if (i > 2) {
        throw std::out_of_range("");
    }
    return up_[i];
}


double Camera::angleXY() {
    double dx = center_[0] - eye_[0];
    double dy = center_[1] - eye_[1];
    return atan2(dy, dx);
}

double Camera::angleZ() {
    double dx = center_[0] - eye_[0];
    double dy = center_[1] - eye_[1];
    double dz = center_[2] - eye_[2];
    double length_xy = sqrt(dx*dx + dy*dy);
    return atan2(dz, length_xy);
}


void Camera::move(double dx, double dy, double dz) {
    eye_[0] += dx;
    center_[0] += dx;
    eye_[1] += dy;
    center_[1] += dy;
    eye_[2] += dz;
    center_[2] += dz;
}

void Camera::moveFB(double amount) {
    double angle_xy = this->angleXY();
    double angle_z = this->angleZ();
    double dx = amount * cos(angle_xy) * cos(angle_z);
    double dy = amount * sin(angle_xy) * cos(angle_z);
    double dz = amount * sin(angle_z);
    this->move(dx, dy, dz);
}

void Camera::moveUD(double amount) {
    this->move(0.0, 0.0, amount);
}

void Camera::strafeLR(double amount) {
    double angle_xy = this->angleXY();
    angle_xy -= M_PI/2.0;
    double dx = amount * cos(angle_xy);
    double dy = amount * sin(angle_xy);
    this->move(dx, dy, 0.0);
}

void Camera::strafeUD(double amount) {
    double angle_xy = this->angleXY();
    double angle_z = this->angleZ();
    double dx = -amount * cos(angle_xy) * sin(angle_z);
    double dy = -amount * sin(angle_xy) * sin(angle_z);
    double dz = amount * cos(angle_z);
    this->move(dx, dy, dz);
}


void Camera::lookLR(double amount) {
    double angle_xy = this->angleXY();
    angle_xy -= amount;
    double old_dx = center_[0] - eye_[0];
    double old_dy = center_[1] - eye_[1];
    double length_xy = sqrt(old_dx*old_dx + old_dy*old_dy);
    double new_dx = length_xy * cos(angle_xy);
    double new_dy = length_xy * sin(angle_xy);
    center_[0] = eye_[0] + new_dx;
    center_[1] = eye_[1] + new_dy;
}

void Camera::lookUD(double amount) {
    double angle_z = this->angleZ();
    angle_z += amount;
    if (-M_PI/2.0 >= angle_z or angle_z >= M_PI/2.0) {
        return;
    }
    double angle_xy = this->angleXY();
    double old_dx = center_[0] - eye_[0];
    double old_dy = center_[1] - eye_[1];
    double old_dz = center_[2] - eye_[2];
    double length = sqrt(old_dx*old_dx + old_dy*old_dy + old_dz*old_dz);
    double new_length_xy = length * cos(angle_z);
    double new_dx = new_length_xy * cos(angle_xy);
    double new_dy = new_length_xy * sin(angle_xy);
    double new_dz = length * sin(angle_z);
    center_[0] = eye_[0] + new_dx;
    center_[1] = eye_[1] + new_dy;
    center_[2] = eye_[2] + new_dz;
}


void Camera::rotateWorldLR(double amount) {
    double angle = atan2(eye_[1], eye_[0]);
    double length = sqrt(eye_[0]*eye_[0] + eye_[1]*eye_[1]);
    angle -= amount;
    eye_[0] = length * cos(angle);
    eye_[1] = length * sin(angle);
    angle = atan2(center_[1], center_[0]);
    length = sqrt(center_[0]*center_[0] + center_[1]*center_[1]);
    angle -= amount;
    center_[0] = length * cos(angle);
    center_[1] = length * sin(angle);
}

void Camera::rotateWorldUD(double amount) {
    double angle_z = this->angleZ();
    angle_z += 2*amount;
    if (-M_PI/2.0 >= angle_z or angle_z >= M_PI/2.0) {
        return;
    }
    double angle_xy = atan2(eye_[1], eye_[0]);

    // rotate x axis
    double factor = sin(angle_xy);

    double angle = atan2(eye_[2], eye_[1]);
    double length = sqrt(eye_[2]*eye_[2] + eye_[1]*eye_[1]);
    angle -= amount * factor;
    eye_[1] = length * cos(angle);
    eye_[2] = length * sin(angle);

    angle = atan2(center_[2], center_[1]);
    length = sqrt(center_[2]*center_[2] + center_[1]*center_[1]);
    angle -= amount * factor;
    center_[1] = length * cos(angle);
    center_[2] = length * sin(angle);

    // rotate y axis
    factor = cos(angle_xy);

    angle = atan2(eye_[2], eye_[0]);
    length = sqrt(eye_[2]*eye_[2] + eye_[0]*eye_[0]);
    angle -= amount * factor;
    eye_[0] = length * cos(angle);
    eye_[2] = length * sin(angle);

    angle = atan2(center_[2], center_[0]);
    length = sqrt(center_[2]*center_[2] + center_[0]*center_[0]);
    angle -= amount * factor;
    center_[0] = length * cos(angle);
    center_[2] = length * sin(angle);
}

} }
