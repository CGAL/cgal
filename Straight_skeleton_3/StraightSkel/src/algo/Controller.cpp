/**
 * @file   algo/Controller.cpp
 * @author Gernot Walzl
 * @date   2012-02-23
 */

#include "algo/Controller.h"

#include "util/Configuration.h"
#include <string>

namespace algo {

Controller::Controller() {
    pause_ = false;
    next_step_ = false;
    time_sleep_ms_ = 3000;
    time_sleep_ms_before_ = time_sleep_ms_;
    screenshot_ = false;
    screenshot_on_wait_ = false;
    util::ConfigurationSPtr config = util::Configuration::getInstance();
    if (config->isLoaded()) {
        std::string section("algo_Controller");
        double sleep = config->getDouble(section, "time_sleep");
        if (sleep != 0.0) {
            time_sleep_ms_ = (int)(sleep * 1000);
        }
        screenshot_on_wait_ = config->getBool(section, "screenshot_on_wait");
    }
}

Controller::~Controller() {
}

ControllerSPtr Controller::create() {
    return ControllerSPtr(new Controller());
}

void Controller::togglePause() {
    Lock l(mutex_);
    pause_ = !pause_;
}

void Controller::nextStep() {
    Lock l(mutex_);
    next_step_ = true;
}

void Controller::skip() {
    Lock l(mutex_);
    time_sleep_ms_before_ = time_sleep_ms_;
    time_sleep_ms_ = 0;
    pause_ = false;
}

void Controller::haltSkip() {
    Lock l(mutex_);
    time_sleep_ms_ = time_sleep_ms_before_;
}

void Controller::wait() {
    if (screenshot_on_wait_) {
        screenshot();
    }
    unsigned int sleep_ms = 10;
    unsigned int current = 0;
    while (pause_ || current <= time_sleep_ms_) {
        if (next_step_) {
            next_step_ = false;
            break;
        }
        thread_sleep(sleep_ms);
        current += sleep_ms;
    }
}

void Controller::screenshot() {
    Lock l(mutex_);
    screenshot_ = true;
}

void Controller::setDispPolygon(data::_2d::PolygonSPtr polygon) {
    Lock l(mutex_);
    polygon_ = polygon;
}

void Controller::setDispSkel2d(data::_2d::skel::StraightSkeletonSPtr skel_2d) {
    Lock l(mutex_);
    skel_2d_ = skel_2d;
}

void Controller::setDispPolyhedron(data::_3d::PolyhedronSPtr polyhedron) {
    Lock l(mutex_);
    polyhedron_ = polyhedron;
}

void Controller::setDispSkel3d(data::_3d::skel::StraightSkeletonSPtr skel_3d) {
    Lock l(mutex_);
    skel_3d_ = skel_3d;
}

void Controller::setDispSphericalPolygon(data::_3d::SphericalPolygonSPtr sphericalpolygon) {
    Lock l(mutex_);
    sphericalpolygon_ = sphericalpolygon;
}

void Controller::setDispSphericalSkel(data::_3d::skel::SphericalSkeletonSPtr sphericalskel) {
    Lock l(mutex_);
    sphericalskel_ = sphericalskel;
}

data::_2d::PolygonSPtr Controller::getDispPolygon() {
    Lock l(mutex_);
    data::_2d::PolygonSPtr result = polygon_;
    polygon_ = data::_2d::PolygonSPtr();
    return result;
}

data::_2d::skel::StraightSkeletonSPtr Controller::getDispSkel2d() {
    Lock l(mutex_);
    data::_2d::skel::StraightSkeletonSPtr result = skel_2d_;
    skel_2d_ = data::_2d::skel::StraightSkeletonSPtr();
    return result;
}

data::_3d::PolyhedronSPtr Controller::getDispPolyhedron() {
    Lock l(mutex_);
    data::_3d::PolyhedronSPtr result = polyhedron_;
    polyhedron_ = data::_3d::PolyhedronSPtr();
    return result;
}

data::_3d::skel::StraightSkeletonSPtr Controller::getDispSkel3d() {
    Lock l(mutex_);
    data::_3d::skel::StraightSkeletonSPtr result = skel_3d_;
    skel_3d_ = data::_3d::skel::StraightSkeletonSPtr();
    return result;
}

data::_3d::SphericalPolygonSPtr Controller::getDispSphericalPolygon() {
    Lock l(mutex_);
    data::_3d::SphericalPolygonSPtr result = sphericalpolygon_;
    sphericalpolygon_ = data::_3d::SphericalPolygonSPtr();
    return result;
}

data::_3d::skel::SphericalSkeletonSPtr Controller::getDispSphericalSkel() {
    Lock l(mutex_);
    data::_3d::skel::SphericalSkeletonSPtr result = sphericalskel_;
    sphericalskel_ = data::_3d::skel::SphericalSkeletonSPtr();
    return result;
}

bool Controller::getScreenshot() {
    Lock l(mutex_);
    bool result = screenshot_;
    screenshot_ = false;
    return result;
}

}
