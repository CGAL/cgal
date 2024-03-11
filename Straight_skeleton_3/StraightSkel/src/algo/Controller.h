/**
 * @file   algo/Controller.h
 * @author Gernot Walzl
 * @date   2012-02-23
 */

#ifndef ALGO_CONTROLLER_H
#define ALGO_CONTROLLER_H

#include "typedefs_thread.h"
#include "algo/ptrs.h"
#include "data/2d/ptrs.h"
#include "data/2d/skel/ptrs.h"
#include "data/3d/ptrs.h"
#include "data/3d/skel/ptrs.h"

namespace algo {

class Controller {
public:
    virtual ~Controller();

    static ControllerSPtr create();

    /**
     * Methods for KeyboardAdapter
     */
    void togglePause();
    void nextStep();
    void skip();
    void haltSkip();

    /**
     * Methods for Algorithms
     */
    void wait();
    void screenshot();

    void setDispPolygon(data::_2d::PolygonSPtr polygon);
    void setDispSkel2d(data::_2d::skel::StraightSkeletonSPtr skel_2d);
    void setDispPolyhedron(data::_3d::PolyhedronSPtr polyhedron);
    void setDispSkel3d(data::_3d::skel::StraightSkeletonSPtr skel_3d);
    void setDispSphericalPolygon(data::_3d::SphericalPolygonSPtr sphericalpolygon);
    void setDispSphericalSkel(data::_3d::skel::SphericalSkeletonSPtr sphericalskel);

    /**
     * Methods for OpenGLWindow
     */
    data::_2d::PolygonSPtr getDispPolygon();
    data::_2d::skel::StraightSkeletonSPtr getDispSkel2d();
    data::_3d::PolyhedronSPtr getDispPolyhedron();
    data::_3d::skel::StraightSkeletonSPtr getDispSkel3d();
    data::_3d::SphericalPolygonSPtr getDispSphericalPolygon();
    data::_3d::skel::SphericalSkeletonSPtr getDispSphericalSkel();

    bool getScreenshot();

protected:
    Controller();

    mutable RecursiveMutex mutex_;

    bool pause_;
    bool next_step_;
    unsigned int time_sleep_ms_;
    unsigned int time_sleep_ms_before_;

    bool screenshot_;
    bool screenshot_on_wait_;

    data::_2d::PolygonSPtr polygon_;
    data::_2d::skel::StraightSkeletonSPtr skel_2d_;
    data::_3d::PolyhedronSPtr polyhedron_;
    data::_3d::skel::StraightSkeletonSPtr skel_3d_;
    data::_3d::SphericalPolygonSPtr sphericalpolygon_;
    data::_3d::skel::SphericalSkeletonSPtr sphericalskel_;
};

}

#endif /* ALGO_CONTROLLER_H */
