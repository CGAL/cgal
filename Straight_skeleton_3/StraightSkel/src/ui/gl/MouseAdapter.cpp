/**
 * @file   ui/gl/MouseAdapter.cpp
 * @author Gernot Walzl
 * @date   2011-12-22
 */

#include "ui/gl/MouseAdapter.h"

#include "util/Configuration.h"
#include "ui/gl/Camera.h"

namespace ui { namespace gl {

MouseAdapter::MouseAdapter(CameraSPtr camera) {
    this->camera_ = camera;

    this->b_look_ = 1;
    this->b_rotate_ = 0;
    this->b_pan_ = 2;
    this->b_forward_ = 3;
    this->b_backpedal_ = 4;

    this->mouse_look_ = false;
    this->mouse_rotate_ = false;
    this->mouse_pan_ = false;

    this->sensitivity_look_ = 0.01;
    this->sensitivity_rotate_ = 0.01;
    this->sensitivity_pan_ = 0.01;
    this->sensitivity_move_ = 1.0;

    loadConfig(util::Configuration::getInstance());

    this->prev_pos_[0] = 0;
    this->prev_pos_[1] = 0;
}

MouseAdapter::~MouseAdapter() {
    this->camera_.reset();
}

MouseAdapterSPtr MouseAdapter::create(CameraSPtr camera) {
    MouseAdapterSPtr result = MouseAdapterSPtr(new MouseAdapter(camera));
    return result;
}

int MouseAdapter::toButton(const std::string& button) {
    int result = -1;
    if (button.compare("MOUSE1") == 0) {
        result = 0;
    } else if (button.compare("MOUSE2") == 0) {
        result = 1;
    } else if (button.compare("MOUSE3") == 0) {
        result = 2;
    } else if (button.compare("MWHEELUP") == 0) {
        result = 3;
    } else if (button.compare("MWHEELDOWN") == 0) {
        result = 4;
    }
    return result;
}

void MouseAdapter::loadConfig(util::ConfigurationSPtr config) {
    if (!config->isLoaded()) {
        return;
    }
    std::string section("ui_gl_MouseAdapter");
    int button = -1;
    button = toButton(config->getString(section, "b_look"));
    if (button != -1) b_look_ = button;
    button = toButton(config->getString(section, "b_rotate"));
    if (button != -1) b_rotate_ = button;
    button = toButton(config->getString(section, "b_pan"));
    if (button != -1) b_pan_ = button;
    button = toButton(config->getString(section, "b_forward"));
    if (button != -1) b_forward_ = button;
    button = toButton(config->getString(section, "b_backpedal"));
    if (button != -1) b_backpedal_ = button;

    double sensitivity = 0.0;
    sensitivity = config->getDouble(section, "sensitivity_look");
    if (sensitivity != 0.0) sensitivity_look_ = sensitivity;
    sensitivity = config->getDouble(section, "sensitivity_rotate");
    if (sensitivity != 0.0) sensitivity_rotate_ = sensitivity;
    sensitivity = config->getDouble(section, "sensitivity_pan");
    if (sensitivity != 0.0) sensitivity_pan_ = sensitivity;
    sensitivity = config->getDouble(section, "sensitivity_move");
    if (sensitivity != 0.0) sensitivity_move_ = sensitivity;
}

void MouseAdapter::pressed(int button, int x, int y) {
    this->prev_pos_[0] = x;
    this->prev_pos_[1] = y;
    if (button == this->b_look_) this->mouse_look_ = true;
    if (button == this->b_rotate_) this->mouse_rotate_ = true;
    if (button == this->b_pan_) this->mouse_pan_ = true;
    if (button == this->b_forward_) this->camera_->moveFB(sensitivity_move_);
    if (button == this->b_backpedal_) this->camera_->moveFB(-sensitivity_move_);
}

void MouseAdapter::released(int button, int x, int y) {
    if (button == this->b_look_) this->mouse_look_ = false;
    if (button == this->b_rotate_) this->mouse_rotate_ = false;
    if (button == this->b_pan_) this->mouse_pan_ = false;
}

void MouseAdapter::dragged(int x, int y) {
    int dx = x - this->prev_pos_[0];
    int dy = y - this->prev_pos_[1];
    this->prev_pos_[0] = x;
    this->prev_pos_[1] = y;
    if (this->mouse_look_) {
        this->camera_->lookLR(dx * this->sensitivity_look_);
        this->camera_->lookUD(-dy * this->sensitivity_look_);
    }
    if (this->mouse_rotate_) {
        this->camera_->rotateWorldLR(dx * this->sensitivity_rotate_);
        this->camera_->rotateWorldUD(-dy * this->sensitivity_rotate_);
    }
    if (this->mouse_pan_) {
        this->camera_->strafeLR(-dx * this->sensitivity_pan_);
        this->camera_->strafeUD(dy * this->sensitivity_pan_);
    }
}

} }
