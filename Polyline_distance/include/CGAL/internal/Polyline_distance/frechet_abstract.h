#pragma once

#include "geometry_basics.h"
#include "certificate.h"
#include "curve.h"

#include <array>

class FrechetAbstract
{
public:
	virtual ~FrechetAbstract() {}
	virtual bool lessThan(distance_t distance, Curve const& curve1, Curve const& curve2) = 0;
	virtual Certificate&  computeCertificate() = 0;

	// yes, this is ugly...
	virtual void setRules(std::array<bool,5> const& enable) {}
	virtual void setPruningLevel(int pruning_level) {};
};
