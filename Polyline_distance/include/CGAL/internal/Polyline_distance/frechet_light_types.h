#pragma once

#include "geometry_basics.h"
#include "id.h"
#include "curve.h"

#include <vector>

//
// Box
//

struct Box {
	PointID min1;
	PointID max1;
	PointID min2;
	PointID max2;

	Box(PointID min1, PointID max1, PointID min2, PointID max2)
		: min1(min1), max1(max1), min2(min2), max2(max2) {}

	bool isCell() const {
		return max1 - min1 == 1 && max2 - min2 == 1;
	}
};
using Boxes = std::vector<Box>;

//
// Inputs
//

struct Inputs {
	CIntervals::iterator begin1;
	CIntervals::iterator end1;
	CIntervals::iterator begin2;
	CIntervals::iterator end2;

	bool haveDownInput() const { return begin1 != end1; }
	bool haveLeftInput() const { return begin2 != end2; }

	bool downContains(PointID point_id) const {
		for (auto it = begin1; it != end1; ++it) {
			if (it->begin <= point_id && it->end >= point_id) {
				return true;
			}
		}
		return false;
	}
	bool leftContains(PointID point_id) const {
		for (auto it = begin2; it != end2; ++it) {
			if (it->begin <= point_id && it->end >= point_id) {
				return true;
			}
		}
		return false;
	}
};

//
// Outputs
//

struct Outputs {
	CIntervalsID id1;
	CIntervalsID id2;
};

//
// QSimpleInterval
//

struct QSimpleInterval
{
	QSimpleInterval() : valid(false) {}
	QSimpleInterval(CPoint const& begin, CPoint const& end)
		: valid(true), free(begin, end) {}

	void setFreeInterval(CPoint const& begin, CPoint const& end) {
		free = {begin, end};
		//valid = true;
	}
	void setFreeInterval(PointID begin, PointID end) {
		free = {begin, 0., end, 0.};
		//valid = true;
	}
	void invalidateFreeInterval() {
		free.begin = CPoint{};
		free.end = CPoint{};
	}

	void setOuterInterval(CPoint const& outerbegin, CPoint const& outerend) {
		outer.begin = outerbegin;
		outer.end = outerend;
	}
	void setOuterInterval(PointID begin, PointID end) {
		setOuterInterval(CPoint(begin, 0.), CPoint(end, 0.));
	}

	CInterval const& getFreeInterval() {
		return free;
	}
	CInterval const& getOuterInterval() {
		return outer;
	}

	void setLastValidPoint(PointID const& point) {
		last_valid_point = point;
		//valid = false;
	}
	PointID const& getLastValidPoint() const {
		return last_valid_point;
	}
	bool hasPartialInformation() const {
		return last_valid_point.valid();
	}

	void validate() { valid = true; }
	void invalidate() { valid = false; }
	void clamp(CPoint const& min, CPoint const& max) {
		free.clamp(min, max);
		outer.clamp(min,max);
	} 

	bool is_empty() const { return free.is_empty(); }
	bool is_valid() const { return valid; }


private:
	bool valid;
	CInterval free;
	PointID last_valid_point;

	CInterval outer;
};

using QSimpleIntervals = std::vector<QSimpleInterval>;
using QSimpleID = ID<QSimpleInterval>;

//
// QSimpleOutputs
//

struct QSimpleOutputs
{
	QSimpleID id1;
	QSimpleID id2;
};

//
// BoxData
//

struct BoxData {
	Box box;
	Inputs inputs;
	Outputs outputs;
	QSimpleOutputs qsimple_outputs;
};
