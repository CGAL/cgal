//A general purpose file reader.
//Note that if we define more file formats,
//it may be beneficial to split this into multiple plugin classes/files.
//Copyright (C) 2013  INRIA - Sophia Antipolis
//
//This program is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//This program is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// Author(s):      Thijs van Lankveld


#ifndef READ_FILE
#define READ_FILE

#include <fstream>
#include <iostream>

/*#include <osgDB/ReaderWriter>
#include <osgDB/ReadFile>*/

#include "quad.h"

const unsigned int LINE_SIZE = 1024;

void readLine(std::ifstream& fin, char* line) {
	fin.getline(line, LINE_SIZE);
	while (fin && line[0] == '#') {
		// Comment line.
		std::cout << "Comment: " << line << std::endl;
		fin.getline(line, LINE_SIZE);
	}
}


template < class Kernel, class PointOutput, class NormalOutput >
bool loadXYZN(char* file, PointOutput po, NormalOutput no) {
	typedef typename Kernel::Point_3	Point;
	typedef typename Kernel::Vector_3	Normal;

	char line[LINE_SIZE];

	std::ifstream fin(file);
	if (!fin) return false;

	unsigned int num = 0;

	// Collect the points.
	while (fin) {
		readLine(fin, line);

		if (strlen(line) > 0) {
			++num;

			// Collect the point and normal.
			double x, y, z;
			double nx, ny, nz;
			int a = sscanf(line,"%lf %lf %lf %lf %lf %lf", &x, &y, &z, &nx, &ny, &nz);

			if (a > 2)
				*po++ = Point(x, y, z);
			if (a > 5)
				*no++ = Normal(nx, ny, nz);
		}
	}

	fin.close();

	return true;
}


/*class VectorCollecter: public osg::ConstValueVisitor {
	inline bool comp(double d1, double d2) const {return abs(d1 - d2) < 0.0001;}

public:
	VectorCollecter(): osg::ConstValueVisitor() {}

	inline void apply(const osg::Vec3d& vec) {_vec3 = osg::Vec3d(vec[0],vec[1],vec[2]);}
	inline void apply(const osg::Vec3& vec) {_vec3 = osg::Vec3d(vec[0],vec[1],vec[2]);}
	inline void apply(const osg::Vec4d& vec) {_vec4 = osg::Vec4d(vec[0],vec[1],vec[2],vec[3]);}
	inline void apply(const osg::Vec4& vec) {_vec4 = osg::Vec4d(vec[0],vec[1],vec[2],vec[3]);}
	inline void apply(const osg::Vec4ub& vec) {_vec4 = osg::Vec4d(vec[0]/256.f, vec[1]/256.f, vec[2]/256.f, vec[3]/256.f);}

	inline bool compareColor(const osg::Vec4d& c) {return comp(_vec4[0], c[0]) && comp(_vec4[1], c[1]) && comp(_vec4[2], c[2]);}

	osg::Vec3d _vec3;
	osg::Vec4d _vec4;

	template < class ToType >
	ToType convert3() const {return ToType(_vec3[0], _vec3[1], _vec3[2]);}
	template < class ToType >
	ToType convert4() const {return ToType(_vec4[0], _vec4[1], _vec4[2], _vec4[3]);}
};

template < class Kernel, typename Color, class PointOutput, class NormalOutput, class ColorOutput >
class PointCollecter: public osg::NodeVisitor {
	typedef typename Kernel::Point_3			Point;
	typedef typename Kernel::Vector_3			Normal;

	PointOutput _po;
	NormalOutput _no;
	ColorOutput _co;
public:
	PointCollecter(PointOutput po, NormalOutput no, ColorOutput co): osg::NodeVisitor(osg::NodeVisitor::TRAVERSE_ALL_CHILDREN), _po(po), _no(no), _co(co) {}

	void apply(osg::Geode& geode) {
		VectorCollecter vc;
		for (unsigned int i = 0; i < geode.getNumDrawables(); ++i) {
			osg::Geometry* geom = geode.getDrawable(i)->asGeometry();

			// Only collect points.
			if (geom && geom->getNumPrimitiveSets() > 0 && geom->getPrimitiveSet(0)->getMode() == osg::PrimitiveSet::POINTS) {
				// Collect locations.
				osg::Array* vertices = geom->getVertexArray();
				for (unsigned int j = 0; j < vertices->getNumElements(); ++j) {
					vertices->accept(j, vc);
					*_po++ = vc.convert3<Point>();
				}

				// Collect normals.
				osg::Array* normals = geom->getNormalArray();
				if (geom->getNormalBinding() == osg::Geometry::BIND_OVERALL) {
					normals->accept(0, vc);
					Normal n = vc.convert3<Normal>();
					for (unsigned int j = 0; j < vertices->getNumElements(); ++j)
						*_no++ = n;
				}
				else if (normals) {
					for (unsigned int j = 0; j < normals->getNumElements() && j < vertices->getNumElements(); ++j) {
						normals->accept(j, vc);
						*_no++ = vc.convert3<Normal>();
					}
					for (unsigned int j = 0; j < (vertices->getNumElements() - normals->getNumElements()); ++j) {
						*_no++ = Normal(0, 0, 0);
					}
				}

				// Collect colors.
				osg::Array* colors = geom->getColorArray();
				if (geom->getColorBinding() == osg::Geometry::BIND_OVERALL) {
					colors->accept(0, vc);
					Color c = vc.convert4<Color>();
					for (unsigned int j = 0; j < vertices->getNumElements(); ++j)
						*_co++ = c;
				}
				else if (colors) {
					for (unsigned int j = 0; j < colors->getNumElements() && j < vertices->getNumElements(); ++j) {
						colors->accept(j, vc);
						*_co++ = vc.convert4<Color>();
					}
					for (unsigned int j = 0; j < (vertices->getNumElements() - colors->getNumElements()); ++j) {
						*_co++ = Color(0, 0, 0, 0);
					}
				}
			}
		}
	}
};

template < class Kernel, typename Color, class PointOutput, class NormalOutput, class ColorOutput >
bool loadOSG(char* file, PointOutput po, NormalOutput no, ColorOutput co) {
	osg::ref_ptr<osg::Node> scene = osgDB::readNodeFile(file);
	PointCollecter<Kernel, Color, PointOutput, NormalOutput, ColorOutput> pc(po, no, co);
	scene->accept(pc);
	return scene.valid();
}*/

#endif // READ_FILE