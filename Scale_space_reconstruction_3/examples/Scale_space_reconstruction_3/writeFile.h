//A general purpose file writer.
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


#ifndef WRITE_FILE
#define WRITE_FILE

#include <map>

/*#include <osg/CullFace>
#include <osg/Material>
#include <osg/Switch>

#include <osgDB/ReaderWriter>
#include <osgDB/WriteFile>

template < class Point >
osg::Vec3 toOSG(const Point& p) {
	return osg::Vec3(CGAL::to_double(p.x()), CGAL::to_double(p.y()), CGAL::to_double(p.z()));
}

template < class PointIterator, class NormalIterator, class ColorIterator >
bool savePointsIVE(osg::Switch& shape,
				   PointIterator points_start, PointIterator points_end,
				   NormalIterator normals_start, NormalIterator normals_end,
				   ColorIterator colors_start, ColorIterator colors_end) {
	osg::ref_ptr<osg::Geode> geode = new osg::Geode;
	osg::ref_ptr<osg::Geometry> geom = new osg::Geometry;

	osg::ref_ptr<osg::Vec3Array> vert = new osg::Vec3Array;
	for (PointIterator pit = points_start; pit != points_end; ++pit)
		vert->push_back(toOSG(*pit));
	geom->setVertexArray(vert.get());

	osg::ref_ptr<osg::Vec3Array> norm = new osg::Vec3Array;
	norm->reserve(vert->getNumElements());
	for (NormalIterator nit = normals_start; nit != normals_end; ++nit)
		norm->push_back(toOSG(*nit));
	if (norm->getNumElements() != vert->getNumElements())
		return false;
	geom->setNormalBinding(osg::Geometry::BIND_PER_VERTEX);
	geom->setNormalArray(norm.get());

	osg::ref_ptr<osg::Vec4Array> colors = new osg::Vec4Array;
	if (colors_start == colors_end) {
		colors->reserve(1);
		colors->push_back(osg::Vec4(0.5,0.5,0.5,1.0));
		geom->setColorBinding(osg::Geometry::BIND_OVERALL);
	}
	else {
		colors->reserve(vert->getNumElements());
		for (ColorIterator cit = colors_start; cit != colors_end; ++cit)
			colors->push_back(osg::Vec4(cit->r(), cit->g(), cit->b(), cit->a()));
		geom->setColorBinding(osg::Geometry::BIND_PER_VERTEX);
	}
	geom->setColorArray(colors.get());

	geom->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::POINTS, 0, vert->size()));

	geode->addDrawable(geom.get());
	shape.addChild(geode.get());

	return true;
}

template < class FacetIterator, class NormalIterator, class ColorIterator >
bool saveTrianglesIVE(osg::Switch& shape,
					  FacetIterator facets_start, FacetIterator facets_end,
					  NormalIterator normals_start, NormalIterator normals_end,
					  ColorIterator colors_start, ColorIterator colors_end) {
	size_t size = std::distance(facets_start, facets_end);
	if (size != std::distance(normals_start, normals_end))
		return false;

	osg::ref_ptr<osg::Geode> geode = new osg::Geode;
	osg::ref_ptr<osg::Geometry> geom = new osg::Geometry;

	osg::ref_ptr<osg::Vec3Array> vert = new osg::Vec3Array;
	vert->reserve(3 * size);
	osg::ref_ptr<osg::Vec3Array> norm = new osg::Vec3Array;
	norm->reserve(size);
	NormalIterator nit = normals_start;
	for (FacetIterator tit = facets_start; tit != facets_end; ++tit, ++nit) {
		vert->push_back(toOSG(tit->first->vertex((tit->second+1)&3)->point()));
		vert->push_back(toOSG(tit->first->vertex((tit->second+2)&3)->point()));
		vert->push_back(toOSG(tit->first->vertex((tit->second+3)&3)->point()));
		norm->push_back(toOSG(*nit));
	}
	
	osg::ref_ptr<osg::Vec4Array> colors = new osg::Vec4Array;
	if (colors_start == colors_end) {
		colors->reserve(1);
		colors->push_back(osg::Vec4(0.5,0.5,0.5,0.5));
		geom->setColorBinding(osg::Geometry::BIND_OVERALL);
	}
	else {
		colors->reserve(vert->getNumElements());
		for (ColorIterator cit = colors_start; cit != colors_end; ++cit)
			colors->push_back(osg::Vec4(cit->r(), cit->g(), cit->b(), cit->a()));
		geom->setColorBinding(osg::Geometry::BIND_PER_VERTEX);
	}

	geom->setVertexArray(vert.get());
	geom->setNormalArray(norm.get());
	geom->setNormalBinding(osg::Geometry::BIND_PER_PRIMITIVE);
	geom->setColorArray(colors.get());
	geom->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::TRIANGLES, 0, vert->size()));
	geode->addDrawable(geom.get());

	osg::StateSet* state = geode->getOrCreateStateSet();
	state->setAttributeAndModes(new osg::CullFace(osg::CullFace::BACK));
//	osg::Material* material = new osg::Material;
//	material->setDiffuse(osg::Material::FRONT_AND_BACK, osg::Vec4(0.2, 0.2, 0.2, 0.25));
//	state->setAttributeAndModes(material);
	shape.addChild(geode.get());

	return true;
}

bool saveIVE(osg::Node* node, const char* file) {
	osg::ref_ptr<osgDB::ReaderWriter::Options> options = new osgDB::ReaderWriter::Options("precision 10");
	return osgDB::writeNodeFile(*node, file, options.get());
}*/

template < class TriangleIterator, class NormalIterator, class ColorIterator >
bool saveOFF(const char* file,
			 TriangleIterator triangles_start, TriangleIterator triangles_end,
			 NormalIterator normals_start, NormalIterator normals_end,
			 ColorIterator colors_start, ColorIterator colors_end) {
	typedef CGAL::Point_3<typename NormalIterator::value_type::R>	Point;

	// Construct a mapping from the points to their indices.
	std::map<Point, unsigned int> map;
	unsigned int ind = 0, tri = 0;
	for (TriangleIterator tit = triangles_start; tit != triangles_end; ++tit) {
		for (unsigned int i = 0; i < 4; ++i)
			if (map.find(tit->vertex(i)) == map.end())
				map[tit->vertex(i)] = ind++;
		++tri;
	}

	if (tri != std::distance(normals_start, normals_end))
		return false;

	std::ofstream fout(file);
	fout << "OFF" << std::endl;
	fout << map.size() << " " << tri << " 0" << std::endl;

	// Write the points.
	ColorIterator cit = colors_start;
	for (std::map<Point, unsigned int>::const_iterator pit = map.begin(); pit != map.end(); ++pit) {
		fout << pit->first;
		if (cit != colors_end)
			fout << " " << *cit++;
		fout << std::endl;
	}

	// Write the facets.
	for (TriangleIterator tit = triangles_start; tit != triangles_end; ++tit) {
		fout << "3 " << map[tit->vertex(0)]
			 << " " << map[tit->vertex(1)]
			 << " " << map[tit->vertex(2)] << std::endl;
	}

	return true;
}

#endif // WRITE_FILE