/*
 *   RapCAD - Rapid prototyping CAD IDE (www.rapcad.org)
 *   Copyright (C) 2010-2018 Giles Bathgate
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifdef USE_CGAL
#include "cgalpolygon.h"
#include "cgalprimitive.h"
#include <CGAL/normal_vector_newell_3.h>

CGALPolygon::CGALPolygon(CGALPrimitive* p) :
	Polygon(p),
	hasPlane(false),
	hole(false)
{
}

QList<CGAL::Point3> CGALPolygon::getPoints() const
{
	QList<CGAL::Point3> points;
	auto* pr=static_cast<CGALPrimitive*>(parent);
	QList<CGAL::Point3> parentPoints=pr->getPoints();
	for(auto i: indexes)
		points.append(parentPoints.at(i));
	return points;
}

QList<CGAL::Point2> CGALPolygon::getXYPoints() const
{
	QList<CGAL::Point2> points;
	auto* pr=static_cast<CGALPrimitive*>(parent);
	QList<CGAL::Point3> parentPoints=pr->getPoints();
	for(auto i: indexes) {
		CGAL::Point3 p3=parentPoints.at(i);
		points.append(CGAL::Point2(p3.x(),p3.y()));
	}
	return points;
}

CGAL::Vector3 CGALPolygon::getNormal() const
{
	return plane.orthogonal_vector();
}

void CGALPolygon::calculatePlane()
{
	if(hasPlane) return;
	CGAL::Vector3 v;
	QList<CGAL::Point3> points=getPoints();
	CGAL::normal_vector_newell_3(points.begin(),points.end(),v);
	plane=CGAL::Plane3(points.first(), v);
	hasPlane=true;
}

CGAL::Plane3 CGALPolygon::getPlane() const
{
	return plane;
}

void CGALPolygon::setPlane(const CGAL::Plane3& p)
{
	plane=p;
	hasPlane=true;
}

bool CGALPolygon::getHole() const
{
	return hole;
}

void CGALPolygon::setHole(bool value)
{
	hole = value;
}

#endif
