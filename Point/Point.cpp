//Title: point.cpp
//Description: Class of points -- location updates
//Platform(s): gcc (Ubuntu/Linaro 4.6.3-1ubuntu5) 4.6.3
//Date:     19/05/2017
//Revision: 19/05/2017

#include "Point.h"

Point::Point()
{
	x   = 0.0;
	y   = 0.0;
	t   = 0.0;
	v   = 0.0;
	oid = 0;
	qid 	  = 0;
	qid_cheat = 0;
}

Point::Point(double _x, double _y, double _t, unsigned int _oid, unsigned int _qid)
{
	x   = _x;
	y   = _y;
	t   = _t;
	v   = 0.0;    //Initially velocity is zero
	oid = _oid;
	qid 	  = _qid;
	qid_cheat = _qid;
}

Point::Point(const Point& p)
{
	x   = p.x;
	y   = p.y;
	t   = p.t;
	v   = p.v;
	oid = p.oid;
	qid 	  = p.qid;
	qid_cheat = p.qid_cheat;
}

Point::~Point()
{
	//cout << "Point deleted" << endl;
}

void Point::printPoint()
{
	cout << "<" << t << ", " << oid << ", " << x << ", " << y << ">" << endl;   
}

Point* Point::copyPoint(Point& p)
{
	Point *q = new Point(p);

	return q;
}