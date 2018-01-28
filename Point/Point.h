//Title: point.h
//Description: Class of points -- location updates
//Platform(s): gcc (Ubuntu/Linaro 4.6.3-1ubuntu5) 4.6.3
//Date:     19/05/2017
//Revision: 11/06/2017

#ifndef POINT_H
#define POINT_H

#include "../Utilities/Structures.h"

class Point
{
	public:
		Point();
		Point(double x, double y, double t, unsigned int oid, unsigned int qid);
		Point(const Point& p);
		~Point();

		void printPoint();
		Point* copyPoint(Point& p);

		double x;
		double y;
		double t;
		double v;
		unsigned int oid;
		unsigned int qid;
		//---------------------
		unsigned int qid_cheat; //This variable keeps the original qid of the incoming location update
		//---------------------
};

#endif /*POINT_H*/