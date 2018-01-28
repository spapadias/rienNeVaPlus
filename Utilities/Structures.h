//Title: structures.h
//Description: Basic data structures and function for manipulating spatial grid partitioning.
//Platform(s): gcc (Ubuntu/Linaro 4.6.3-1ubuntu5) 4.6.3
//Date:     29/12/2007
//Revision: 11/05/2017

#ifndef STRUCTURESGRID_H
#define STRUCTURESGRID_H

#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iterator>
#include <iomanip>
#include <list>
#include <map>
#include <queue>
#include <set>
#include <string.h>
#include <sstream> 
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <utility>
#include <unistd.h>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <random>

using namespace std;

#define DELIMITER " "				
#define SEPARATOR ","
//Rome box
#define X_MIN 1757000
#define Y_MIN 4608800
#define X_MAX 1822000
#define Y_MAX 4673800

fstream fout; //Output file for processing results

// Query result ---------------------
typedef struct qRes
{
    double max;                       //Max of velocities
    double sum;						  //Sum of velocities
    unsigned int count;
    unsigned int qid;
    unordered_set<unsigned int> oids; //Object ids that contribute to this results

    qRes(unsigned int a)
    {
        max   = 0.0;
        sum   = 0.0;
        count = 0;
        qid   = a;
        oids.clear();
    }

    ~qRes() 
    {
        oids.clear();
    }
} queryRes;

// Boxes -------------------
typedef struct boxInfo 
{
	double x_min;
	double y_min;
	double x_max;
	double y_max;

	boxInfo(const double _x_min, const double _y_min, const double _x_max, const double _y_max) 
	{
        x_min = _x_min;
        y_min = _y_min;
        x_max = _x_max;
        y_max = _y_max;
	}

	boxInfo() {};
} box;

// Vertices ---------------
typedef struct vertex2D 
{
    double x;
    double y;

    vertex2D(double _x, double _y)
    {
        x = _x;
        y = _y;
    }
} vertex;

// Polygons --------------
typedef struct poly2D 
{											//ASSUMPTION: Νo holes inside polygons
    vector<vertex *> vertices;				//Vertices defining the external ring of this polygon
} polygon;

// Queries ---------------
typedef struct Query_tuple 
{
    double t;					
	unsigned int qid;		
    unsigned int cid;
    unsigned int range;
    unsigned int slide;
    //Query Buffers
	vertex*  startV;
	vertex*  endV;
	polygon* qArea;
} TQuery;

double max(double a, double b)
{
    if (a >= b) return a;
    else return b;
}

double min(double a, double b)
{
    if (a <= b) return a;
    else return b;
}

double min3(double f1, double f2, double f3)
{
    double minf = (f1 < f2) ? f1 : f2;
    minf = (minf < f3) ? minf : f3;
    return minf;
}

//Checks whether an object's location is contained within the given rectangle.
bool pointInRect(double x, double y, box rect)
{
	return ( ((x >= rect.x_min) && (x <= rect.x_max) && (y >= rect.y_min) && (y <= rect.y_max)) ? true : false);
}

//Checks if two segments overlap along the same axis, i.e. they have a common interval.
bool segmentOverlap(double a_min, double a_max, double b_min, double b_max)
{
	return ( ((a_max < b_min) || (b_max < a_min)) ? false : true);
}

//Intersection exists if only these rectangles have overlapping extents over both axes.
bool rectIntersect(box rect1, box rect2)
{
	return (segmentOverlap(rect1.x_min, rect1.x_max, rect2.x_min, rect2.x_max) && segmentOverlap(rect1.y_min, rect1.y_max, rect2.y_min, rect2.y_max));
}

//Checks whether the first rectangle is FULLY contained within the second one.
bool rectContain(box rect1, box rect2)
{
	return ((rect1.x_min >= rect2.x_min) && ((rect1.y_min >= rect2.y_min) && (rect1.x_max <= rect2.x_max) && (rect1.y_max <= rect2.y_max)) ? true : false);
}

#endif /* STRUCTURESGRID_H */

// // Point ----------------------------
// typedef struct point
// {
//     double x;
//     double y;
//     double t;
//     double v;
//     unsigned int oid;
//     unsigned int qid;

//     point(double _x, double _y, double _t, unsigned int _oid, unsigned int _qid)
//     {
//         x   = _x;
//         y   = _y;
//         t   = _t;
//         v   = 0.0;    //Initially velocity is zero
//         oid = _oid;
//         qid = _qid;
//     }

//     point(const Point& p)
//     {
//         x   = p.x;
//         y   = p.y;
//         t   = p.t;
//         v   = p.v;
//         oid = p.oid;
//         qid = p.qid;
//     }

//     ~point()
//     {
//         cout << "Point deleted" << endl;
//     }
// } Point;

// void quickSortVector(vector<Point>& candidate, int l, int r)
// {
//     int i, j;
//     double mid;
//     Point tmp;
//     i = l; j = r; mid = candidate[(i + j) / 2].value;
//     while (i <= j)
//     {
//         while (candidate[i].value < mid) i++;
//         while (candidate[j].value > mid) j--;
//         if (i <= j)
//         {
//             tmp = candidate[i]; candidate[i] = candidate[j]; candidate[j] = tmp;
//             i++; j--;
//         }
//     }
//     if (i<r) quickSortVector(candidate, i, r);
//     if (j>l) quickSortVector(candidate, l, j);
// }

// void quickSort(double *value, int *x, int *y, int l, int r)
// {
//     int i, j;
//     double mid;
//     double tmp;
//     int tmpInt;
//     i = l; j = r; mid = value[(i + j) / 2];
//     while (i <= j)
//     {
//         while (value[i] < mid) i++;
//         while (value[j] > mid) j--;
//         if (i <= j)
//         {
//             tmp = value[i]; value[i] = value[j]; value[j] = tmp;
//             tmpInt = x[i]; x[i] = x[j]; x[j] = tmpInt;
//             tmpInt = y[i]; y[i] = y[j]; y[j] = tmpInt;
//             i++; j--;
//         }
//     }
//     if (i<r) quickSort(value, x, y, i, r);
//     if (j>l) quickSort(value, x, y, l, j);
// }