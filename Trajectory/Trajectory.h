//Title: Trajectory.h
//Description:  Class that defines the Trajectory of a moving object
//Author: Kostas Patroumpas
//Platform(s): gcc (Ubuntu/Linaro 4.6.3-1ubuntu5) 4.6.3
//Date:     03/03/2010
//Revision: 20/05/2017

#ifndef TRAJECTORY_H
#define TRAJECTORY_H
//------------------------------------------------------------
#include "../Utilities/Structures.h"
#include "../Point/Point.h"
#include "../RegularGrid/RegularGrid.h"
#include "../Synopsis_DEPRECATED/Synopsis.h"
//------------------------------------------------------------
//class Synopsis;
//------------------------------------------------------------

class Trajectory
{

public:

    Trajectory(
                RegularGrid *, 
                unsigned int, 
                double, 
                unordered_map<unsigned int, Point*>       *,
                unordered_map<unsigned int, Point*>       *,  
                unordered_map<unsigned int, TQuery*>      *,
                unordered_map<unsigned int, unsigned int> *
              );

    ~Trajectory();

    int findQid(Point *, unsigned int);
    int findVelocity(Point *);
    void insertPoint(Point *);
    unsigned int evictPoint(double);
    double getDistance(Point *, Point *);
    unsigned int countPoints(); 
    void reportState();
    double measureTraveledDistance();

    class CompareCos
    {
        public:
            bool operator()(const std::pair<unsigned int, double> &left, const std::pair<unsigned int, double> &right) 
            {   
                if ( left.second > right.second )
                    return true;
                else
                    return false;
            }    
    };

    //-----------------
    RegularGrid  *grid;
    Synopsis     *synopsis;
    //---------------------
    list<Point *> trajList;                         //The last and preLast location updates could be missing from the trajList!
    unordered_map<unsigned int, Point*> *preLast;   //p_preprev location updates
    unordered_map<unsigned int, Point*> *last;      //p_prev location updates
    unordered_map<unsigned int, TQuery*> *Q;        //Valid queries in the system
    unordered_map<unsigned int, unsigned int> *U;
    unsigned int oid;
    double mrut;                                    //Most Recent Update Time
    double distance;                                //Total traveled distance

private:    
    
};

#endif /* TRAJECTORY_H */