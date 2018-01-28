//Title: Trajectory.h
//Description: Create a sliding windowing construct over trajectory features.
//Author: Kostas Patroumpas
//Platform(s): gcc (Ubuntu/Linaro 4.6.3-1ubuntu5) 4.6.3
//Date: 3/3/2010
//Revision: 4/3/2017

#ifndef TRAJECTORY_H
#define TRAJECTORY_H //Inclusion guards

#include "structuresTraj.h"
#include "RegularGrid.h"
#include "Synopsis.h"

using namespace std;

//class Synopsis; /* Forward declaration of Synopsis */

class Trajectory
{

public:
    Trajectory(unsigned int, double, RegularGrid *);      
    ~Trajectory();
    double insertPoint(point *, unordered_map< unsigned int, unsigned int > *, int, int, unordered_map< unsigned int, unsigned int > *, unsigned int, unordered_set<unsigned int> *);
    unsigned int check4ExpiringPoints(double);
    unsigned int findQid(point *, unsigned int);    
    double findVelocity(point *, unsigned int);     
    double findVelocitySpecial(point *, point *);   
    double scoringOfTrajectory(unordered_map< unsigned int, unsigned int > *, unsigned int);
    unsigned int countPoints();
    double measureTraveledDistance();
    void printState(ofstream &);
    void reportState();

    //Class to compare the priority queue elements as we insert them into the priority queue
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

    sequence *seq;         //The sequence of points of the trajectory  
    Synopsis *synopsis;    //Pointer to the synopsis of this trajectory
    RegularGrid *grid;     //Pointer to the grid that we use
    unsigned int oid;
    double lastUpdateTime;
    double travelDistance;

private:    
    
};

#endif /* TRAJECTORY_H */