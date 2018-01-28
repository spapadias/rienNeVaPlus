//Title: Trajectory.cpp 
//Description: Create a sliding windowing construct over trajectory features.
//Author: Kostas Patroumpas, Serafeim Papadias
//Platform(s): gcc (Ubuntu/Linaro 4.6.3-1ubuntu5) 4.6.3
//Date: 3/3/2010
//Revision: 4/3/2017

#include "Trajectory.h"

//Constructor for the trajectory referring to specific object
Trajectory::Trajectory(unsigned int pOid , double t0 , RegularGrid *pgrid)
{
    grid = pgrid;
    synopsis = new Synopsis(this);

    oid = pOid;                               
    lastUpdateTime = t0;                      //Time of last refresh; initially, it coincides with the time that the window is firstly being applied
    
    seq = new sequence(oid);                 

    travelDistance = 0.0;
}

//Destructor
Trajectory::~Trajectory()
{
}

//Insert a new point into this trajectory
double Trajectory::insertPoint(point *p, unordered_map< unsigned int, unsigned int > *U, int t, int s, unordered_map< unsigned int, unsigned int > *lastQID, unsigned int trajExists, unordered_set<unsigned int> *randomQ)
{
    double ret = synopsis->updateSynopsis(p, U, t, s, lastQID, randomQ);

    if ( trajExists == 0 )
        seq->segPoints.push_back(p);
    else //if ( trajExists == 1 )
    {
        if ( (*U)[ oid ] == 1 )                     //Maintain trajectories for important synopses only
            if ( p->qid != 0 )                      //Do not include no-in-query tuples into trajectories 
                seq->segPoints.push_back(p);
    }

    //Update the last timestamp that this object received a location
    lastUpdateTime = p->t;

    return ret;
}

//Check oldest points for expiration from the window state according to the given time bound
unsigned int Trajectory::check4ExpiringPoints(double t_rear)
{    
    unsigned int num = 0;
    
    //Multiple points may be expiring at each round
    point *q;
          
    //Remove expiring positions from the sequence
    while ( ( !seq->segPoints.empty() ) && ( (q= seq->segPoints.front())->t < t_rear ) )                     
    {   
        seq->segPoints.pop_front();                      
        delete q;                                    //Deallocates the memory block pointed by q       
        num++;
    }
       
    //Update time of last refresh
    if ( !seq->segPoints.empty() )
        lastUpdateTime = seq->segPoints.back()->t;   //lastUpdateTime is the timestamp of last update 

    //Update the trajectory synopsis as well
    //synopsis->cropSynopsis( t_rear );
    
    return num;                                      //Number of expired locations
}

//Assigns qID to incoming tuple before included to trajectory
unsigned int Trajectory::findQid(point *p, unsigned int mode)
{
    unsigned int cellID;
    unsigned int flagPNP = 0;                                   //Flag variable for PointInPolygon Algorithm
    int result;                                                 //The result of PNP : 1 for inside , 0 for not inside
    map<unsigned int , vector<TQuery*> >::iterator qIt;         //Iterator to the input map structure
    double previous_x,previous_y;
    list<unsigned int> listQueries;                             //List to keep all the queries that a new tuple falls into
    vector<TQuery*>::iterator vecIt;                            //Iterator to the parts of each unique qID/arterial road
    double scalarProd;                                          //Here we store the scalarProduct that we calculate between orientation vectors
    double lengthV1, lengthV2;                                  //The length of the vectors we take their scalar/dot product
    double cosine;                                              //The cosΘ of the angle θ between the two vertices 

    //We steal the qID and read it from the input, so save the right qID for this tuple/point here.
    unsigned int correctID = p->qid;
    //Define a priority queue
    priority_queue< pair<unsigned int, double>, vector< pair< unsigned int, double> >, CompareCos> priorityQueueCos;

    if (mode == 1)
    {
        //Here we "steal" and assign to the first's tuple of each and every new created trajectory qID, the qID that we read from the Input file
        return 1;                                               //Return True
    }
    else if (mode == 2)
    {
        //Trajectory found!
        //Find the x and y coordinates of the last point of the trajectory
        previous_x = seq->segPoints.back()->x;
        previous_y = seq->segPoints.back()->y; 

        cellID = grid->HashLocation( p->x , p->y );

        if ( grid->cell[cellID].QryInfo.size() != 0 ) 
        {
            flagPNP = 0;                                        //Initialize the flag
            for ( qIt = grid->cell[cellID].QryInfo.begin() ; qIt != grid->cell[cellID].QryInfo.end() ; qIt++ )
            {
                for ( vecIt = qIt->second.begin(); vecIt != qIt->second.end(); vecIt++ )
                {
                    //We assume that the points of RandomTQuery are given clock or anti-clockwise                       
                    result = grid->pnpoly( (*vecIt)->qArea->vertices.size() , (*vecIt)->qArea , p->x , p->y ); 
                    //cout << "PNP result: " << result << endl;
                    if ( result == 1 )
                    {
                        //The points is inside this box. Although we must determine whether to assign it here or not
                        lengthV1 = sqrt( pow( (*vecIt)->endV->x - (*vecIt)->startV->x, 2 ) + pow( (*vecIt)->endV->y - (*vecIt)->startV->y, 2 ) ); //cout << lengthV1 << " ";
                        lengthV2 = sqrt( pow( p->x - previous_x, 2 ) + pow( p->y - previous_y, 2 ) ); //cout << lengthV2 << " ";
                        scalarProd = ( (*vecIt)->endV->x - (*vecIt)->startV->x )*( p->x - previous_x ) + ( (*vecIt)->endV->y - (*vecIt)->startV->y )*( p->y - previous_y ); //cout << scalarProd << " ";  
                        cosine = scalarProd/( lengthV1 * lengthV2 ); //cout << cosine << endl;

                        //Construct a priority queue with pairs (qID,cosCalculated) and use a minimum sort function to top() the min cosine
                        priorityQueueCos.push( make_pair( qIt->first , (1 - cosine) ) ); //Push the pair into the queue
                        flagPNP = 1; //Make flag True
                    }
                }
            }
            //If flagPNP = 1 then we have found the qID for the new incoming tuple
            if ( flagPNP == 0 ) {
                p->qid = 0 ;                                    //Tuple does not belong to any query - flagPNP = 0
            }
            if ( flagPNP == 1 ) {
                p->qid = priorityQueueCos.top().first;
                //cout << priorityQueueCos.top().first << " " << priorityQueueCos.top().second << endl;
            }
        }
        else {
            p->qid = 0 ;                                        //Tuple does not belong to any query
        }

        //Check if the qid found is correct
        if ( p->qid == correctID )
        {
            //cout << "TRUE" << endl;
            p->qid = correctID;                                 //Assign the correct qID to this tuple to be totally correct
            return 1;                                           //Return True - Correctly asssigned to that road
        }
        else
        {
            //cout << "FALSE" << endl;
            p->qid = correctID;                                 //Assign the correct qID to this tuple to be totally correct
            return 0;                                           //Return False - Incorrectly assigned to that road
        }
    }
    else 
        return 0;                                               //Return the qid that was assigned to the incoming tuple
}

//Assigns velocity of tuple before included to trajectory list
double Trajectory::findVelocity(point *p, unsigned int mode)
{
    point *previous_point;

    if (mode == 1)
    {
        p->veloPoint = 0.0;                          //First point => //CAUTION: Wrong to assign "zero velocity" to first tuples
    }
    else if (mode == 2)
    {
        /* Find the velocity between this and the previous point - do this for all points */
        previous_point = seq->segPoints.back();
        //Adjust the velocity of the new incoming point before you insert it
        //1 m/sec = 3.6 km/hour
        p->veloPoint = 3.6 * sqrt( pow((p->x - previous_point->x),2) + pow((p->y - previous_point->y),2) ) / (p->t - previous_point->t);
    }

    //Correct the constant velocity of the first tuple in queries for each trajectory segment in each queries
    if ( seq->segPoints.size() == 2 )
        if ( seq->segPoints.front()->veloPoint < 0.00001 )
        {
            //cout << "fixed velocity!" << endl;
            seq->segPoints.front()->veloPoint = p->veloPoint;
        }

    return p->veloPoint;
}

//Finds correctly the instant velocity of a newly incoming point, even if the transition is inactive ==> active
double Trajectory::findVelocitySpecial(point *p, point *lastPointKept)
{
    //1 m/sec = 3.6 km/hour
    p->veloPoint = 3.6 * sqrt( pow((p->x - lastPointKept->x), 2) + pow((p->y - lastPointKept->y), 2) ) / ( p->t - lastPointKept->t );

    //Correct the constant velocity of the first tuple in queries for each trajectory segment in each queries
    if ( seq->segPoints.size() == 2 )
        if ( seq->segPoints.front()->veloPoint < 0.00001 )
        {
            //cout << "fixed velocity!" << endl;
            seq->segPoints.front()->veloPoint = p->veloPoint;
        }

    return p->veloPoint;                           
}

//Function which returns the score of the trajectory
double Trajectory::scoringOfTrajectory( unordered_map< unsigned int, unsigned int > *U, unsigned int boolWeightedScore )
{
    //Check the trajectories that have not been selected from the first criterion
    if ( (*U)[ oid ] == 0 )
        return synopsis->scoring(boolWeightedScore);
    else
        return -2.0; //-2.0 -- Trajectory already selected by the first criterion
}

//Count point locations participating in this trajectory (referring to a specific object)
unsigned int Trajectory::countPoints()
{
    return seq->segPoints.size();
}

//Calculate distance along the sequence of points
double Trajectory::measureTraveledDistance()
{
    double dist = 0.0f;
    list<point *>::iterator pit;                         //Previous point
    list<point *>::iterator it;                          //Current point
    it = seq->segPoints.begin();
    pit = seq->segPoints.begin();
    it++;

    while (it != seq->segPoints.end())	
    {
        dist += getDistance(*pit, *it);
        pit++;
        it++;  
    }

    travelDistance = dist;                               //Distance in the same units as the original coordinates
    return dist;
}

//Print to the given file the current state of this trajectory, i.e., all points retained in this instantiation
void Trajectory::printState(ofstream &sout)
{
    list<point *>::iterator it;
    for (it=seq->segPoints.begin(); it != seq->segPoints.end(); it++ )
        sout << "\r\n" << (*it)->t << DELIMITER << setprecision(10) << fixed << (*it)->x << DELIMITER << setprecision(10) << fixed << (*it)->y;   
}

//Report (= print to the standard output) the current state of this trajectory, i.e., the sequence of points retained in this instantiation
void Trajectory::reportState()
{
    list<point *>::iterator it;
    cout << "{";
    for (it=seq->segPoints.begin(); it != seq->segPoints.end(); it++ )
        //cout << "(" << setprecision(6) << (*it)->t << SEPARATOR << (*it)->x << SEPARATOR << (*it)->y << SEPARATOR << (*it)->veloPoint << SEPARATOR << (*it)->qid << ") -> ";   
        cout << "(" << setprecision(6) << (*it)->qid << SEPARATOR << setprecision(2) << (int)(*it)->t << setprecision(6) << ") -> ";   
    cout << "}";
    //And after that print the synopsis' info
    cout << endl;
    synopsis->printSynopsis();
}