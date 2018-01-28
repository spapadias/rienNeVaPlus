//Title: Trajectory.cpp 
//Description: Create a sliding windowing construct over trajectory features.
//Author: Kostas Patroumpas, Serafeim Papadias
//Platform(s): gcc (Ubuntu/Linaro 4.6.3-1ubuntu5) 4.6.3
//Date:     03/03/2010
//Revision: 20/05/2017

#include "Trajectory.h"

Trajectory::Trajectory(
                        RegularGrid *_grid, 
                        unsigned int _oid, 
                        double _t0, 
                        unordered_map<unsigned int, Point*>       *_preLast, 
                        unordered_map<unsigned int, Point*>       *_last, 
                        unordered_map<unsigned int, TQuery*>      *_Q,
                        unordered_map<unsigned int, unsigned int> *_U
                      )
{
    grid     = _grid;   
    oid      = _oid ;                               
    mrut     = _t0  ;    //Time of latest refresh; initially, it coincides with the time that the window is firstly being applied
    synopsis = new Synopsis(this->grid, this);
    distance = 0.0  ;
    preLast  = _preLast; //Map with the p_preprev for all objects
    last     = _last;    //Map with the p_prev for all objects    -- Total space O(2N) 
    Q        = _Q;
    U        = _U;
}

Trajectory::~Trajectory()
{
    //Delete the pointer to the sequence of location updates
    for ( auto i = trajList.begin() ; i != trajList.end() ; i++ )
        delete (*i);
    trajList.clear();
    //Delete the synopsis of the trajectory
    synopsis->~Synopsis(); //Destructor of Synopsis class
//    cout << "Trajectory of " << oid << " deleted" << endl;
}

double Trajectory::getDistance(Point *p1, Point *p2)
{
    return sqrt(pow((p1->x - p2->x),2) + pow((p1->y - p2->y),2));
} 

int Trajectory::findQid(Point *p, unsigned int mode)
{
    unsigned int cid;//, flagPNP = 0;                                   
    double last_x, last_y, scalar, len1, len2, cosine; 
    priority_queue< pair<unsigned int, double>, vector< pair< unsigned int, double> >, CompareCos> queueCos;

    //--------------------------------
    unsigned int correct_qid = p->qid;
    //--------------------------------

    if      ( mode == 1 )
        return 1; //Classify the first tuple of each moving object correctly
    else if ( mode == 2 )
    {
        //Fetch the x & y coordinates
        last_x = (*last)[oid]->x;
        last_y = (*last)[oid]->y;
        cid = grid->HashLocation(p->x, p->y); //Hash the location to grid
        
        if ( grid->cell[cid].QryList.size() != 0 )
        {
            for ( auto i = grid->cell[cid].QryList.begin() ; i != grid->cell[cid].QryList.end() ; i++ )
            {
                if ( grid->pnpoly( i->second->qArea->vertices.size(), i->second->qArea, p->x, p->y ) == 1 )
                {
                    len1 = sqrt( 
                                pow(i->second->endV->x - i->second->startV->x, 2) 
                                + 
                                pow(i->second->endV->y - i->second->startV->y, 2) 
                               );
                    len2 = sqrt( 
                                pow( p->x - last_x, 2) 
                                +
                                pow( p->y - last_y, 2) 
                               );
                    scalar = (i->second->endV->x - i->second->startV->x) * (p->x - last_x) 
                             + 
                             (i->second->endV->y - i->second->startV->y) * (p->y - last_y);
                    cosine = ( (double)scalar) / ( len1 * len2 );
                    queueCos.push( make_pair(i->first, (1 - cosine)) );
                }
            }

            if ( queueCos.empty() )
                p->qid = 0;
            else
            {
                p->qid = queueCos.top().first;
                //Check qid for validity in Q
                if ( Q->find(p->qid) == Q->end() )
                    p->qid = 0;
            }
        }
        else
            p->qid = 0;

        // //If qid is not in the valid query set then set its qit to 0 -----
        // if ( Q->find(p->qid) == Q->end() )
        //     p->qid = 0;
        // //----------------------------------------------------------------

        //--------------
        //For Statistics
        if ( p->qid != 0 )
        {
            //Check if the qid found is correct
            if ( p->qid == correct_qid )
            { 
                p->qid = correct_qid; //CHEATING!                              
                return 1;                                       //Return True  - Correctly asssigned to that road
            }
            else
            { 
                p->qid = correct_qid; //CHEATING!
                return 0;                                       //Return False - Incorrectly assigned to that road
            }    
        }
        else //If qid == 0 return -1000
            return -1000;                                       //Non-query code

        cout << "Never goes through here" << endl;
        return 1000;
    }

    //---------------------------------------------------------------------------
    //---------------------------------------------------------------------------
    cout << "Never goes through here" << endl;
    return 1001;                                                //Wrong mode code
}

int Trajectory::findVelocity(Point *p)
{
    //Calculate the instantenuous velocity
    p->v = 3.6 * getDistance((*last)[this->oid], p) / (p->t - (*last)[this->oid]->t); //1 m/sec = 3.6 km/hour

    if ( p->v > 130.0 )
        return -1;      //Not well map matched point due to the low sampling rate -- Map matching algorithm issue

    //Normal velocity -- point to be inserted
    return 1;
}

void Trajectory::insertPoint(Point *p)
{
    //Update the synopsis of this trajectory before inserting the new incoming point "p" to the list of points
    this->synopsis->updateC(p);
    this->synopsis->updateQ(p);

    this->trajList.push_back(p);

    //---
    unsigned int cid = this->grid->HashLocation(p->x, p->y);  //Hash the location to grid
    this->grid->cell[cid].HisObjList.insert(p->oid);          //Insert this oid to the corresponding cell
    //---
}

unsigned int Trajectory::evictPoint(double t_rear)
{    
    unsigned int num = 0;
    
    //Multiple points may be expiring at each round
    Point *q;
  
    //---
    unsigned int last_cid = 0;
    //---
    //Keep only location updates in (t_rear,t_now], where t_rear = t_now - win_range             
    while ((!(trajList.empty())) && ((q= trajList.front())->t < t_rear))
    {     
        //---
        //---
        unsigned int cur_cid = this->grid->HashLocation(trajList.front()->x, trajList.front()->y);
        //---
        if ( cur_cid != last_cid )
        {
            //---
            if ( this->grid->cell[last_cid].HisObjList.find(this->oid) != this->grid->cell[last_cid].HisObjList.end() )
                this->grid->cell[last_cid].HisObjList.erase(this->oid);
            //---
            last_cid = cur_cid;
        }
        //---
        //---

        trajList.pop_front();

        delete q;          
        num++;
    } //----------------------------------------------------------------------------
  
    //Update mrut-----------------
    if (!trajList.empty())
        mrut = trajList.back()->t;
    //----------------------------

    //Refine the synopsis of this trajectory
    this->synopsis->cropQ(t_rear);
    this->synopsis->cropC(t_rear);
    //--------------------------------------

    return num;    //Number of expired locations
}

unsigned int Trajectory::countPoints()
{
    return trajList.size();
}

double Trajectory::measureTraveledDistance()
{
    double dist = 0.0f;
    list<Point *>::iterator pit;                         //Previous point
    list<Point *>::iterator it;                          //Current point
    it  = trajList.begin();
    pit = trajList.begin();
    it++;

    while (it != trajList.end())	
    {
        dist += getDistance(*pit, *it);
        pit++;
        it++;  
    }

    this->distance = dist;                               //Distance in the same units as the original coordinates
    return dist;
}

void Trajectory::reportState()
{
    cout << "{";
    for (auto it = trajList.begin(); it != trajList.end(); it++ )
        cout << "(" << setprecision(6) << (*it)->qid << SEPARATOR << setprecision(2) << (int)(*it)->t << setprecision(6) << ") -> ";   
    cout << "}" << endl;
}