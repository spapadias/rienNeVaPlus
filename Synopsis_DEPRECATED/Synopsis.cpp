//Title: Synopsis.cpp
//Description: Synopsis of a given trajectory of a moving object on a given grid.
//Platform(s): gcc (Ubuntu/Linaro 4.6.3-1ubuntu5) 4.6.3
//Date:     4/4/2015
//Revision: 4/3/2017

#include "Synopsis.h"

Synopsis::Synopsis(RegularGrid *_grid, Trajectory *_traj) 
{
    grid     = _grid;
	traj     = _traj;     
    duration =   0.0;                      
}

Synopsis::~Synopsis() 
{
    Qynopsis.clear();
    Cynopsis.clear();
    //cout << "Synopsis of " << traj->oid << " deleted" << endl;
}

void Synopsis::updateQ(Point *pnt)
{ 
    if ( Qynopsis.empty() )
    {
        duration += 1.0;
        synPartQ *part = new synPartQ();
        part->qid = pnt->qid;        //We passed the queryID
        part->t1  = pnt->t;
        part->t2  = pnt->t;  
        Qynopsis.push_back(*part); 
    }
    else if ( Qynopsis.back().qid != pnt->qid ) //object changed qid
    {
        duration += pnt->t - Qynopsis.back().t2 + 1.0;
        Qynopsis.back().t2 = pnt->t;
        synPartQ *part = new synPartQ();
        part->qid = pnt->qid;                   //We passed the queryID
        part->t1 = pnt->t;
        part->t2 = pnt->t;                          //When we first construct the synopsis , t1 = t2 , if we have only one element
        Qynopsis.push_back((*part));
    }
};

void Synopsis::updateC(Point *pnt)
{ 
    if ( Cynopsis.empty() )
    {
        duration += 1.0;
        synPartC *part = new synPartC();
        part->cid = this->grid->HashLocation(pnt->x, pnt->y); //Returns the cid of the current location update
        part->t1  = pnt->t;
        part->t2  = pnt->t;  
        Cynopsis.push_back(*part); 

        //---
        //Insert the oid to the History list in the Regular Grid
        // cout << "up_1" << endl; 
        this->grid->cell[part->cid].HisObjList.insert(pnt->oid);
        // cout << "up_2" << endl;
        //---
    }
    else //Synopsis exists
    {
        duration +=  pnt->t - Cynopsis.back().t2 + 1.0;
        Cynopsis.back().t2 = pnt->t;        
        unsigned int cur_cid = this->grid->HashLocation(pnt->x, pnt->y);
        if ( Cynopsis.back().cid != cur_cid )
        {
            synPartC *part = new synPartC();
            part->cid = cur_cid;         
            part->t1  = pnt->t;
            part->t2  = pnt->t;
            Cynopsis.push_back(*part); 

            //---
            //Insert the oid to the History list in the Regular Grid
            // cout << "up_11" << endl;
            this->grid->cell[part->cid].HisObjList.insert(pnt->oid);
            // cout << "up_22" << endl;
            //---
        }
    }
};

void Synopsis::cropQ(double t_rear)  //t_rear = t_NOW - maxWinRange
{
    while ( (!Qynopsis.empty()) && (Qynopsis.front().t2 < t_rear) ) 
    {          
        duration -= Qynopsis.front().getInt() + 1.0;
        Qynopsis.pop_front();                                              
    }

    if (!Qynopsis.empty())
    {
        if ( Qynopsis.front().t1 < t_rear )
        {
            duration -= t_rear - Qynopsis.front().t1 + 1.0;
            Qynopsis.front().t1 = t_rear;
        }
    }
}

void Synopsis::cropC(double t_rear)  //t_rear = t_NOW - maxWinRange
{
    while ( (!Cynopsis.empty()) && (Cynopsis.front().t2 < t_rear) ) 
    {          
        duration -= Cynopsis.front().getInt() + 1.0;

        //---
        //Erase the oid from the cell that the synopsis part belongs to and is going to be evicted
        // cout << "crop_1" << endl;
        this->grid->cell[Cynopsis.front().cid].HisObjList.erase(this->traj->oid);
        // cout << "crop_2" << endl;
        //---

        Cynopsis.pop_front();                                              
    }

    if (!Cynopsis.empty())
    {
        if ( Cynopsis.front().t1 < t_rear )
        {
            duration -= t_rear - Cynopsis.front().t1 + 1.0;
            Cynopsis.front().t1 = t_rear;
        }
    }
}

double Synopsis::scoringC()
{
    double score = 0.0;                   
    
    for ( auto i = Cynopsis.rbegin(); i != Cynopsis.rend(); i++ )
    {
        if ( this->grid->cell[i->cid].QryList.size() > 0 )
        {
            //The this part shall contribute something
            score += this->grid->cell[i->cid].QryList.size() * i->getInt();
        }
    }

    return score;                                           
}

int Synopsis::sizeQ()
{
    int size = Qynopsis.size();
    return size;                    
}

int Synopsis::sizeC()
{
    int size = Cynopsis.size();
    return size;                    
}

void Synopsis::print()
{
    cout << "--> Qynopsis <--" << endl;
    for ( auto i = Qynopsis.begin() ; i != Qynopsis.end() ; i++ )
        cout  << "(" << i->qid << ")" << "[" << i->t1 << ", " << i->t2 << "]" << "--->";
    cout << endl;

    cout << "--> Cynopsis <--" << endl;
    for ( auto i = Cynopsis.begin() ; i != Cynopsis.end() ; i++ )
        cout  << "(" << i->cid << ")" << "[" << i->t1 << ", " << i->t2 << "]" << "--->";
    cout << endl;
}

// double Synopsis::scoring(unsigned int boolWeightedScore)
// {
//     //boolWeightedScore = 0 -- No weight decay is applied
//     //boolWeightedScore = 1 -- Weight    decay is applied 
//     double score = 0.0;                   
//     double decWeight = 1.0/2.0;            
    
//     double normalizedTimeFraction = 1.0;
//     //list< synopart >::reverse_iterator i;

//     for ( auto i = Cynopsis.rbegin(); i != Cynopsis.rend(); i++ )
//     {
//         normalizedTimeFraction = ((double)i->getInt()) / duration;
//         if ( i->qid != 0 )
//             if ( boolWeightedScore == 0 )
//                 score += normalizedTimeFraction;
//             else if ( boolWeightedScore == 1 )
//                 score += decWeight * normalizedTimeFraction;

//         //Update decWeight for applying weight decaying
//         decWeight *= 1.0/2.0;                                   
//     }                                                       

//     return score;                                           
// }