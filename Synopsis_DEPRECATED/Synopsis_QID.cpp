//Title: Synopsis.cpp
//Description: Synopsis of a given trajectory of a moving object on a given grid.
//Platform(s): gcc (Ubuntu/Linaro 4.6.3-1ubuntu5) 4.6.3
//Date:     4/4/2015
//Revision: 4/3/2017

#include "Synopsis.h"

Synopsis::Synopsis(Trajectory *_traj, unsigned int _type) 
{
	traj     = _traj;     
    type     = _type;
    duration =   0.0;                      
}

Synopsis::~Synopsis() 
{
    Qynopsis.clear();
    Cynopsis.clear();
    cout << "Synopsis of " << traj->oid << " deleted" << endl;
}

double Synopsis::update(point *pnt, unordered_map< unsigned int, unsigned int > *U, int t, int s, unordered_map< unsigned int, unsigned int > *lastQID, unordered_set<unsigned int> *randomQ)
{
    unordered_multimap< unsigned int, synopart >::iterator iter;

    //Synopsis is empty -- does not exist
    if ( synopsis.empty() )
    {
        //New part of synopsis. Create and append it to the synopsis
        synopart *part = new synopart();
        //If the pnt belongs to a non-query then input in synopart qid = 0
        if ( ( pnt->qid == 0 ) || ( randomQ->find(pnt->qid) == randomQ->end() ) ) //And the remaining of the code is untouched!
        {
            part->qid = 0; //Insert it with qid zero
            //CAUTION: Changing the id of the non-zero BUT non-query points
            pnt->qid  = 0; //For main as well            
        }
        else
            part->qid = pnt->qid;              
        part->t1  = pnt->t;
        part->t2  = pnt->t;  
        //Update the duration     
        total_duration += 1.0;
        //Insert the pointers to list & hash
        synopsis.push_back( *part ); 
        synopsisHash.insert( {pnt->qid, *part} );

        return 1.0; //Return the duration of this object in that road in the time interval [t-r, t)
    }
    else //Synopsis exists
    {
        //Calculate the slide duration -- both for active and inactive trajectories
        double slide_duration = 0.0;
        if ( (pnt->t - synopsis.back().t2) <= s ) //CAUTION: Check again the condition!
        //if ( (t - synopsis.back().t2) <= s )
        {
            slide_duration = pnt->t - synopsis.back().t2 + 1.0;
            //cout << "1" << endl;            
        } 
        else
        {
            //cout << "2" << endl;
            slide_duration = pnt->t - (t-s) + 1.0; //CAUTION: pnt->t <> (t = t_NOW) . Not always equal. Batch stream processing.
        }

        //Update total duration
        total_duration +=  pnt->t - synopsis.back().t2 + 1.0;
        unsigned int temp_t2 = synopsis.back().t2;
        //Close the time interval for the last synopsis part
        synopsis.back().t2 = pnt->t;

        //Update the synopsisHash
        iter = synopsisHash.find( synopsis.back().qid );
        if ( iter != synopsisHash.end() )
        {
            if ( synopsisHash.count( synopsis.back().qid ) > 1 )
            {
                //cout << "More than one:\t" << synopsisHash.count( synopsis.back().qid ) << endl;
                //Find the specific synopart entry
                auto range = synopsisHash.equal_range( synopsis.back().qid );
                for ( auto it = range.first ; it != range.second ; it++ )
                    if ( ( fabs( it->second.t1 - synopsis.back().t1 ) <= 0.00001 ) && (  fabs( it->second.t2 - temp_t2 ) <= 0.00001 ) )
                    {
                        it->second.t2 = pnt->t;                        
                        break; //Entry found
                    }
            }
            else
            {
                //cout << "One element" << endl;
                iter->second.t2 = pnt->t;
            }
        }
        else
            cout << "SKATA UPDATE" << endl;

        //Distinguish 
        if ( synopsis.back().qid != pnt->qid )
        {
            //New part of synopsis. Create and append it to the synopsis
            synopart* part = new synopart();
            part->qid = pnt->qid;         
            part->t1  = pnt->t;
            part->t2  = pnt->t;
            //Insert the new synopart to list & hash
            synopsis.push_back( *part ); 
            synopsisHash.insert( {pnt->qid, *part} );

            //Change the qid of the lastQID entry
            //(*lastQID)[ traj->oid ] = pnt->qid;
        }

        return slide_duration;
    }

    return 0.0;
};

//Function that deletes the obsolete parts of the synopsis
void Synopsis::cropSynopsis(double t_rear)                                  //t_rear = t_NOW - maxWinRange
{
    unordered_multimap< unsigned int, synopart >::iterator iter;

    while ( (!synopsis.empty()) && (synopsis.front().t2 < t_rear ) ) 
    {          
        total_duration -= synopsis.front().getTimeDif() + 1.0;              //Update the duration of the synopsis

        //End of update of duration
        iter = synopsisHash.find( synopsis.front().qid );
        if ( iter != synopsisHash.end() )
        {
            if ( synopsisHash.count( synopsis.front().qid ) > 1 )
            {
                //cout << "More than one to DELETE:\t" << synopsisHash.count( synopsis.front().qid ) << endl;
                //Find the specific synopart entry
                auto range = synopsisHash.equal_range( synopsis.back().qid );
                for ( auto it = range.first ; it != range.second ; it++ )
                    if ( ( fabs(it->second.t2 - synopsis.front().t2 ) <= 0.00001 ) && ( fabs( it->second.t1 - synopsis.front().t1 ) <= 0.00001 ) )
                    {
                        synopsisHash.erase( it );
                        break; //Entry found & deleted
                    }
            }
            else
                if ( iter->second.t2 < t_rear )
                    synopsisHash.erase( iter );
        }
        // else
        //     cout << "SKATA DELETE while loop" << endl;

        //Delete the front synopart  
        synopsis.pop_front();                                              
    }

    if ( !synopsis.empty() )
    {
        if ( synopsis.front().t1 < t_rear )
        {
            total_duration -= t_rear - synopsis.front().t1 + 1.0;
            unsigned int temp_t1 = synopsis.front().t1;
            synopsis.front().t1 = t_rear;

            //End of update of duration
            iter = synopsisHash.find( synopsis.front().qid );
            if ( iter != synopsisHash.end() )
            {
                if ( synopsisHash.count( synopsis.front().qid ) > 1 )
                {
                    //cout << "More than one to DELETErefFront:\t" << synopsisHash.count( synopsis.front().qid ) << endl;
                    //Find the specific synopart entry
                    auto range = synopsisHash.equal_range( synopsis.back().qid );
                    for ( auto it = range.first ; it != range.second ; it++ )
                        if ( ( fabs( it->second.t2 - synopsis.front().t2 <= 0.00001 ) ) && ( fabs( it->second.t1 - temp_t1 ) <= 0.00001 )  )
                        {
                            it->second.t1 = t_rear;
                            break; //Entry found & deleted
                        }
                }
                else
                    iter->second.t1 = t_rear;
            }
            else
                cout << "SKATA DELETE front refinement" << endl;
        }
    }
}

//Finds the time interval at which the object was moving on q_id - O( log(synopsis.size()) )
double Synopsis::findTimeIntervalSynopsisFast(unsigned int q_id)
{
    double interval = -2.0;
    if ( synopsisHash.find( q_id ) != synopsisHash.end() )
    {
        if ( synopsisHash.count( q_id ) > 1 )
        { 
            interval = 0.0;
            //Find the specific synopart entry
            auto range = synopsisHash.equal_range( q_id );
            for ( auto it = range.first ; it != range.second ; it++ )
            {
                interval += it->second.getTimeDif();
            }
            //cout << "many" << endl;
        }
        else
        {
            interval = synopsisHash.find( q_id )->second.getTimeDif();
            //cout << "one" << endl;
        }
    }

    return interval;
}

//Assigns a score at the trajectory with decremented weights: 1/2, 1/4, 1/8, 1/16 ... etc etc
double Synopsis::scoring(unsigned int boolWeightedScore)
{
    //boolWeightedScore = 0 -- No weight decay is applied
    //boolWeightedScore = 1 -- Weight    decay is applied 
    double score = 0.0;                   
    double decWeight = 1.0/2.0;            
    
    double normalizedTimeFraction = 1.0;
    //list< synopart >::reverse_iterator i;

    for ( auto i = synopsis.rbegin(); i != synopsis.rend(); i++ )
    {
        normalizedTimeFraction = ((double)i->getTimeDif()) / total_duration;
        if ( i->qid != 0 )
            if ( boolWeightedScore == 0 )
                score += normalizedTimeFraction;
            else if ( boolWeightedScore == 1 )
                score += decWeight * normalizedTimeFraction;

        //Update decWeight for applying weight decaying
        decWeight *= 1.0/2.0;                                   
    }                                                       

    return score;                                           
}

// bool perform = ( (*it)->loc->qid != 0 ) && ( randomQ.find( (*it)->loc->qid ) != randomQ.end() ); //Only for queries in the center
//Calculates the TF score of the trajectory
double Synopsis::scoringTF(list< map < unsigned int, long double >* >::iterator acc)
{
    double score = 0.0;
    //list< synopart >::reverse_iterator i; //Easy to apply decay weights here!!!           

    for ( auto i = synopsis.rbegin(); i != synopsis.rend(); i++ )
        if ( i->qid != 0 ) //Every non-query contributes 0 to the synopsis score
            score += ((double)i->getTimeDif()) / (*(*acc))[ i->qid ];

    return score;
}

//Size of the synopsis in the interval [t_NOW - maxWinRange , t_NOW] - O(1)
int Synopsis::sizeSynopsis()
{
    int size = synopsisHash.size();
    return size;                    
}

//Prints the synopsis of this moving objects
void Synopsis::printSynopsis()
{
    //list< synopart >:: iterator i;
    cout << endl;
    for ( auto i = synopsis.begin() ; i != synopsis.end() ; i++ )
        cout  << i->qid << "[" << i->t1 << "," << i->t2 << "]" << " -> ";
    cout << endl;
}