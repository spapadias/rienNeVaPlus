//Title: main.cpp 
//Description: Main class of the Load Shedding Trajectory System
//Author: Serafeim Papadias
//Platform(s): gcc (Ubuntu/Linaro 4.6.3-1ubuntu5) 4.6.3
//Date:     04/04/2015
//Revision: 15/03/2017

#include "Trajectory.cpp"
#include "RegularGrid.cpp"
#include "Scan.cpp"
#include "parsingFunctions.cpp"
#include "winCPU.cpp"
#include "LoadShedder.cpp"
#include "linearRegression.cpp"
#include "ReadQueries.cpp"

using namespace std;

//Jaccard similarity 
double Jaccard(vector<int>)
{
    return 0.0;
}

//Global Variables
unordered_map<unsigned int, Trajectory* > T;  
unordered_map<unsigned int, unsigned int> U;            
unordered_map<unsigned int, unsigned int> lastQID;
unordered_map<unsigned int, MyTuple*> lastTupleOfTraj;
unordered_map<unsigned int, TQuery*> queryMap;
unordered_map<unsigned int, queryRes*> results;
unordered_map< unsigned int, unsigned int > invIndex;
Trajectory  *traj;
LoadShedder *loadS;

bool firstTime = true;

//Statistics of the program
unsigned int correctPIP = 0, totalPIP = 0; 
unsigned int cycles_total = 0, cycles_ls = 0;
unsigned int totalNumTrajSoFar = 0, active = 0, newbies = 0, profitable = 0;
unsigned int tup_total = 0, tup_processed = 0, tup_deleted = 0, tup_sofar = 0;
unsigned int boolW = 0;
int alpha = 0, delta = 0;

//LoadShedder *loadS;                                                       //A pointer to a LoadShedder instance
double L = 0.0, L_expected = 0.0, L_naive = 0.0;                            //The actual load of our system
double Lnaive_n_1_DIF_L_actual_n = 0.0; 
double tempLnaiveDIFF = 0.0;
bool belowCplusTheta  = true;       
bool firstTimeRunningLoadShedding = true;                                   //Bool for the first time load shedding

//Time Measurement in CPU cycles
double upStart, upEnd, evStart, evEnd, lsStart, lsEnd;
double upSum = 0.0, evSum = 0.0, lsSum = 0.0;
unsigned int upCount = 1, evCount = 1, lsCount = 1;
double cur_up_load = 0.0, cur_ev_load = 0.0, cur_ls_load = 0.0;             //Loads for update, evaluation and load shedding for each timestamp

int main(int argc, char *argv[])
{
    if (argc != 16)  
    {
        //EXAMPLE:  ./name_of_executable VALID obj3.txt 3 0 10 5 0 10 10 
        cout << "Usage: " << argv[0] <<     "[mode]"                              << 
                                            "[{real | synthetic}]"                <<
                                            "[input-file]"                        << 
                                            "{[stream-rate] | [timestamp-attr]}"  << 
                                            "[start-time]"                        <<
                                            "[query-file]"                        << 
                                            "[randomQ-file]"                      <<
                                            "[Granularity_X]"                     <<
                                            "[Granularity_Y]"                     <<
                                            "[Capacity(C)_code]"                  <<
                                            "[obj_count_expected]"                <<
                                            "[range]"                             <<
                                            "[slide]"                             <<
                                            "[LS_scheme]"                         <<
                                            "[operation {tim | acc}]"             << endl;

        exit(0);
    }

    cout.setf(ios::boolalpha);                                              
    char *modeTime;                
    modeTime = argv[1]; 

    //If equal to 0 => syntetic else if equal to 1 => real
    unsigned int isRealDataset = atoi(argv[2]);

    char *fileName;                
    fileName = argv[3];                                                     

    bool mode = (strcmp(modeTime, "VALID") == 0);                           
    //mode = 0 for VALID timestamps , else mode = 1 for TRANS timestamps
    //Scan operator over the incoming stream file (with object locations))
    Scan * scanStream = new Scan(fileName, mode);
    //Depending on the mode of reading, the fourth argument specifies ...
    if (mode)
        scanStream->setTimeAttribute(atoi(argv[4]));                        
    else
        scanStream->setRate(atoi(argv[4]));                                
    //Specifies the time when the window is being applied
    int t0 = atoi(argv[5]);

    char *QueryFileName;                
    QueryFileName = argv[6];

    char * randomQfile;
    randomQfile = argv[7];

    //Read the amount of incoming tuples per execution cycle -- from NO_LS experiment
    fstream finTup;
    finTup.open( randomQfile, ios::in );
    //Hash table which keeps the number of incoming tuples for each execution cycle
    unordered_set<unsigned int> randomQ;
    unsigned int temp_q_id;
    //cout << "Size:\t" << randomQ.size() << endl;
    
    while ( finTup >> temp_q_id )
        randomQ.insert( temp_q_id );
    finTup.close();

    //Granularity at X-axis
    unsigned int gx = atoi(argv[8]);                                        //Number of grid cells along x-dimension
    //Granularity at Y-axis
    unsigned int gy = atoi(argv[9]);                                        //Number of grid cells along y-dimension 

    scanStream->curTime = t0;

    /*const*/ unsigned int C_code = atoi(argv[10]);

    //Total number of moving objects that we expect to serve
    unsigned int obj_cnt = atoi(argv[11]);

    //Global Sliding Window settings
    int global_range = atoi(argv[12]);
    int global_slide = atoi(argv[13]);

    //Capacity in terms of the number of trajectories to process
    const unsigned int C = (unsigned int)( (C_code / 100.0) * (double)obj_cnt ); //e.g. 50% * 100000 = 50000 trajectories at most!

    //LS schema to use
    unsigned int schemaLS = atoi(argv[14]);

    //Type of operation
    //1 -- time experiment, 2 -- accuracy experiment
    unsigned int op = atoi(argv[15]);

    //Create grid partitioning
    box athensBox = box(450000 , 4175000,  500000, 4225000);     //Parameter box for Athens' synthetic dataset
    box romeBox   = box(1700000, 4500000, 1900000, 4700000);     //Parameter box for Rome's real dataset
    RegularGrid *grid;
    if ( isRealDataset == 1 )
        grid = new RegularGrid( gx, gy, romeBox );
    else
        grid = new RegularGrid( gx, gy, athensBox );
    //Allocate the space
    if ( !grid->Allocate() ) 
    {
        cout << "Memory Allocation Problem -- Grid partitioning!" << endl; 
        return 1;
    }

    //Query Reading , Query Hashing to the grid
    queryMap.clear();                                                                            
    ReadQueries *qRead = new ReadQueries(&global_range, &global_slide, &queryMap, grid, &randomQ);
    if ( isRealDataset == 1 )
        qRead->queryRead(QueryFileName, 1);                     
    else
        qRead->queryRead(QueryFileName, 0);                     

    //cout << "#queries:\t" << queryMap.size() << endl;
    //grid->printGridState();
    //exit(0); //Terminate the program

    vector<MyTuple*> inTuples;  

    //Initialize U -- CAUTION: For all expected queries 
    for (unsigned int i = 1; i < (obj_cnt + 1); i++ )
    {
        U.insert( {i, 1} );
        //lastTupleOfTraj[ i ];
    }

    //Time, Loads
    if ( op == 1 )
    {
        cout << "t" << ";" << "L_actual" << ";" << "C" << ";" <<
        //newby, profitable, active_traj, total_traj 
        "#active" << ";" << "#newbies" << ";" << "#profitable" << ";" << "#total" << ";" << 
        //#tuple_processed, #tuples total
        "#tup_processed" << ";" << "#tup_total" << ";" << "#tup_sofar" << ";" << "#tup_deleted" << ";" << 
        //#UNLS, #OVLS, #totalLS, #cyclesTotal 
        "#cyclesLS" << ";" << "#cyclesTotal" << ";" <<  "LS:_si_o_no?" << ";" << 
        //#queries_answered, case_of_queries 
        "#queries_ans" << ";" <<
        //Times in CPU cycles
        "up_time" << ";" << "ev_time" <<  ";" << "ls_time" << endl;
    }

    //Global clock -- t_NOW
    int t = 0;

    //STREAM INPUT - Keep processing data file until it gets exhausted
    while (scanStream->exhausted == false)
    {

        // unsigned int ones = 0;
        // for ( auto i = U.begin() ; i != U.end() ; i++ )
        //     if ( i->second == 1 )
        //         ones++;
        // cout << "U.size:\t" << ones << endl; 

        //Read streaming data -- Batch until timestamp t
        inTuples = scanStream->consumeInput(t);          

        //----------------------------------------------------------------------------------------------
        //UPDATE PHASE: Create new tuples for the current timestamp value ------------------------------
        //----------------------------------------------------------------------------------------------
        upStart = get_cpu_time();

        cycles_total++;
        newbies       = 0;
        tup_total     = 0;
        tup_processed = 0;
        tup_deleted   = 0;

        for ( auto it = inTuples.begin(); it != inTuples.end(); it++ )
        {   
            tup_total++;

            auto itT = T.find( (*it)->oid );
            if ( itT == T.end() )
            {
                totalNumTrajSoFar++;
                newbies++;
                tup_processed++;

                traj = new Trajectory( (*it)->oid, t, grid );

                //---------------------------------------------------------------------------------
                int retPIP = traj->findQid((*it)->loc, 1, &lastTupleOfTraj);
                if ( retPIP != -1000 )
                {
                    if ( retPIP == 1 )
                        correctPIP++;
                    totalPIP++;
                }
                
                //Adjust the qid of non-query tuples with non-zero qid to 0
                if ( ((*it)->loc->qid == 0) || (randomQ.find( (*it)->loc->qid ) == randomQ.end()) )
                    (*it)->loc->qid = 0;
                //----------------------------------------------------------------------------------

                //Velocity calculation
                //traj->findVelocity( (*it)->loc, lastTupleOfTraj[(*it)->oid]->loc, 1 );

                traj->insertPoint( (*it)->loc, &U, t, 0 );
                T.insert( { (*it)->oid, traj } ); //Create trajectory even if first tuples qid = 0
            }
            else 
            {
                if ( U[ (*it)->oid ] == 1 )
                {
                    tup_processed++;

                    //---------------------------------------------------------------------------------
                    int retPIP = itT->second->findQid((*it)->loc, 2, &lastTupleOfTraj); 
                    if ( retPIP != -1000 )
                    {
                        if ( retPIP == 1 )
                            correctPIP++;
                        totalPIP++;
                    }
                 
                    //Adjust the qid of non-query tuples with non-zero qid to 0
                    if ( ((*it)->loc->qid == 0) || (randomQ.find( (*it)->loc->qid ) == randomQ.end()) )
                    {
                        (*it)->loc->qid = 0;

                        //Even if we do NOT include this location update into trajectory -- We update the timeLastUpdate
                        if ( T.find((*it)->oid)!= T.end() )
                            T[ (*it)->oid ]->lastUpdateTime = (*it)->loc->t;
                        //Maintaint the lastQID & the last tuple for each object
                        lastQID[ (*it)->oid ] = (*it)->loc->qid;
                        lastTupleOfTraj[ (*it)->oid ] = (*it);
                        
                        continue;               //Proceed to next tuple
                    }
                    //---------------------------------------------------------------------------------

                    //------------------------------------------------------------------------------------------------------
                    //if ( itT->second->seq->segPoints.size() > 2 )
                        if ( itT->second->findVelocity( (*it)->loc, lastTupleOfTraj[(*it)->oid]->loc, (*it)->oid ) == -1 )
                        {
                            //Even if we do NOT include this location update into trajectory -- We update the timeLastUpdate
                            if ( T.find((*it)->oid)!= T.end() )
                                T[ (*it)->oid ]->lastUpdateTime = (*it)->loc->t;
                            //Maintaint the lastQID & the last tuple for each object
                            lastQID[ (*it)->oid ] = (*it)->loc->qid;
                            lastTupleOfTraj[ (*it)->oid ] = (*it);

                            continue;  
                        }
                    //------------------------------------------------------------------------------------------------------

                    itT->second->insertPoint( (*it)->loc, &U, t, 1 );
                }
            }

            //If traj exists -- update the last time that the current object sent a  location update
            if ( T.find((*it)->oid)!= T.end() )
                T[ (*it)->oid ]->lastUpdateTime = (*it)->loc->t;
            //Maintaint the lastQID & the last tuple for each object
            lastQID[ (*it)->oid ] = (*it)->loc->qid;
            lastTupleOfTraj[ (*it)->oid ] = (*it);
        }

        //Calculate the number of active trajectories -- the ones that are still sending location updates
        active = 0;
        for ( auto itT = T.begin() ; itT != T.end() ; itT++ )
            if ( itT->second->lastUpdateTime >= (t - (double) (1.0 * global_slide)) )
                active++;
        
        upEnd = get_cpu_time();

        //Times per Execution Cycle
        cur_up_load = upEnd - upStart;

        if ( active > C )
        {  
            //Total update time
            upSum += cur_up_load;                     //Update Time of this cycle in CPU cycles
            upCount++;                                //Increase update counter
        }

        //---------------------------------------------------------------------------------------------------------------
        //EVALUATION PHASE ----------------------------------------------------------------------------------------------
        //---------------------------------------------------------------------------------------------------------------
        evStart = get_cpu_time();

        //Condition for query answering
        if ( (t >= global_range) && (t % global_slide == 0) )
        {
            //Initialize the inverted index structure -- O(N)
            invIndex.clear();
            for ( auto i = T.begin(); i != T.end(); i++ )
                invIndex.insert( {i->first, 0} );

            //Initialize results map structure
            results.clear();
            //Time Complexity -- O(M)
            for ( auto i = queryMap.begin() ; i != queryMap.end() ; i++ )
            {
                queryRes *res = new queryRes(i->first);
                results.insert( {i->first, res} );
            }

            //Time Complexity -- O(N x trajectory.size() )
            for ( auto i = T.begin() ; i != T.end() ; i++ )
            {
                //Output the number of deleted location updates in this execution cycle
                tup_deleted += i->second->check4ExpiringPoints( t - (double) global_range );

                //Partial results for this trajectory
                unordered_map<unsigned int, queryRes*> partialResults;
                partialResults.clear();

                //Scan this trajectory and update the instant Vel aggregates
                for ( auto j = i->second->seq->segPoints.rbegin() ; j != i->second->seq->segPoints.rend() ; j++ )
                {
                    if ( std::isnan((*j)->v) )
                    {
                        continue; //Ingnore the filthy location updates
                        //cout << "malakia" << endl;
                    }

                    //Include only non-zero qids to the result
                    if ( (*j)->qid != 0 )
                    {
                        //First time encoutering this qid while traversing the trajectory
                        if ( partialResults.find( (*j)->qid ) == partialResults.end() )
                        {
                            //Increase the #queries that this trajectory contributes
                            invIndex[ i->first ]++;

                            queryRes *partial = new queryRes( (*j)->qid );

                            //cout << std::fixed;
                            //cout << setprecision(3);
                            //cout << "oid:\t" << (partial->qid == (*j)->qid) << " sum:\t" << partial->sum << " count:\t" << partial->count << endl; 

                            partial->sum += (*j)->v;
                            partial->count++;

                            //Insert the new partial result into the map
                            partialResults.insert( {(*j)->qid, partial} );    
    //                            cout << "qid:\t" << partial->qid /*(partial->qid == (*j)->qid)*/ << " sum:\t" << partial->sum << " count:\t" << partial->count << endl; 
                        }
                        else
                        {
                            //Update the existing partial result
                            partialResults[ (*j)->qid ]->sum += (*j)->v;
                            partialResults[ (*j)->qid ]->count++;
                        }
                    }
                    // else
                    //     cout << "proto stigma -- oid:\t" << i->first << " qid:\t" << (*j)->qid << endl;
                }

                //partialResults now contain all the contribution of this trajectory -- Insert them into results
                for ( auto j = partialResults.begin() ; j != partialResults.end() ; j++ )
                {
                    if ( queryMap.find( j->first ) != queryMap.end() )
                        if ( j->second->count > 1 ) //Already initialized final results for this execution cycle
                        {
                            //Here we insert the average velocities for each trajectory in each query
                            results[ j->first ]->sum += j->second->sum / (double)(j->second->count);

                            results[ j->first ]->count++;
                            results[ j->first ]->oids.insert( i->first ); //Include this objects id to the set of oids contributing to this qid
                        }
                        // else
                        //     cout << "COUNTER == 0" << endl;
                }
              
            }

            //Now output the final velocity results by traversing the results map
            for ( auto i = results.begin() ; i != results.end() ; i++ )
            {
                if      ( op == 1 )
                    ; //Do nothing with speed measurements
                else if ( op == 2 )
                    if ( i->second->count == 0 )
                        cout << t << "," << i->first << ",-1000," << i->second->count << endl;
                    else  
                        cout << t << "," << i->first << "," << i->second->sum / (double) i->second->count << "," << i->second->count << endl;
            }

        }

        //Calculate the number of profitable trajectories in this execution cycle
        profitable = 0;
        for ( auto i = invIndex.begin(); i != invIndex.end(); i++ )
            if ( invIndex[ i->first ] > 0 )
                profitable++;            

        //Average evaluation time measurement
        evEnd = get_cpu_time();

        //Evaluation Time for this execution cycle
        cur_ev_load = evEnd - evStart;

        if ( active > C )
        {  
            evSum += cur_ev_load;         //Total Sum of Evaluation Time in CPU cycles
            evCount++;                    //Increase evaluation counter
        }

        //------------------------------------------------------------------------------------------------
        //LOAD SHEDDING PHASE ----------------------------------------------------------------------------
        //------------------------------------------------------------------------------------------------
        //lsStart = get_cpu_time(); 

        //Calculate the load L of this timestamp - t
        L = cur_up_load + cur_ev_load + cur_ls_load;

        //Number of trajectories to keep in next execution cycle
        alpha = C;

        bool happenedLS = false;
        //NO LS condition!
        if ( schemaLS != 0 )
        {
            if ( active > C )         
            {   
                happenedLS = true;

                cycles_ls++;

                //Initialize the bitmap U
                for (unsigned int i = 1; i < (obj_cnt + 1); i++ )
                    U[ i ] = 0;

                loadS = new LoadShedder( &T, &U, &invIndex, &results, &alpha, boolW ); //CAUTION: We might need to move this one to LS section later
                lsStart = get_cpu_time();
                //------------------------------- 
                int retValue = 0;
                if      ( schemaLS == 1 )
                    retValue = loadS->RandomLS(); 
                else if ( schemaLS == 2 )                          
                    retValue = loadS->TF();
                //-------------------------------
                lsEnd = get_cpu_time();
            }
            // //CAUTION: This might be bad for initial model in the paper. CPU time.
            // else
            // {
            //     //Fix the bitmap U -- If no LS happening use all available trajectories
            //     for (unsigned int i = 1; i < (obj_cnt + 1); i++ )
            //         U[ i ] = 1;
            // }
        }  

        //lsEnd = get_cpu_time();

        if ( schemaLS != 0 )
            cur_ls_load = lsEnd - lsStart;
        else
            cur_ls_load = 0.0;

        if ( active > C )
        {  
            lsSum += cur_ls_load;
            lsCount++;
        }

        //Resultados -- Time, Loads
        if ( op == 1 )
        {
            cout << std::fixed << t << ";" << setprecision(3) << L << ";" << setprecision(3) << C << ";" <<
            //newby, profitable, active_traj, total_traj 
            active << ";" << newbies << ";" << profitable << ";" << T.size() << ";" <<
            //#tuple_processed, #tuples total
            tup_processed << ";" << tup_total << ";" << tup_sofar << ";" << tup_deleted << ";" << 
            //#UNLS, #OVLS, #totalLS, #cyclesTotal 
            cycles_ls << ";" << cycles_total << ";" << happenedLS << ";" <<
            //#queries_answered, case_of_queries 
            results.size() << ";" <<
            //Times in CPU cycles
            setprecision(3) << cur_up_load << ";" << setprecision(3) <<  cur_ev_load << ";" << setprecision(3) << cur_ls_load << endl;
        }

        inTuples.clear(); //Clear the tuples you read in the previous step

        //Maintain the total number of tuples so far in the experiment
        tup_sofar += tup_total;

        t = t + global_slide; //Here δ = global_slide

    } //End of while loop - input file exhausted
  
    cout << endl;
    cout << "**********************" << endl;    
    cout << scanStream->recCount << " records processed in total." << endl;

    cout << "Execution completed successfully!" << endl;

    double up_time = upSum / ((double)(upCount - 1));
    double ev_time = evSum / ((double)(evCount - 1));
    double ls_time = lsSum / ((double)(lsCount - 1));
    double total_time = up_time + ev_time + ls_time;

    cout << endl;
    cout << "*** Time Statistics in CPU cycles ***" << endl;
    cout << "Average Update time    : " << up_time << " CPU cycles." << endl;
    cout << "Average Evaluation time: " << ev_time << " CPU cycles." << endl;
    cout << "Average LS time        : " << ls_time << " CPU cycles." << endl;
    //Total percentage
    cout << "Average Update Percentage    : " << (100 * up_time) / ( total_time ) << " %" << endl;
    cout << "Average Evaluation Percentage: " << (100 * ev_time) / ( total_time ) << " %" << endl;
    cout << "Average LS Percentage        : " << (100 * ls_time) / ( total_time ) << " %" << endl;

    cout << "REMARK\nStatistics of efficiency where gathered only from execution cycles where LS happened. #executionCyclesLS:\t" << (upCount - 1) << endl;
    cout << "up_count = ev_count = ls_count:\t" << ( (upCount == lsCount) && ( upCount == evCount ) ) << endl;

    cout << "Percentage of correct PIP executions:\t" << (correctPIP / (double) totalPIP) * 100 << endl;

    //Memory Deallocation -------------------------
    delete scanStream;
    //    delete grid;
    for ( auto i = T.begin() ; i != T.end() ; i++ )
        delete i->second;
    T.clear();
    U.clear();
    queryMap.clear();
    lastQID.clear();
    lastTupleOfTraj.clear();
    invIndex.clear();
    //----------------------------------------------

    return 0;
}