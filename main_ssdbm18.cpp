//Title: main.cpp 
//Description: Main class of the Load Shedding Trajectory System
//Author: Serafeim Papadias
//Platform(s): gcc (Ubuntu/Linaro 4.6.3-1ubuntu5) 4.6.3
//Date:     04/04/2015
//Revision: 20/01/18 --> NEW Version of code

#include "./Utilities/Structures.h"             
#include "./Utilities/Functions.cpp"
#include "./Utilities/winCPU.cpp"
#include "./Utilities/LinearRegression.cpp"
#include "./Point/Point.cpp"
#include "./RegularGrid/RegularGrid.cpp"
#include "./ReadQueries/ReadQueries.cpp"
#include "./Scan/Scan.cpp"
#include "./Trajectory/Trajectory.cpp"
#include "./LoadShedder/LoadShedder.cpp"
#include "./Synopsis_DEPRECATED/Synopsis.cpp"

using namespace std;

int main(int argc, char *argv[])
{
    if (argc != 25)  
    {
        cout << "Usage: " << argv[0] <<     "[mode]"                              <<
    //--- Dataset (object + query file) 
                                            "[{real | synthetic}]"                <<
                                            "[obj-count]"                         <<
                                            "[object-file]"                       <<
                                            "[query-count]"                       <<
                                            "[query-file]"                        <<
                                            "[ranQuery-file]"                     <<
                                            "{[stream-rate] | [timestamp-attr]}"  << 
                                            "[start-time]"                        <<
    //--- Grid Setting
                                            "[Granularity_X]"                     <<
                                            "[Granularity_Y]"                     <<    
    //--- Load Shedding Setting
                                            "[Lpeak]"                             <<
                                            "[Capacity(C)_code]"                  <<
                                            "[LS_scheme]"                         <<
                                            "[thresLS]"                           <<
    //--- Linear Regression Setting                      
                                            "[prediction]"                        <<
                                            "[num_Loads]"                         <<
    //--- Sliding Window Setting                
                                            "[{global | variableWindows}]"        <<
                                            "[range]"                             <<
                                            "[slide]"                             <<
    //--- System Setting
                                            "δ-period"                            <<
    //--- Type of Query (AVG / MAX)         
                                            "{AVG | MAX}"                         <<
    //--- Type of operation (time / accuracy)
                                            "[operation {tim | acc}]"             <<
                                            "[firstCritPolicy]"                   << endl;
        exit(0);
    }

    cout.setf(ios::boolalpha); char *modeTime; modeTime = argv[1]; 
    unsigned int isRealDataset = atoi(argv[2]);     //If equal to 0 => syntetic else if equal to 1 => real
    unsigned int obj_cnt = atoi(argv[3]);
    char *objectFile; objectFile = argv[4]; 
    unsigned int qry_cnt = atoi(argv[5]);           //Number of queries to define
    char *queryFile; queryFile = argv[6];
    char * randomQfile; randomQfile = argv[7];
    //-----------------------------
    bool mode = (strcmp(modeTime, "VALID") == 0);   //mode = 0 for VALID timestamps, else mode = 1 for TRANS timestamps
    Scan *scanStream = new Scan(objectFile, mode);
    if (mode)
        scanStream->setTimeAttribute(atoi(argv[8]));                       
     else
        scanStream->setRate(atoi(argv[8]));   
    int t0 = atoi(argv[9]);                         //Time when the window starts being applied
    scanStream->curTime = t0;
    //-----------------------------
    unsigned int gx = atoi(argv[10]);
    unsigned int gy = atoi(argv[11]);
    const double Lpeak = atof(argv[12]);        
    unsigned int C_code = atoi(argv[13]);       //10, 25, 33, 50, 67, 80% Lpeak
    //Pick a Capacity (w.r.t. CPU cycles or w.r.t. #trajectories)
    const double C = (double)( (double)C_code * Lpeak ) / 100.0;
    // const unsigned int C = (unsigned int)( (C_code / 100.0) * (double)obj_cnt ); //e.g. 50% * 100000 = 50000 trajectories at most!
    unsigned int schemaLS = atoi(argv[14]);
    double thresLS = atoi(argv[15]); double theta = ( (double)thresLS * Lpeak ) / 100.0;
    unsigned int prediction = atoi(argv[16]); //Default 1 --> linear regression, 2 --> zero-knowledge model (Lpred = L)
    unsigned int numLoads = atoi(argv[17]); //#Loads to take under consideration for the prediction
    unsigned int isGlobalWin = atoi(argv[18]); // 1 if we use global window, 0 if each query has its own window
    unsigned int global_range = atoi(argv[19]);
    unsigned int global_slide = atoi(argv[20]);
    unsigned int period = atoi(argv[21]);
    unsigned int queryType = atoi(argv[22]); //1 ==> AVG , 2 ==> MAX
    unsigned int op = atoi(argv[23]); //1 -- time experiment, 2 -- accuracy experiment
    unsigned int firstCritPolicy = atoi(argv[24]);

    //Create the randomQ structure with the qids that need answer
    fstream finTup; finTup.open( randomQfile, ios::in );
    vector<unsigned int> randomQ;
    unsigned int temp_qid;
//    cout << "ranQ_before:\t" << randomQ.size() << endl;   
    while (finTup >> temp_qid)
        randomQ.push_back(temp_qid);
    finTup.close();

    // // Here maintain the qids randomly
    // auto rng = std::default_random_engine {};
    // shuffle(begin(randomQ), end(randomQ), rng);

    // // Keep only qry_cnt queries at random
    // while (randomQ.size() > (unsigned long)qry_cnt)
    //     randomQ.pop_back();

//    cout << "ranQ__after:\t" << randomQ.size() << endl;
    
    //---------------------------------------------------------------------------------------------------------------------------------//
    //------------------ Global variables ---------------------------------------------------------------------------------------------//
    //---------------------------------------------------------------------------------------------------------------------------------//
    Trajectory  *traj  = NULL;
    LoadShedder *loadS = NULL;    
    vector<Point*> inTuples;
    //------------------------//
    //------------------------//
    //This must be unsigned int
    unsigned int t = 0;       // Global Clock -- t_now
    //------------------------//
    unsigned int notAnumber = 0;
    //------------------------//
    unordered_map<unsigned int, Point* >      O; //Last Location Update of an Object
    unordered_map<unsigned int, Trajectory* > T; //Trajectory of Objects
    unordered_map<unsigned int, unsigned int> U;            
    unordered_map<unsigned int, TQuery*>      Q;
    //unordered_map<unsigned int, unsigned int> lastQID;
    unordered_map<unsigned int, Point*>     lastTuple;
    unordered_map<unsigned int, Point*>     preLastTuple;
    unordered_map<unsigned int, queryRes*>    results;
    unordered_map<unsigned int, unsigned int> invIndex;
    map<int, double> trendData; //CAUTION: Do not use unordered_map here. The elements are not stored sequentialy by key
                                //The elements are not stored sequentialy by key BUT we want them to be so
    
    //PIP Statistics
    unsigned int correctPIP = 0, totalPIP = 0;    
    //--- Simulate PIP -----------------
    double sumPIPtime = 0.0;
    unsigned int countPIP = 0;
    double pipStart = 0.0, pipEnd = 0.0;
    //----------------------------------

    //Statistics of the program
    unsigned int cycles_total = 0, cycles_ls = 0, cycles_over = 0, cycles_under = 0;
    unsigned int active = 0, newbies = 0, profitable = 0, traj_sofar = 0;
    unsigned int tup_total = 0, tup_processed = 0, tup_deleted = 0, tup_sofar = 0;
    //unsigned int boolW = 0;
    
    //Load Shedding 
    double L = 0.0, L_pred = 0.0;                            //The actual load of the system & the predicted one
    bool firstTimeLS = true, activeUnderLS = false;
    int traj_used = 0, n = 0, delta = 0;                     //#traj used in this execution cycle & n - number of traj to be used
    vector<double> x_reg, y_reg;      
    int  retValue = 0;
    bool happenedLS = false;      
    unsigned int numQans = 0;
    unsigned int tempNumOfLoads4pred = 0;

    // Useful for statistics
    double maxLoad = 0.0;                    

    //Time Measurement in CPU cycles
    double upStart = 0.0, upEnd = 0.0, evStart = 0.0, evEnd = 0.0, lsStart = 0.0, lsEnd = 0.0;
    double upSum = 0.0, evSum = 0.0, lsSum = 0.0;
    unsigned int upCount = 1, evCount = 1, lsCount = 1;
    double cur_up_load = 0.0, cur_ev_load = 0.0, cur_ls_load = 0.0; //Loads for update, evaluation and load shedding for each timestamp
    //Time Measurement in MILliseconds
    double upStartMIL = 0.0, upEndMIL = 0.0, evStartMIL = 0.0, evEndMIL = 0.0, lsStartMIL = 0.0, lsEndMIL = 0.0;
    double upSumMIL = 0.0, evSumMIL = 0.0, lsSumMIL = 0.0;
    double cur_up_loadMIL = 0.0, cur_ev_loadMIL = 0.0, cur_ls_loadMIL = 0.0; //Loads for update, evaluation and load shedding for each timestamp

    //---------------------------------------------------------------------------------------------------------------------------------//
    //------------------ Create grid partiioning --------------------------------------------------------------------------------------//
    //---------------------------------------------------------------------------------------------------------------------------------//
    box athensBox = box(450000 , 4175000,  500000, 4225000);  //Parameter box for Athens' synthetic dataset
    box romeBox   = box(1757000, 4608800, 1822000, 4673800);  //Parameter box for Rome's real dataset
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
    
// cout << "Grid Patritioning OK" << endl;

    //---------------------------------------------------------------------------------------------------------------------------------//
    //------------------ Query Reading , Query Hashing to the grid --------------------------------------------------------------------//
    //---------------------------------------------------------------------------------------------------------------------------------//
    Q.clear();                                                                            
    ReadQueries *qRead = new ReadQueries(isGlobalWin, global_range, global_slide, &Q, grid, &randomQ);
    if ( isRealDataset == 1 )
        qRead->queryRead(queryFile, 1);                     
    else
        qRead->queryRead(queryFile, 0);
    //Clear the randomQ set structure
    randomQ.clear(); 
    //-----------------------------------------------------------------------------------//

//    cout << "#queries:\t" << Q.size() << endl;
//    grid->printGridState();

//exit(0);

// cout << "Query Reading and Query Hashing OK" << endl;
    
    //---------------------------------------------------------------------------------------------------------------------------------//
    //------------------ Initialize U -- CAUTION: For all expected queries ------------------------------------------------------------//
    //---------------------------------------------------------------------------------------------------------------------------------//
    for (unsigned int i = 1; i < (obj_cnt + 10); i++ )
    {
        U.insert( {i, 1} );
        lastTuple[i];
    //    preLastTuple[i];
    }

// cout << "U initiallization OK" << endl;
//exit(0);

    //---------------------------------------------------------------------------------------------------------------------------------//
    //------------------ Time, Loads --------------------------------------------------------------------------------------------------//
    //---------------------------------------------------------------------------------------------------------------------------------//
    if ( op == 1 )
    {
        cout << "t" << ";" << "L_actual" << ";" << "L_pred" << ";" << "C" << ";" << "C-theta" << ";" << "C+theta" << ";" <<
        //newby, profitable, active_traj, total_traj 
        "#active" << ";" << "#newbies" << ";" << "#profitable" << ";" << "#total" << ";" << 
        //traj used, difference, trajs to be used
        "n" << ";" << "s" <<";" << "n'" << ";" <<
        //#tuple_processed, #tuples total
        "#tup_processed" << ";" << "#tup_total" << ";" << "#tup_sofar" << ";" << "maintained" << ";" << "#tup_deleted" << ";" << 
        //#UNLS, #OVLS, #totalLS, #cyclesTotal 
        "cycles_under" << ";" << "cycles_over" << ";" << "#cycles_ls" << ";" << "#cyclesTotal" << ";" <<  "LS:_si_o_no?" << ";" << 
        //#queries_answered, case_of_queries 
        "#queries_ans" << ";" <<
        //Times in CPU cycles
        "up_time" << ";" << "ev_time" <<  ";" << "ls_time" << ";" <<
        //Times in CPU cycles
        "MIL_up_time" << ";" << "MIL_ev_time" <<  ";" << "MIL_ls_time" << ";" << "#Loads4pred" << endl;
    } 


// cout << "Preprocessing OK" << endl;
//exit(0);
    //auto rng = std::default_random_engine {};

    unsigned int maintained = 0; // Number of tuples maintained in total from the beginning 

    std::random_device rd;     // only used once to initialise (seed) engine
    std::mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)
    std::uniform_real_distribution<> dis(0.0, 1.0);

    //---------------------------------------------------------------------------------------------------------------------------------//
    //------------------ STREAM INPUT - Keep processing data file until it gets exhausted ---------------------------------------------//
    //---------------------------------------------------------------------------------------------------------------------------------//
    while ( scanStream->exhausted == false )
    {
// cout << "Initiating Execution Cycle" << endl;

        inTuples = scanStream->consumeInput(t); //Batch Processing until timestamp t

        //-----------------------------------------------------------------------------------------------------------------------------//
        // Trajectory Update phase ----------------------------------------------------------------------------------------------------//
        //-----------------------------------------------------------------------------------------------------------------------------//
        upStart = get_cpu_time(); upStartMIL = get_wall_time();
        //-----------------------------------------------------  
        cycles_total++;
        newbies       = 0;
        tup_total     = 0;
        tup_processed = 0;
        tup_deleted   = 0;

        // --- Gather the active tids for this execution cycle
        vector<unsigned int> active_tids;

        for (auto i = U.begin() ; i != U.end() ; i++)
            if (i->second == 1)
                active_tids.push_back(i->first); // Gather the trajectory ids with U[id]=1
        // --- --- ---

// cout << "Before tuple for loop" << endl;

        for ( auto i = inTuples.begin(); i != inTuples.end(); i++ )
        {  
            tup_total++;

            //***
            // RandomTupleLS --> Reservoir sampling on the actual incoming tuples --> not the trajectory ids!!!
            if ( schemaLS == 6 )
            {
                // CAUTION: Pilot testify
               if (!firstTimeLS)           // CAUTION!!!
                {
                    double acceptance_probability = (1.5 * n) / ((double) (tup_sofar + tup_total));
                    if (dis(rng) < acceptance_probability) // With probability 3n/tupleSoFar --> Process the tuple
                        maintained++;
                    else 
                        continue;                 
                }
            }  // ---
            //***

            auto j = T.find( (*i)->oid );           
            if ( j == T.end() )
            {
                newbies++;
                traj_sofar++;
                tup_processed++;

                traj = new Trajectory( grid, (*i)->oid, (*i)->t, &preLastTuple, &lastTuple, &Q, &U );

                //----------------------------------
                pipStart   = get_wall_time();
                //PIP
                int retPIP = traj->findQid( *i, 1 );
                //---
                pipEnd     = get_wall_time();
                sumPIPtime = pipEnd - pipStart;
                countPIP++;
                //----------------------------------
                if (retPIP != -1000)
                {
                    if (retPIP == 1)
                        correctPIP++;
                    totalPIP++;
                } //How many out of the <> 0 qids have been correctly classified?
                
                //Velocity
                //traj->findVelocity( *i, 1 ); -- We do not calculate the first location update's instant velocity
                //Insertion

            //----------------------------
                if (isRealDataset == 1)
                {
                    if (t < 4900) //BUG SOLVER
                        traj->insertPoint(*i);
                }
                else // artificial dataset
                {
                    if (t < 4100)
                        traj->insertPoint(*i);
                }
            //----------------------------
// cout << "test0" << endl;
                //Insertion traj
                T.insert( {traj->oid, traj} ); //Create trajectory even if first tuples qid = 0
// cout << "test1" << endl;
    //***

                // ReservoirSampleLS --> After processing the 1rst tuple of a new trajectory -- Reservoir its trajectory id!!!
                if ( schemaLS == 3 )
                {
                    // CAUTION: Pilot testify
                   if (!firstTimeLS)           // CAUTION!!!
                    {
                        loadS = new LoadShedder(grid, &T, &U, &invIndex, &results, &Q, &n, obj_cnt, -1, t, global_range); // isUnder=-1 --> means we do not care about that now
                        loadS->ReservoirSampleLS_tuple((*i)->oid, &active_tids);                            
                    }
                }  // ---
    //***
// cout << "test2" << endl;

            }
            else if ( U[(*i)->oid] == 1 )
            {
// cout << "test3" << endl;
               tup_processed++;

                //----------------------------------
                pipStart   = get_wall_time();
                //PIP
                int retPIP = traj->findQid( *i, 2 );
                //---
                pipEnd     = get_wall_time();
                sumPIPtime = pipEnd - pipStart;
                countPIP++;
// cout << "test4" << endl;
                //----------------------------------
                if (retPIP != -1000)
                {   
                    if (retPIP == 1)
                        correctPIP++;
                    totalPIP++;
                } //How many out of the <> 0 qids have been correctly classified?
                else
                {
                    T[(*i)->oid]->mrut   = (*i)->t;
        T[(*i)->oid]->synopsis->updateQ((*i));
        T[(*i)->oid]->synopsis->updateC((*i));
                    lastTuple[(*i)->oid] = (*i);
                    grid->updateObject(*i);
                    grid->ts             = (*i)->t;
                    continue; //Proceed to next tuple
                    //CAUTION: The location updates with qid == 0 are NOT included to the trajectory
                }
// cout << "test5" << endl;
                //Velocity
                if ( j->second->findVelocity( *i ) == -1 )
                {
                    T[(*i)->oid]->mrut   = (*i)->t;
        T[(*i)->oid]->synopsis->updateQ((*i));
        T[(*i)->oid]->synopsis->updateC((*i));
                    lastTuple[(*i)->oid] = (*i);
                    grid->updateObject(*i);
                    grid->ts             = (*i)->t;
                    continue; //Proceed to next tuple
                }
// cout << "test6" << endl;
                //Insertion
                j->second->insertPoint( *i );
// cout << "test7" << endl;
            }
            // ---
            // ---
            // ---
            // else // Case where oid is unimportant and all trajectories have been added
            // {
            //     if (schemaLS == 3)
            //     {
            //         if ((unsigned long)n >= T.size())
            //         {
            //             // Include this trajectory-id with some probability --> 1 / unimportant
            //             if ((unsigned long)n != active_tids.size())
            //             { 
            //                 double probability_thres = 1 / (n - active_tids.size());
            //                 double min = 0.0;
            //                 double max = 1.0;
            //                 // std::random_device rd;                           // only used once to initialise (seed) engine
            //                 // //std::mt19937 rnggg(rd());                        // random-number engine used (Mersenne-Twister in this case)
            //                 // std::default_random_engine re;
            //                 // std::uniform_int_distribution<double> unif(min,max); // guaranteed unbiased
            //                 // auto output = unif(re);

            //                 std::random_device rd;  //Will be used to obtain a seed for the random number engine
            //                 std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
            //                 std::uniform_real_distribution<> dis(min, max);
            //                 auto output = dis(gen);

            //                 if (output > probability_thres)
            //                     U[(*i)->oid] = 1; // Validate the trajectory with the aforementioned probability
            //             }  
            //         }
            //     }
            // }
            // ---
            // ---
            // ---
// cout << "test8" << endl;

            //Necessary Updates -------------
            T[(*i)->oid]->mrut   = (*i)->t;
            //----------------------------
            if (isRealDataset == 1)
            {
                if (t < 4900) //BUG SOLVER
                {
                    T[(*i)->oid]->synopsis->updateQ((*i));
                    T[(*i)->oid]->synopsis->updateC((*i));
                }
            }
            else // artificial dataset
            {
                if (t < 4100)
                {
                    T[(*i)->oid]->synopsis->updateQ((*i));
                    T[(*i)->oid]->synopsis->updateC((*i));
                }
            }
            //----------------------------
            lastTuple[(*i)->oid] = (*i);
            //Grid Update -------------------
            grid->updateObject(*i);
            grid->ts             = (*i)->t;
            //-------------------------------
// cout << "test9" << endl;
        
        //---------------------------------------    
        }

// cout << "After tuple for loop -- Almost DONE" << endl;        
        
        //Calculate the number of active trajectories -- O(N)
        active = 0;
        for ( auto i = T.begin() ; i != T.end() ; i++ )
            if ( i->second->mrut > ((double)t - (double)global_slide) ) // The objs that sent at least a location update in [t_now - range, t_now)
                active++;
        //-------------------------------------------------------------

        //-------------------------------------------------
        upEndMIL = get_wall_time(); upEnd = get_cpu_time(); 
        //-------------------------------------------------

        //Times per Execution Cycle
        cur_up_load = upEnd - upStart; cur_up_loadMIL = upEndMIL - upStartMIL; //Update Time of this cycle in CPU cycles
        //Total update time
        upSum      += cur_up_load;      //Total Update Time in CPU cycles
        upSumMIL   += cur_up_loadMIL;   //Total Update Time in MILliseconds
        upCount++;           

// -->
        // if (T.find(666) != T.end())
        //     T[666]->synopsis->print(); 

    // grid->printGridState();
    // cout << "cycle:\t" << t << endl;        
    // inTuples.clear();
    // tup_sofar += tup_total;
    // t = t + period;
    // continue;
    // cout << "Barrier" << endl;

// cout << "Update Phase OK -- " << t << endl;        

        //-----------------------------------------------------------------------------------------------------------------------------//
        // Query Evaluation phase -----------------------------------------------------------------------------------------------------//
        //-----------------------------------------------------------------------------------------------------------------------------//
        evStart = get_cpu_time(); evStartMIL = get_wall_time();
        //-----------------------------------------------------

        if ( true ) //Condition
        {
// cout << "0" << endl;

            //Initialize the inverted index structure -- O(N)
            invIndex.clear();
            for ( auto i = T.begin(); i != T.end(); i++ )
                invIndex.insert( {i->first, 0} ); //(oid , #queries_involved_in_cur_tim_win)
            //------------------------------------------------------------------------------

            //Memory Clearance ----------------------------------------
            for ( auto i = results.begin() ; i != results.end() ; i++ )
                delete i->second;
            results.clear(); //----------------------------------------
            //---------------------------------------------------------
// cout << "1" << endl;
            //Time Complexity -- O(M)
            for ( auto i = Q.begin() ; i != Q.end() ; i++ )
            {
                queryRes *res = new queryRes(i->first);
                results.insert( {i->first, res} );
            } //---------------------------------------------
// cout << "2" << endl;

            //Time Complexity -- O(N x trajectory.size_in_cur_window() ) ----------------
            for ( auto i = T.begin() ; i != T.end() ; i++ )
            {
                //Output the number of deleted location updates in this execution cycle -
                tup_deleted += i->second->evictPoint( (double)t - (double)global_range );
                //-----------------------------------------------------------------------

                //Partial results for this trajectory ----------------
                unordered_map<unsigned int, queryRes*> partialResults;
                partialResults.clear(); //----------------------------
                //----------------------------------------------------
                
                //Scan trajectory i and update the instant velocity aggregates
                for ( auto j = i->second->trajList.rbegin() ; j != i->second->trajList.rend() ; j++ )
                {
                    //Ignore the filthy location updates
                    if ( std::isnan((*j)->v) )
                    {
                        notAnumber++;
                        continue;
                    } //--------------------------------

                    //Include only non-zero qids to the result --------------------------
                    if ( (*j)->qid != 0 ) //---------------------------------------------
                    {
                        //First time encoutering this qid while traversing the trajectory
                        if ( partialResults.find((*j)->qid) == partialResults.end() )
                        {
                            //Increase the #queries that this trajectory contributes
                            invIndex[i->first]++;

                            queryRes *partial = new queryRes( (*j)->qid );
                            //cout << std::fixed << setprecision(3) << "oid:\t" << (partial->qid == (*j)->qid) << " sum:\t" << partial->sum << " count:\t" << partial->count << endl; 
                            partial->sum += (*j)->v;
                            partial->count++;

                            //Insert the new partial result into the map
                            partialResults.insert( {(*j)->qid, partial} );
                            //cout << "qid:\t" << partial->qid /*(partial->qid == (*j)->qid)*/ << " sum:\t" << partial->sum << " count:\t" << partial->count << endl; 
                        }
                        else
                        {
                            //Update the existing partial result
                            partialResults[(*j)->qid]->sum += (*j)->v;
                            partialResults[(*j)->qid]->count++;
                        }
                    }
                    // else
                    //     cout << "proto stigma -- oid:\t" << i->first << " qid:\t" << (*j)->qid << endl;
                //----------------------------------------------------------------------------------------------  
                }

                //The structure partialResults now contain all the contribution of this trajectory -- Insert them into global results
                for ( auto j = partialResults.begin() ; j != partialResults.end() ; j++ )
                {
                    if ( (Q.find(j->first) != Q.end()) && (j->second->count >= 1) ) //Already initialized final results for this execution cycle
                    {
                        //Here we insert the average velocities for each trajectory in each query
                        //-------------------
                        if      ( queryType == 1 )
                        {
                            results[j->first]->sum += j->second->sum / (double)(j->second->count);
                        }
                        else if ( queryType == 2 )
                        {
                            double new_partial_speed = j->second->sum / (double)(j->second->count);
                            results[j->first]->max   =  max(results[j->first]->max, new_partial_speed);
                        }
                        //-------------------------
                        results[j->first]->count++;
                        results[j->first]->oids.insert(i->first); //Include this oid to the set of oids contributing to this qid
                    }
                    // else cout << "COUNTER == 0" << endl;
                }

                //Memory Clearance ------------------------------------------------------
                for ( auto k = partialResults.begin() ; k != partialResults.end() ; k++ )
                    delete k->second; //-------------------------------------------------
                partialResults.clear();
                //-----------------------------------------------------------------------
            }
// cout << "3" << endl;

            //Now output the final velocity results by traversing the results map
            for ( auto i = results.begin() ; i != results.end() ; i++ )
            {
                if      ( op == 1 )
                    {;} //Do nothing with speed measurements
                else if ( op == 2 )
                {  
                    // Condition of processing queries
                    bool conditionara = true;
                    if (isGlobalWin == 1)
                        conditionara = (t >= (float)(Q[i->second->qid]->t + Q[i->second->qid]->range)) && (((t - (unsigned int)Q[i->second->qid]->t) % Q[i->second->qid]->slide) == 0);
                    else // Case with global window
                        conditionara = (t >= (unsigned int)(Q[i->second->qid]->t + global_range)) && ((((int)t - (int)Q[i->second->qid]->t) % global_slide) == 0);
                    // Condition readily calculated
                    if (conditionara)
                    { 
                        if (i->second->count == 0)
                        {     
                            cout << t << "," << i->first << ",-1," << i->second->count << endl;
                        }
                        else  
                        {
                            //-------------------
                            if (queryType == 1)
                            { 
                                if ( i->second->sum > 0.001 ) //If sum <> 0
                                    cout << t << "," << i->first << "," << i->second->sum / (double) i->second->count << "," << i->second->count << endl;
                                else
                                    cout << t << "," << i->first << ",-1," << i->second->count << endl;
                            }
                            else if (queryType == 2)
                            {
                                cout << t << "," << i->first << "," << i->second->max << "," << i->second->count << endl;
                            }
                            //-------------------
                        }
                    } // End of the condition for outputing a query result
                }
            }
// cout << "4" << endl;

            //Save the number of queries answered
            numQans = results.size();
            // //-----------------------------------
            // //Memory Clearance ----------------------------------------
            // for ( auto i = results.begin() ; i != results.end() ; i++ )
            //     delete i->second;
            // results.clear(); //----------------------------------------
            // //---------------------------------------------------------

        } //-------------------------------------------------------------------------------------------------------------------------------------

        //Calculate the number of profitable trajectories in this execution cycle
        profitable = 0;
        for ( auto i = invIndex.begin(); i != invIndex.end(); i++ )
            if ( invIndex[i->first] > 0 )
                profitable++;
        //-----------------------------------------------------------------------           

        //-------------------------------------------------
        evEndMIL = get_wall_time(); evEnd = get_cpu_time(); 
        //-------------------------------------------------

        //Query Evaluation Time ------
        cur_ev_load = evEnd - evStart; cur_ev_loadMIL = evEndMIL - evStartMIL; //Evaluation Time of this cycle in CPU cycles & MILliseconds
        evSum      += cur_ev_load;         //Total Evaluation Time in CPU cycles
        evSumMIL   += cur_ev_loadMIL;      //Total Evaluation Time in MILliseconds
        evCount++;                         //Increase evaluation counter
        //------------------------------------------------------------------------

// cout << "Evaluation Phase OK -- " << t << endl;        

        //-----------------------------------------------------------------------------------------------------------------------------//
        // Load Adaption phase --------------------------------------------------------------------------------------------------------//
        //-----------------------------------------------------------------------------------------------------------------------------//
        lsStart = get_cpu_time(); lsStartMIL = get_wall_time();
        //-----------------------------------------------------
        L    = cur_up_load + cur_ev_load + cur_ls_load;    //Calculate the total load L for t_now
        //------------------------------------------
        maxLoad = max(maxLoad, L); //Always keep the maximum observed Load in maxLoad

        retValue = 0;
        happenedLS = false;
        //----------------------
        if ( (schemaLS != 0) ) //Condition
        {
            //Naive Linear Regression ---------------------------------------------------------------
            trendData.insert( {t, (L - C)} ); //-----------------------------------------------------

            //If trendData.size() > global_range ==> erase
            if ( t > ( (numLoads -1) * global_slide) )
            {   //cout << "numLoads:\t" << numLoads; //cout << " before:\t" << trendData.size();
                while ( (trendData.size() != 0) && (trendData.begin()->first < (const int)(t - ( (numLoads - 1) * global_slide) )) )
                    trendData.erase( trendData.begin() );  //cout << " after:\t" << trendData.size() << endl;
            }

            tempNumOfLoads4pred = trendData.size();

            //----------------------------------------------------------
            for ( auto i = trendData.begin(); i != trendData.end(); i++)
            {
                x_reg.push_back( i->first  );
                y_reg.push_back( i->second );
            } //--------------------------------------------------------

            //Linear Regression Prediction
            double L_linreg  = giveEstimate(linearReg(x_reg, y_reg), t + period) + C;

            //Clean the registers
            x_reg.clear();
            y_reg.clear();

            //-------------------------------------------
            if (prediction == 1)
                L_pred = L_linreg;      // Find the expected load in the next execution cycle by linear intepolating on the last few noted loads
            else if (prediction == 2)
                L_pred = L;             // Assume that we expect (more or less) the same amount of load for next execution cycle
            else 
            {
                cout << "Wrong prediction INPUT" << endl;
                return 1;
            } //-----------------------------------------

            //Load Adaption Model
            if (fabs(L_pred - C) >= theta)         
            {   
                if (L_pred > C + theta) //Overload
                {
                    happenedLS = true;
                    
                    //Activate ----------
                    if ( !activeUnderLS )
                        activeUnderLS = true; //Activate UnderLoad LS every-
                    
                    //Fist time Load Adaption ------------------------------
                    if ( firstTimeLS )
                    {
                        n = T.size();         // Here the first time -- newbies are included ==> it is old + new trajectories
                        firstTimeLS = false;  //First time running LS boolean
                    }

                    traj_used = n;
                    //delta = round( ((C_destination - L_pred) / (double)L) * n ); // Here we must include the newbies as well
 /*sigspatial16*/   //delta = round( ((C_destination - L_pred) / (double)L) * (n /*+ newbies*/) ); // Here we must include the newbies as well
        /*mdm18*/   //delta = round(((L_pred - L)/(double)L) * n); 
        /*mdm18*/   //delta = round(((L_pred - L)/(double)max(L,L_pred)) * n);

                    // Function like walfare paper -- vldbj16
        /*mdm18*/   double z = 0.0;
                    double ratio = (abs(L_pred-C)*100.0)/C; // (L_pred - C)/C
                    if (ratio > 1.0)
                        z = ratio;
                    else
                        z = 1.0;
                    // Calculate the δ
                    delta = round( ((log2(z+1) * abs(L_pred-C)) * n) / (double)z ); //(L_pred - C)
                    // Subtract !!!
                    n -= delta;
            
//    cout << "delta_over : " << delta << endl;
    
                    cycles_ls++;
                    cycles_over++;

                    // for (unsigned int i = 1; i < (obj_cnt + 10); i++)
                    //     U[i] = 0;

                    loadS = new LoadShedder(grid, &T, &U, &invIndex, &results, &Q, &n, obj_cnt, 0, t, global_range);
                    if      ( schemaLS == 1 ) // SelectFirst
                        retValue = loadS->SelectFirstUnorderedLS(); 
                    else if ( schemaLS == 2 )
                        retValue = loadS->SelectFirstOrderedLS(); 
                    else if ( schemaLS == 3 )
                        retValue = loadS->PermutationSampleLS();
                    else if ( schemaLS == 4 ) // Reservoir Sample
                        retValue = loadS->ReservoirSampleLS_over();
                    else if ( schemaLS == 5 )
                        retValue = loadS->SemanticLS(loadS->FewObjectsFunc(firstCritPolicy));
                    else if ( schemaLS == 6 )
                    {
                        ;//retValue = loadS->RandomTupleLS_over();
                        //cout << "n = " << n << " * tuplesSoFar = " << tup_sofar + tup_total << " maintained = " << maintained << endl;

                    } //Dont do nothing for now! //retValue = loadS->RandomTupleLS();
                    // ---
                }
                else if ( (L_pred < C - theta) && activeUnderLS ) //Underload
                {
                    if (!firstTimeLS)
                    {
                        happenedLS = true;

                        traj_used = n;
                        //delta = round( ((C_destination - L_pred) / (double)L) * n );
    /*sigspatial16*/    //delta = round( ((C_destination - L_pred) / (double)L) * (n /*+ newbies*/) ); // Here we must include the newbies as well
          /*mdm18*/     delta = round( abs((L_pred - C)/(double)C) * n);
          /*mdm18*/     //delta = round(((L_pred - L)/(double)max(L,L_pred)) * n);

                        // Function like walfare paper -- vldbj16
          // /*mdm18*/     double z = 0.0;
          //               double ratio = (abs(L_pred-C)*100.0)/C; // (L_pred - C)/C
          //               if (ratio > 1.0)
          //                   z = ratio;
          //               else
          //                   z = 1.0;
          //               // Calculate the δ
          //               delta = round( ((log2(z+1) * abs(L_pred-C)) * n) / (double)z ); //(L_pred - C)

                        // Add !!!
                        n += delta; //CAUTION: signed delta needed

//    cout << "delta_under : " << delta << endl;

                        // settled
                        // if ((unsigned long)n >= T.size())
                        // {
                        //     n = T.size();
                        //     retValue = -1;
                        // }

                        cycles_ls++;
                        cycles_under++;

                        // for (unsigned int i = 1; i < (obj_cnt + 10); i++)
                        //     U[i] = 0;

                        loadS = new LoadShedder(grid, &T, &U, &invIndex, &results, &Q, &n, obj_cnt, 1, t, global_range);
                        if      ( schemaLS == 1 ) // SelectFirst
                            retValue = loadS->SelectFirstUnorderedLS(); 
                        else if ( schemaLS == 2 )
                            retValue = loadS->SelectFirstOrderedLS(); 
                        else if ( schemaLS == 3 )
                            retValue = loadS->PermutationSampleLS();
                        else if ( schemaLS == 4 ) // Reservoir Sample
                            retValue = loadS->ReservoirSampleLS_under();
                        else if ( schemaLS == 5 )
                            retValue = loadS->SemanticLS(loadS->FewObjectsFunc(firstCritPolicy));
                        else if ( schemaLS == 6 )
                        {
                            ;//retValue = loadS->RandomTupleLS_under();
                            //cout << "n = " << n << " * tuplesSoFar = " << tup_sofar + tup_total << " maintained = " << maintained << endl;

                        } // Do not do nothing for now! //retValue = loadS->RandomTupleLS();
                        // ---
                    }
                }

                //Check if the under LS should be deactivated
                if ( activeUnderLS )
                    if ( !firstTimeLS )
                        //If Load is low enough DEACTIVATE UnderLoadLS
                        if ( (L < C - theta) && (L_pred < C - theta) )
                            if ( retValue == -1 )
                            {
                                    activeUnderLS = false;         
                                    traj_used = n;   //Deactivate the UNDERload shedding
                            }
            }
            //-----------------------------------------------------------------------------------------------------
        }

        //-------------------------------------------------
        lsEndMIL = get_wall_time(); lsEnd = get_cpu_time(); 
        //-------------------------------------------------

        if ( schemaLS != 0 )
        {
            cur_ls_load    = lsEnd    - lsStart;
            cur_ls_loadMIL = lsEndMIL - lsStartMIL;
        }
        else
        {
            cur_ls_load    = 0.0;
            cur_ls_loadMIL = 0.0;
        }
        //Load Adaption Time ---
        lsSum    += cur_ls_load;    //Load Adaption Time in CPU cycles
        lsSumMIL += cur_ls_loadMIL; //Load Adaption Time in MILliseconds
        lsCount++;                  //Increase load adaption counter

// cout << "Load Shedding Phase OK -- " << t << endl;        

        //Resultados -- Time, Loads
        if ( op == 1 )
        {
            cout << std::fixed << t << ";" << setprecision(3) << L << ";" << setprecision(3) << L_pred << ";" << 
            setprecision(3) << C << ";" << setprecision(3) << C - theta << ";" << setprecision(3) << C + theta << ";" <<
            //newby, profitable, active_traj, total_traj 
            active << ";" << newbies << ";" << profitable << ";" << T.size() << ";" <<
            //n: #traj used in this exe cycle, n': #traj TO BE USED in the next exe cycle
            traj_used << ";" << delta << ";" << n << ";" <<
            //#tuple_processed, #tuples total
            tup_processed << ";" << tup_total << ";" << tup_sofar << ";" << maintained << ";" << tup_deleted << ";" << 
            //#UNLS, #OVLS, #totalLS, #cyclesTotal 
            cycles_under << ";" << cycles_over << ";" << cycles_ls << ";" << cycles_total << ";" << happenedLS << ";" <<
            //#queries_answered, case_of_queries 
            /*results.size() */ numQans << ";" <<
            //Times in CPU cycles
            setprecision(3) << cur_up_load << ";" << setprecision(3) <<  cur_ev_load << ";" << setprecision(3) << cur_ls_load << ";" << 
            //Times in MILliseconds
            setprecision(3) << cur_up_loadMIL << ";" << setprecision(3) <<  cur_ev_loadMIL << ";" << setprecision(3) << cur_ls_loadMIL << ";" << tempNumOfLoads4pred << endl;
   
        }

        inTuples.clear(); //Clear the tuples you read in the previous step

        //Maintain the total number of tuples so far in the experiment
        tup_sofar += tup_total;

        t = t + period; //Here δ = global_slide

    } //End of while loop - input file exhausted

    if ( op == 1 )
    {
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

        cout << "PIP correctness:\t" << (correctPIP / (double) totalPIP) * 100 << " %" << endl;
        cout << "#notAnum:\t" << notAnumber << endl;
        cout << "MAX Load:\t" << maxLoad << endl;
    }

    //Memory Deallocation -------------------------
    delete scanStream;
    //    delete grid;
    for ( auto i = T.begin() ; i != T.end() ; i++ )
        delete i->second;
    T.clear();
    U.clear();
    Q.clear();
    //lastQID.clear();
    preLastTuple.clear();
    lastTuple.clear();
    invIndex.clear();
    //----------------------------------------------

    return 0;
}