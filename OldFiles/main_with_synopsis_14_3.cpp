//Title: main.cpp 
//Description: Main class of the Load Shedding Trajectory System
//Author: Serafeim Papadias
//Platform(s): gcc (Ubuntu/Linaro 4.6.3-1ubuntu5) 4.6.3
//Date: 4/4/2015
//Revision: 22/2/2017

#include "Trajectory.cpp"
#include "RegularGrid.cpp"
#include "Synopsis.cpp"
#include "Scan.cpp"
#include "EvalStruct.cpp"
#include "parsingFunctions.cpp"
#include "winCPU.cpp"
#include "LoadShedder.cpp"
#include "linearRegression.cpp"
#include "numQueriesAtCycle_t.cpp"

using namespace std;

//Global Variables
/*unordered_*/map< unsigned int, Trajectory* > T;  
unordered_map< unsigned int, unsigned int > U;            
unordered_map< unsigned int, unsigned int > preU;         
unordered_map< unsigned int, map<unsigned int , EvalStruct*> > K;
unordered_set< unsigned int > Sk;
unordered_map< unsigned int, unsigned int > lastQID;
unordered_map< unsigned int, MyTuple* > lastTupleOfTraj;
map< unsigned int , TQuery* > queryMap;

//Keep track of the number of queries/roads each trajectory currently touches
map< unsigned int, unsigned int > numOfQueriesEachObjectIn;

double tmpValue = 0.0;
EvalStruct *evStruct;
list<point*>::iterator tempIt;
unsigned int tInitQry, slideQry, rangeQry;

//Evaluation Phase
double sum;
int counter;

unsigned int numOfCycles = 0, numOfCyclesLS = 0, numOfCyclesUnderLS = 0, numOfCyclesOverLS = 0, numOfProfitableTrajs = 0;
unsigned int totalNumTuplesPerExeCycle = 0, numTuplesProcessedPerExeCycle = 0, totalNumTuplesSoFar = 0, totalNumTrajSoFar = 0, numOfActiveTraj = 0;
int tempNumOfTrajKeptInThisCycle = 0, numOfTrajectoriesLastlyKept = 0, delta = 0;
unsigned int numOfQueriesToRespond = 0;
unsigned int numOfTrajFirst = 0, numOfTrajSecond = 0; //numOfTrajectoriesFirst + numOfTrajectoriesSecond = numOfTrajectoriesLastlyKept

//LoadShedder *loadS;                                                       //A pointer to a LoadShedder instance
double L = 0.0, L_expected = 0.0, L_naive = 0.0;                            //The actual load of our system
double L_eval1 = 0.0, L_eval2 = 0.0, L_eval3 = 0.0, L_eval4 = 0.0;          //The partial loads for PIP, K update and Query Evaluation
double Lnaive_n_1_DIF_L_actual_n   = 0.0; 
double tempLnaiveDIFF = 0.0;
bool belowCplusTheta = true;       
int retValue = 0;  //The return int value of the semantic LS function
bool firstTimeRunningLoadShedding = true;                                   //Bool for the first time load shedding

//Define
map<int, double> trendData;
map<int, double>::iterator itTrend;
vector<double> x_reg, y_reg;                                  

//Time Measurement in CPU cycles
double upStart, upEnd, evStart, evEnd, lsStart, lsEnd, boostStart1, boostEnd1, boostStart2, boostEnd2, boostStart3, boostEnd3, ls_inEv_start, ls_inEv_end;
double upSum = 0.0, evSum = 0.0, lsSum = 0.0, boostSum1 = 0.0, boostSum2 = 0.0, boostSum3 = 0.0;
unsigned int upCount = 1, evCount = 1, lsCount = 1;
double cur_up_load = 0.0, cur_ev_load = 0.0, cur_ls_load = 0.0;             //Loads for update, evaluation and load shedding for each timestamp
//Measure Time For maintaining & updating synopsis
double synoTime_start, synoTime_end;
double synoTime_sum = 0.0;

int main(int argc, char *argv[])
{
    if (argc != 18)  
    {
        //EXAMPLE:  ./synopsis VALID obj3.txt 3 0 10 5 0 10 10 
        cout << "Usage: " << argv[0] <<     "[mode]"                              << 
                                            "[input-file]"                        << 
                                            "{[stream-rate] | [timestamp-attr]}" << 
                                            "[start-time]"                        <<
                                            "[query-file]"                        <<
                                            "[Granularity_X]"                     <<
                                            "[Granularity_Y]"                     <<
                                            "[Capacity(C)_code]"                  <<
                                            "[applyFirstCriterion]"               <<
                                            "[firstCritPolicy]"                   << 
                                            "[WeightDecaying]"                    <<
                                            "[obj_count_expected]"                <<
                                            "[LS_scheme]"                         <<
                                            "[operation {tim | acc}]"             << 
                                            "[range]"                             <<
                                            "[slide]"                             << 
                                            "[randomQ-file]"                      << endl;

        exit(0);
    }

    cout.setf(ios::boolalpha);                                              //Print "true" and "false" instead of 1,0
    char *modeTime;                
    modeTime = argv[1]; 

    char *fileName;                
    fileName = argv[2];                                                     //ASCII file with object tuples of schema < t id x y qid >    

    bool mode = (strcmp(modeTime, "VALID") == 0);                           //True for VALID timestamps
    //mode = 0 for VALID timestamps , else mode = 1 for TRANS timestamps
    //Scan operator over the incoming stream file (with object locations))
    Scan * scanStream = new Scan(fileName, mode);
    //Depending on the mode of reading, the fourth argument specifies ...
    if (mode)
        scanStream->setTimeAttribute(atoi(argv[3]));                        //...the timestamp attribute in the schema of tuples
    else
        scanStream->setRate(atoi(argv[3]));                                 //...the stream arrival rate (in tuples/sec)
    //Specifies the time when the window is being applied
    int t0 = atoi(argv[4]);

    char *QueryFileName;                
    QueryFileName = argv[5];                //ASCII file with query tuples of schema < timestamp query_id qPolygon win_range win_slide >  

    //Granularity at X-axis
    unsigned int gx = atoi(argv[6]);                                        //Number of grid cells along x-dimension
    //Granularity at Y-axis
    unsigned int gy = atoi(argv[7]);                                        //Number of grid cells along y-dimension 

    scanStream->curTime = t0;                                               //Initialization: timestamp value when the window is being applied

    /*const*/ unsigned int C_code = atoi(argv[8]);

    //Boolean variable to keep if the first criterion will be applied or not
    unsigned int applyFirst = atoi(argv[9]);

    //Read the unsigned int which defines the first criterion policy for specifing the threshold
    unsigned int firstCritPolicy = atoi(argv[10]);
                            
    //If 0 -- No weighting scheme  applied , If 1 -- Weighting scheme applied to scores
    unsigned int boolWeightedScore = atoi(argv[11]);

    //Total number of moving objects that we expect to serve
    unsigned int obj_cnt = atoi(argv[12]);

    //Capacity in terms of the number of trajectories to process
    const unsigned int C = (C_code / 100.0) * obj_cnt; //e.g. 50% * 100000 = 50000 trajectories at most!

    //LS schema to use
    unsigned int schemaLS = atoi(argv[13]);

    //Type of operation
    //1 -- time experiment, 2 -- accuracy experiment
    unsigned int op = atoi(argv[14]);

    //Global Sliding Window settings
    int global_range = atoi(argv[15]);
    int global_slide = atoi(argv[16]);

    //Read only the queries in center of Athens
    char *queriesCenterFile;
    queriesCenterFile = argv[17];

    //Read the amount of incoming tuples per execution cycle -- from NO_LS experiment
    fstream finTup;
    finTup.open( queriesCenterFile, ios::in );
    //Hash table which keeps the number of incoming tuples for each execution cycle
    unordered_set<unsigned int> randomQ;
    unsigned int temp_q_id;
    //cout << "Size:\t" << randomQ.size() << endl;
    
    while ( finTup >> temp_q_id )
        randomQ.insert( temp_q_id );
    finTup.close();

    //Create grid partitioning according to the given parameters
    //box parameterBox = box(450000, 4175000, 500000, 4225000);       //Parameter box for Athens' synthetic dataset
    box parameterBox = box(1700000, 4500000, 1900000, 4700000);     //Parameter box for Rome's real dataset
    RegularGrid *grid = new RegularGrid(gx, gy, parameterBox);
    if (!grid->Allocate()) 
    {
        cout << "Memory Allocation Problem -- Grid partitioning!" << endl; 
        return 1;
    }

    //********************* Query Reading , Query Hashing to the grid , define global_range and MinWinSlide  *********
    fstream fin;
    fin.open( QueryFileName, ios::in );                                     //Open query file and read all of it
    unsigned int t_input, qid, range, slide, cid;
    string geom;
    list<unsigned int> ListOfCells;
    string X_Ystart, X_Yend;                                                //String for which has the point of the start and end vectors of each query
    double x1, y1, x2, y2;
    vector<string> orientationVerticesStart, orientationVerticesEnd;        //Define two tmp_vertices for start and end of the orientation vector

    //Scan the query file
    //while(fin >> t_input >> qid >> range >> slide >> X_Ystart >> X_Yend)    //Trick to avoid reading last input line twice!!!
    while(fin >> t_input >> qid >> cid >> X_Ystart >> X_Yend)               //Trick to avoid reading last input line twice!!!
    {
        //Check if inside the randomQ set of queries to receive answer
        if ( randomQ.find(qid) == randomQ.end() )
        {
            std::getline(fin, geom, '\n'); //Read the rest of the line
            continue;                      //and continue to the next iteration
        }

        //Create a new query tuple to store it to queryMap
        TQuery *query = new TQuery();
        //query->id = qid; query->ts = t_input; query->winRange = range; query->winSlide = slide;
        //Modification - Only one global window for all queries
        query->id = qid; query->ts = t_input; query->winRange = global_range; query->winSlide = global_slide; query->cid = cid;

        //Define the v1 = (x1,y1) and v2 = (x2,y2) vertices for calculating the orientation vector of each query
        orientationVerticesStart.clear();
        orientationVerticesStart = split(X_Ystart, SEPARATOR);

        orientationVerticesEnd.clear();
        orientationVerticesEnd = split(X_Yend, SEPARATOR);

        //Initialize min/max values according to coordinates of the first vertex
        x1 = atof(orientationVerticesStart[0].c_str()); y1 = atof(orientationVerticesStart[1].c_str());
        x2 =   atof(orientationVerticesEnd[0].c_str());   y2 = atof(orientationVerticesEnd[1].c_str());

        //Define the pointer to vertices we just created
        vertex *v1 = new vertex(x1, y1); vertex *v2 = new vertex(x2, y2);

        //Assign the values to query tuple that we created
        query->startV = v1; query->endV = v2;
        //Up until here we have read the points of the query orientation vector and we must read the query polygon 

        //Then, all remaining attribute values are considered as a VARYING set of coordinates defining a geometric shape (polygon) 
        std::getline(fin, geom, '\n');                                      //Copy remainder of this line as the geometry string

        std::pair<polygon*, polygon*> res = constructPolygon(geom);         //Create a polygon and its bounding box from the geometry string        
        polygon *p = res.first;                                             //The polygon object
        polygon *MBB = res.second;                                          //The MINIMUM bounding box - MBB

        query->qArea = p;                                                   //Assign the polygon object we read to the query area

        //Create an instance of a box structure (x_min,y_min,x_max,y_max) 
        box *MBBox = new box( MBB->vertices[0]->x , MBB->vertices[0]->y , MBB->vertices[1]->x , MBB->vertices[1]->y );

        //Hash the MBBox to see at which grid cell it falls into
        ListOfCells = grid->HashBox( *MBBox );
        for ( auto ListIt = ListOfCells.begin() ; ListIt != ListOfCells.end() ; ListIt++)
        {
            //Here we have (qID, vector<TQuery*>) format
            auto qryIter = grid->cell[(*ListIt)].QryInfo.find( query->id );      //Search for an entry of this query into the grid's QryInfo structure
            if (qryIter != grid->cell[(*ListIt)].QryInfo.end() )
            {
                //There is an entry already, thus we have already hashed another segment of this road
                grid->cell[(*ListIt)].QryInfo[query->id].push_back(query);  //Push the new pointer to the other part of this query at the end of the vector of query pointers
            }   
            else
            {
                //First time that we find an entry/a part of the road with this qID
                vector<TQuery*> tmp;                                       
                //tmp.clear();
                tmp.push_back(query);
                grid->cell[(*ListIt)].QryInfo.insert(pair< unsigned int, vector<TQuery*> >( query->id , tmp )); //Make and entry for this query at the hash structure of the grid
            }     
        }

        if ( queryMap.find(query->id) == queryMap.end() )
            queryMap.insert( pair<unsigned int , TQuery*>( query->id , query ) );

        if ( K.find(query->id) == K.end() )
            K[ query->id ];         
    }
    fin.close();
    //************** END OF Query Reading , Query Hashing to the grid , define global_range and MinWinSlide *********

    if ( queryMap.size() == 0 )
    {
        cout << "Kanena Erotima" << endl;
        exit(1);
    }

    vector<MyTuple*> inTuples;
    Trajectory *traj;                                                  //Pointer to a trajectory    

    //Initialize the booLS map structure to ones with respect to the expected amount of moving objects
    for (unsigned int i = 1; i < (obj_cnt + 1); i++ )
    {
        U.insert( pair< unsigned int, unsigned int >(i, 1) );
        lastTupleOfTraj[ i ];
    }
    //Now U has 1 to all the elements our system expects to serve
    preU = U;  

    //Resultados
    //Time, Loads
    if ( op == 1 )
    {
        cout << "t" << ";" << "L_actual" << ";" << "C" << ";" << 
        //#trajs from each criterion, total_#trajs
        "#1rst_Criterion" << ";" << "#2nd_Criterion" << ";" << "total_#traj" << ";" <<
        //newby, profitable 
        "#Newby_Updates" << ";" << "#traj_Profitable_This_ExeCycle" << ";" << "#traj_active_StillSendingTuples" << ";" <<
        //#tuple_processed, #tuples total
        "#tuples_Processed" << ";" << "#tuples_Total" << ";" <<
        //#UNLS, #OVLS, #totalLS, #cyclesTotal 
        "#performed_LS" << ";" << "LS: si o no?" << ";" << "#cycles_Total" << ";" <<
        //#queries_answered, case_of_queries 
        "#queries_Answered" << ";" <<
        //Times in CPU cycles
        "Update_time" << ";" << "SynUpdate_time" << ";" << "Eval_1--PIP" << ";" << "Eval_2--EvStrUPdate" << ";" << "Eval_3--Qans" << ";" << "Eval_4--VelCalc" << ";" << "EvalTotal_time" <<  ";" << "LS_time" << endl;
    }

    //Incremental Update of the TF_scores
    list< map < unsigned int, long double >* > slideTables;

    //Global clock -- t_NOW
    int t = 0;

    unsigned int numOfNewbyUpdates = 0;

    LoadShedder *loadS = new LoadShedder( &T, &K, &numOfTrajectoriesLastlyKept, &slideTables, queryMap.size(), boolWeightedScore );

    //STREAM INPUT - Keep processing data file until it gets exhausted
    while (scanStream->exhausted == false)
    {
        //Read streaming data -- Batch until timestamp t
        inTuples = scanStream->consumeInput(t);          

        //**********************************************************************************************
        //UPDATE PHASE: Create new tuples for the current timestamp value ******************************
        //**********************************************************************************************

        upStart = get_cpu_time(); 

        //New slide -- New cycle -- New slide table
        map< unsigned int, long double > *slide_t = new map< unsigned int, long double >;
        //Initialize slide_t table -- O(M)
        for ( auto qIt = queryMap.begin() ; qIt != queryMap.end() ; qIt++ )
            slide_t->insert( pair<unsigned int, long double>( qIt->first, 0.0 ) );
        //slideTables.push_back( slide_t );

        numOfNewbyUpdates = 0;
        numOfCycles++;

        totalNumTuplesPerExeCycle     = 0;
        numTuplesProcessedPerExeCycle = 0;

        boostSum1 = 0.0; boostSum2 = 0.0; boostSum3 = 0.0; synoTime_sum = 0.0;

        for ( auto it = inTuples.begin(); it != inTuples.end(); it++ )
        {   
            //Boolean to handle if the tuple coresponds to query or non-query road
            //bool perform = ( (*it)->loc->qid != 0 ) && ( randomQ.find( (*it)->loc->qid ) != randomQ.end() ); //Only for queries in the center
            //bool perform = ( (*it)->loc->qid != 0 ) && ( randomQ.find( (*it)->loc->qid ) == randomQ.end() ); //Only for queries NOT in the center
            //If perfrom = TRUE then EvalStruct will be performed, else NOT

            totalNumTuplesPerExeCycle++;

            //Initialize the doubles for time tracking to zeros
            boostStart1 = 0.0; boostEnd1 = 0.0; boostStart2 = 0.0; boostEnd2 = 0.0; boostStart3 = 0.0; boostEnd3 = 0.0;
            
            auto itT = T.find( (*it)->oid );
            if (itT == T.end())
            {
                totalNumTrajSoFar++;

                numOfNewbyUpdates++;
                numOfQueriesEachObjectIn[ (*it)->oid ] = 0;
                numTuplesProcessedPerExeCycle++;

                traj = new Trajectory((*it)->oid, t, grid);

                boostStart1 = get_cpu_time();
                traj->findQid((*it)->loc, 1);
                boostEnd1   = get_cpu_time();

                boostStart2 = get_cpu_time();
                traj->findVelocity((*it)->loc, 1);
                boostEnd2   = get_cpu_time();

                //Insert the new incoming tuple into the trajectory
            //    unsigned int qid_before = (*it)->loc->qid; 
                synoTime_start = get_cpu_time();
                double sTableTime = traj->insertPoint((*it)->loc, &U, t, global_slide, &lastQID, 0, &randomQ);
                synoTime_end   = get_cpu_time();
            //    unsigned int qid_after  = (*it)->loc->qid;
            //    if ( qid_before != qid_after )
            //    cout << "oid(" << (*it)->oid << ") -- qid before:\t" << qid_before << "\tqid after:\t" << qid_after << endl;
                synoTime_sum += synoTime_end - synoTime_start;
                if ( (*it)->loc->qid != 0 )
                    (*slide_t)[ (*it)->loc->qid ] += (long double)sTableTime;

                //Always insert the first tuple of every incoming object
                T.insert(pair< unsigned int, Trajectory* >((*it)->oid, traj));

                boostStart3 = get_cpu_time();
                //If this tuple belongs to an existing query...
                if ( (*it)->loc->qid != 0 ) 
                {
                    evStruct = new EvalStruct( (*it)->loc ,  T[(*it)->oid]->seq->segPoints.begin() );       
                    K[ (*it)->loc->qid ][ (*it)->oid ] = evStruct ; //Insert the first entry to the structure
                }
                boostEnd3   = get_cpu_time();
                
                //Important to know the lastQID of the object
                lastQID[ (*it)->oid ] = (*it)->loc->qid;
                //Important to maintain the last tuple of all objects
                lastTupleOfTraj[ (*it)->oid ] = (*it);

            }
            else 
            {
                //( (*it)->loc->qid == 0 ) || ( randomQ.find( (*it)->loc->qid ) == randomQ.end() ); //Only for queries in the center

                boostStart1 = get_cpu_time();

                //Adjust the qid of non-query tuples with non-zero qid to 0
                if ( ((*it)->loc->qid == 0) || (randomQ.find( (*it)->loc->qid ) == randomQ.end()) )
                    (*it)->loc->qid = 0;
                
                //Process only active trajectories & tuples referring to "changing to another road"
                if ( (U[ (*it)->oid ] == 1) || ( (*it)->loc->qid != lastQID[ (*it)->oid ] ) ) //CAUTION!!!!! CHECK IT OUTTTTT!!
                {
                    numTuplesProcessedPerExeCycle++;
                    //Measure time only for PIP
                    itT->second->findQid((*it)->loc, 2); //We must pay for PIP
                }
                boostEnd1 = get_cpu_time();

                boostStart2 = get_cpu_time();
                if ( U[ (*it)->oid ] == 1 )
                    if ( (*it)->loc->qid != 0 ) //If the qid belongs to a query
                        itT->second->findVelocitySpecial( (*it)->loc, lastTupleOfTraj[ (*it)->oid ]->loc );  
                boostEnd2 = get_cpu_time();

                //MALAKAAAAAAA AUTO PAEI PRIN TO EVALSTRUCT UPDATE!!!!!!!
                synoTime_start = get_cpu_time();
                double sTableTime = itT->second->insertPoint((*it)->loc, &U, t, global_slide, &lastQID, 1, &randomQ);
                synoTime_end   = get_cpu_time();
                synoTime_sum += synoTime_end - synoTime_start;
                if ( (*it)->loc->qid != 0 )
                    (*slide_t)[ (*it)->loc->qid ] += (long double)sTableTime; //Update the corresponding slides table time

                boostStart3 = get_cpu_time();
                if ( U[ (*it)->oid ] == 1 )
                    if ( (*it)->loc->qid != 0 )
                    {
                        //Query Processor - Traj already exists
                        auto ItEval = K[ (*it)->loc->qid ].find( (*it)->oid );          

                        if ( ItEval == K[ (*it)->loc->qid ].end() ) //What if many entries??? for the same query??? 
                        {
                            //Did not found an entry - so create one entry
                            tempIt = itT->second->seq->segPoints.end();
                            --tempIt; //Move the iterator(pointer) to the last element of the list
                            evStruct = new EvalStruct( (*it)->loc ,  tempIt );
                            K[ (*it)->loc->qid ][ (*it)->oid ] = evStruct; //Insert the first entry to the structure
                        }
                        else //The object continues to move at the same query - so update the front pointer of the entry (qid,oid) 
                        {
                            tempIt = itT->second->seq->segPoints.end(); //CAUTION: The point must have been inserted!!!!!!!
                            --tempIt; //Move the iterator(pointer) to the last element of the list
                            ItEval->second->update( (*it)->loc , tempIt );
                        }
                    }
                boostEnd3 = get_cpu_time();

                //Important to know the lastQID of the object
                lastQID[ (*it)->oid ] = (*it)->loc->qid;
                //Important to maintain the last tuple of all objects
                lastTupleOfTraj[ (*it)->oid ] = (*it); 
            }

            //Boost1 & Boost2 times for this execution cycle
            boostSum1 += boostEnd1 - boostStart1; //PIP time
            boostSum2 += boostEnd2 - boostStart2; //Velocity Calculation time
            boostSum3 += boostEnd3 - boostStart3; //EvalStruct Update time
        }

        //Insert the current slideTable in the list
        slideTables.push_back( slide_t );
        //cout << "slideTableList size(before):\t" << slideTables.size();
        while ( slideTables.size() > (((double)global_range) / ((double)global_slide)) )
            slideTables.pop_front();
        //cout << "\tslideTableList size(after):\t" << slideTables.size() << endl;

        //Initialize the number of active trajectories
        numOfActiveTraj = 0;
        //Update all the sypnoses of the existing trajectories
        synoTime_start = get_cpu_time();
        for ( auto itT = T.begin() ; itT != T.end() ; itT++ )
        {
            // if (t > 1700)
            // {
            //     if ( U[ itT->second->oid ] == 1 ) cout << "\nU = 1 -- active\n"; else cout << "U = 0 -- inactive\n****\n";  
            //     itT->second->reportState();
            //     //itT->second->synopsis->printSynopsis();
            // }
            //Crop the synopsis
            itT->second->synopsis->cropSynopsis( ((double)t - (double)global_range) );

            if ( itT->second->lastUpdateTime > (t - (double) (2 * global_slide)) )
                numOfActiveTraj++;
        }
        synoTime_end   = get_cpu_time();
        synoTime_sum += synoTime_end - synoTime_start;
        
        upEnd = get_cpu_time(); 

        //Times per Execution Cycle
        cur_up_load = (upEnd - upStart) - boostSum1 - boostSum2 - boostSum3;

        if ( numOfActiveTraj > C )
        {  
            //Total update time
            upSum += cur_up_load;                     //Update Time of this cycle in CPU cycles
            upCount++;                                //Increase update counter
        }

        //***************************************************************************************************************
        //EVALUATION PHASE **********************************************************************************************
        //***************************************************************************************************************

        evStart = get_cpu_time();

        //Initialize the number of profitable trajectories per execution cycle
        numOfQueriesEachObjectIn.clear();
        for ( auto itT = T.begin(); itT != T.end(); itT++ )
            numOfQueriesEachObjectIn[ itT->first ] = 0;
        
        //Temp pointer to the start and the end of each query
        list<point*>::iterator tmp_front,tmp_rear, tmp_preFront;

        Sk.clear(); //Clear the set of distinct object ids in K

        //Starting the evaluation with the help of the K structure
        //*************************************************************
        numOfQueriesToRespond = 0;                               
        for ( auto ItK = K.begin() ; ItK != K.end() ; ItK++ )
        {
            //Take query specifications -- (range, slide, tInit)
            if ( (randomQ.find(ItK->first) == randomQ.end()) || (queryMap.find(ItK->first) == queryMap.end()) ) //CAUTION: Here we also check the constraint that the query that we are answering belongs to the set
            {
                //cout << "Den brisko auto:\t" << ItK->first << endl;
                //cout << "PROBLIMA STIN ANAGNOSI EROTIMATON" << endl;
                continue;
            }
            else
            {
                tInitQry = queryMap.find(ItK->first)->second->ts;          
                slideQry = queryMap.find(ItK->first)->second->winSlide;    
                rangeQry = queryMap.find(ItK->first)->second->winRange;    
            }
            
            //Initialize the total sum and counter for each query
            sum     = 0.0;
            counter =   0;

            //Check query answering condition
            if ( (t >= tInitQry + rangeQry) && ((t - tInitQry) % slideQry == 0) )
            {
                //In this execution cycle increase the number of queries to respond
                numOfQueriesToRespond++;

                if ( ItK->second.size() > 0 )
                {
                    //Query answering condition 
                    for ( auto ItEval = ItK->second.begin() ; ItEval != ItK->second.end() ; )
                    {
                        //unsigned int temp_counter_before = ItEval->second->count;

                        tmp_front = ItEval->second->allParts.front()->front; //Oldest entry
                        tmp_rear  = ItEval->second->allParts.back()->rear;   //Newest entry

                        //Check if map entry is expired and must be deleted
                        if ( (*tmp_rear)->t < (double)(t - rangeQry) ) 
                        {
                            //******
                            //Insert deletion code here
                            //******

                            //Delete entries on the go IF needed
                            ItK->second.erase( ItEval++ );        //Note the post increment! //CAUTION:is that right?? check it!
                            //continue;                           //Proceed to the object ID for this query
                        }
                        else
                        {
                            //cout << "edo";
                            //Update EvalStruct before responding to the queries
                            while ( ( (*tmp_front)->t < (double)(t - rangeQry)  ) && ( distance(tmp_front,tmp_rear) > 1 ) )    //while front.time < NOW - WinRange
                            {
                                tmp_front = ItEval->second->refineFront();
                            }
                            //CAUTION: Maybe reconsider!!!!!!!

                            //Include this object's contribution
                            sum += ItEval->second->evaluate();
                            counter += 1; 

                            // unsigned int temp_counter_after = ItEval->second->count;
                            // if ( temp_counter_after == 0 )
                            //     cout << "counter_before:\t" << temp_counter_before << "\tcounter_after:\t" << temp_counter_after << endl;
                            
                            /*tmpValue = ItEval->second->evaluate(); 
                            if ( tmpValue > 0.0 )
                            {
                                sum     += tmpValue;
                                counter += 1; 
                            }*/

                            //Entry for TFlazy code here...

                            //Increase the entry of numOfQueriesEachObjecIn for this object
                            numOfQueriesEachObjectIn[ ItEval->first ]++;  //Increase by 1
                            Sk.insert( ItEval->first );

                            ++ItEval;
                            //cout << "\tekei\t" << endl;
                        }
                    }
                        
                    if ( op == 1 )
                    {
                        ;//Do nothing with speed measurements
                    }
                    else if ( op == 2 )
                        if ( counter == 0 )
                            //Do 
                            cout << t << "," << ItK->first << ",-1000," << ItK->second.size() << endl;
                            //cout << t << "," << ItK->first << ",0," << ItK->second.size() << endl; 
                            //OK. Return zero to compare with the ground truth!
                        else  
                            cout << t << "," << ItK->first << "," << (sum / (double) counter) << "," << ItK->second.size() << endl;

                    //Entry for TFlazy code here...
                }
                else
                {
                    if ( op == 1 )
                        {
                            ;//Do nothing with speed measurements
                        } 
                    else if ( op == 2 )
                        //timestamp,qid,avg_speed=0,0
                        cout << t << "," << ItK->first << ",-1000," << ItK->second.size() << endl;
                }
            }
        }

        //Initialize the number of profitable trajectories per execution cycle
        numOfProfitableTrajs = 0;
        for ( auto IterNumOfQ = numOfQueriesEachObjectIn.begin(); IterNumOfQ != numOfQueriesEachObjectIn.end(); IterNumOfQ++ )
        {
            if ( numOfQueriesEachObjectIn[ IterNumOfQ->first ] > 0 )
                //If the traj with oid belong to at least 1 query then its 
                numOfProfitableTrajs++; 
        }

        //Average evaluation time measurement
        evEnd = get_cpu_time();

        //Evaluation Time for this execution cycle
        cur_ev_load = (evEnd - evStart) + boostSum1 + boostSum2 + boostSum3;

        if ( numOfActiveTraj > C )
        {  
            evSum += cur_ev_load;         //Total Sum of Evaluation Time in CPU cycles
            evCount++;                    //Increase evaluation counter
        }

        //************************************************************************************************
        //LOAD SHEDDING PHASE ****************************************************************************
        //************************************************************************************************

        //lsStart = get_cpu_time(); 

        //Calculate the load L of this timestamp - t
        //cur_up_load -- from current cycle
        //cur_ev_load -- from current cycle
        //cur_ls_load -- from previous cycle
        L = cur_up_load + cur_ev_load + cur_ls_load;
        //Define the loads
        L_eval1 = boostSum1;                                  //PIP          -- f(Input Stream Rate)
        L_eval2 = boostSum3;                                  //K update     -- f(Input Stream Rate)
        L_eval3 = evEnd - evStart;                            //Query Eval   -- f(Query WorkLoad) 
        L_eval4 = boostSum2;                                  //Velocity calc-- constant      

        //Keep as many trajectories as the capacity specifies
        numOfTrajectoriesLastlyKept = C;

        bool happenedLS = false;
        //NO LS condition!
        if ( schemaLS != 0 )
        {
            if ( numOfActiveTraj > C )         
            {   
                happenedLS = true;
                preU = U;                                 //Copy U's contents

                //Increase the number of cycles that LS segment runs
                numOfCyclesLS++;

                //Set U entries to 0
                // for ( auto itT = T.begin() ; itT != T.end() ; itT++ )
                //     U[ itT->first ] = 0; 
            
                //Initialize the booLS map structure to ones with respect to the expected amount of moving objects
                for (unsigned int i = 1; i < (obj_cnt + 1); i++ )
                    U[ i ] = 0;

                //Does this line here increases the time?? Reasonable since it creates again and again a new loadS class!! CHECK it out!!!!!!!
                //loadS = new LoadShedder( &T, &K, &numOfTrajectoriesLastlyKept, &slideTables, queryMap.size(), boolWeightedScore );
                lsStart = get_cpu_time(); 
                    if      ( schemaLS == 1 )
                    {
                        // if ( firstTimeRunningLoadShedding )
                        // {
                            retValue = loadS->RandomLS( &U, 0 ); 
                        //    firstTimeRunningLoadShedding = false;
                        // }
                        // else
                        // {
                        //     ;//cout << "DEN KANO ALLO RANDOM" << endl;
                        // }
                    }          
                    else if ( schemaLS == 2 )                          
                        retValue = loadS->SemanticLS( loadS->FewObjectsFunc(firstCritPolicy), &U, &numOfTrajFirst, &numOfTrajSecond , applyFirst, 0 );
                    else if ( schemaLS == 3 )
                        retValue = loadS->SemanticBiasedLS( loadS->FewObjectsFunc(firstCritPolicy), &U, &numOfTrajFirst, &numOfTrajSecond , &Sk, applyFirst, 0 );
                    else if ( schemaLS == 4 )
                        retValue = loadS->TFLS( &U, 0 );
                    else if ( schemaLS == 5 )
                        retValue = loadS->TFbiasedLS( &U, &Sk, 0 );
               lsEnd = get_cpu_time();
            }
        }  
 
        //lsEnd = get_cpu_time();

        if ( schemaLS != 0 ) 
            cur_ls_load = lsEnd - lsStart;
        else
            cur_ls_load = 0.0;

        if ( numOfActiveTraj > C )
        {  
            lsSum += cur_ls_load;
            lsCount++;
        }

        //Resultados -- Time, Loads
        if ( op == 1 )
        {
            cout << t << ";" << setprecision(6) << L << ";" << setprecision(6) << C << ";" <<
            //#trajs from each criterion, total_#trajs
            numOfTrajFirst << ";" << numOfTrajSecond << ";" << T.size() << ";" <<
            //newby, profitable 
            numOfNewbyUpdates << ";" << numOfProfitableTrajs << ";" << numOfActiveTraj << ";" <<
            //#tuple_processed, #tuples total
            numTuplesProcessedPerExeCycle << ";" << totalNumTuplesPerExeCycle << ";" <<
            //#UNLS, #OVLS, #totalLS, #cyclesTotal 
            numOfCyclesLS << ";" << happenedLS << ";" << numOfCycles << ";" <<
            //#queries_answered, case_of_queries 
            setprecision(6) << numOfQueriesToRespond << ";" <<
            //Times in CPU cycles
            setprecision(6) << cur_up_load << ";" << setprecision(6) << synoTime_sum << ";" << setprecision(6) << L_eval1 << ";" << 
            setprecision(6) << L_eval2 << ";" << setprecision(6) << L_eval3 << ";" << setprecision(6) << L_eval4 << ";" <<
            setprecision(6) << cur_ev_load << ";" << setprecision(6) << cur_ls_load << endl;
        }     

        inTuples.clear(); //Clear the tuples you read in the previous step

        //Maintain the total number of tuples so far in the experiment
        totalNumTuplesSoFar += totalNumTuplesPerExeCycle;

        t = t + global_slide; //Here δ = global_slide
    } //End of while loop - input file exhausted
    
    cout << endl;
    cout << "**********************" << endl;    
    cout << scanStream->recCount << " records processed in total." << endl;

    delete scanStream;

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
    cout << ( (upCount == lsCount) && ( upCount == evCount ) ) << endl;

    return 0;
}