//Tile: mainGrid.cpp 
//Description: Entry point to management of regular grid partiioning over a point data stream.
//Author: Kostas Patroumpas
//Platform(s): gcc 4.8.3
//Date:     29/12/2007
//Revision: 11/05/2017 (previous 26/03/2017)

#include "./Utilities/structuresGrid.h"
#include "./RegularGrid/RegularGrid.cpp"
#include "./Utilities/structuresTraj.h"
#include "./Scan/Scan.cpp"
#include "./Utilities/winCPU.cpp"
#include "./Utilities/linearRegression.cpp"
#include "./LoadShedder/LoadShedder.cpp"

TObject *ObjData = NULL;

int main(int argc, char* argv[])
{			
	//---------------------------------------------------------------------------------------------------------------------------------//
	//----------- Parse the Arguments -------------------------------------------------------------------------------------------------//
	//---------------------------------------------------------------------------------------------------------------------------------//
	if (argc != 18) 
	{
		cerr << "\nUsage: " << argv[0] <<   "[mode]"                              <<
											"[{real | synthetic}]"                <<
											"[object-file]"                       <<
                                            "{[stream-rate] | [timestamp-attr]}"  << 
											"[start-time]"                        <<
											"[Granulariy_X]"                      <<
											"[Granulariy_Y]"                      <<
											"[Object_count]"                      <<
											"[Query__count]"                      <<
											"[Query_extent]";                          
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

    Scan *scanStream = new Scan(fileName, mode);
    //Depending on the mode of reading, the fourth argument specifies ...
    if (mode)
        scanStream->setTimeAttribute(atoi(argv[4]));                        
    else
        scanStream->setRate(atoi(argv[4]));                                
    
    //Specifies the time when the window is being applied
    int t0 = atoi(argv[5]);
    scanStream->curTime = t0;

    //Granulariy at X-axis
    unsigned int gx = atoi(argv[6]);
    //Granulariy at Y-axis
    unsigned int gy = atoi(argv[7]);

    //Total number of moving objects that we expect to come
    unsigned int obj_cnt = atoi(argv[8]);

    //Number of random queries to define
    unsigned int qry_cnt = atoi(argv[9]);

 `   //The extent of each query
    unsigned int extent  = atoi(argv[10]);

    //Read the input Lpeak observed at the experiment where NOLS took place
    const double Lpeak = atof(argv[11]);

    //Read the input Capaciy C code
    //10, 25, 33, 50, 67, 80% Lpeak
    unsigned int C_code = atoi(argv[12]);
    //Capaciy -- in terms of CPU cycles
    const double C = (double)((double)C_code * Lpeak)/100.0;
    // //Capaciy -- in terms of the #trajectories to process
    // const unsigned int C = (unsigned int)( (C_code / 100.0) * (double)obj_cnt ); //e.g. 50% * 100000 = 50000 trajectories at most!

    //Threshold for the secure band as a percentage of Lpeak
    double thresLS = atoi(argv[13]);

    //Calculate theta value from the saturation band
    double theta   = ((double)thresLS * Lpeak)/100.0;

    //Prediction cases
    //1 -- Use the simplistic linear regression on observed load in each execution cycle, 2 -- Other more sophisticated prediction model
    unsigned int prediction = atoi(argv[14]);

    //The time window [tnow - pred_window, tnow] in which we perform linear regression for Load prediction
    //60 -- 2 loads/point for regression, 120 -- 3 etc etc
    unsigned int pred_window = atoi(argv[15]); 

    // //Global Sliding Window settings
    // int global_range = atoi(argv[16]);
    // int global_slide = atoi(argv[17]);

    //LS schema to use
    unsigned int schemaLS = atoi(argv[16]);

    //Type of operation
    //1 -- time experiment, 2 -- accuracy experiment
    unsigned int op = atoi(argv[17]);

    //---------------------------------------------------------------------------------------------------------------------------------//
    //------------------ Global variables ---------------------------------------------------------------------------------------------//
    //---------------------------------------------------------------------------------------------------------------------------------//
    unsigned int delta = 10; 								//Period between successive execution cycles
    vector<MyTuple*> inTuples;  
    unsigned int t = 0; 									//Global Clock -- t_now
    unordered_map<unsigned int, TObject*> O;				//Map wih location updates for moving objects -- Only last tuple
	unordered_map<unsigned int, unsigned int> U; 
	unordered_map<unsigned int, TQuery*> queryMap; 			//Map wih queries
	unordered_map<unsigned int, MyTuple*> lastTuple;
	unordered_map<unsigned int, unsigned int> invIndex; 	//obj <==> number of queries i contributes to
	map<int, double> trendData; 							//CAUTION: Do not use unordered_map here
															//The elements are not stored sequentialy by key BUT we want them to be so
	//Statistics of the program
	unsigned int cycles_total = 0, cycles_ls = 0, cycles_over = 0, cycles_under = 0;
	unsigned int active = 0, newbies = 0, profiable = 0, traj_sofar = 0;
	unsigned int tup_total = 0, tup_processed = 0, tup_deleted = 0, tup_sofar = 0;
	unsigned int boolW = 0;
	
	//Load Shedding 
//	LoadShedder *loadS;
	double L = 0.0, L_pred = 0.0;	                         //The actual load of our system & the predicted one
	double dif_Lpred = 0.0, tempDIFFpred = 0.0;
	bool firstTimeLS = true, activeUnderLS = false;
	int traj_used = 0, n = 0, delta = 0;                     //#traj used in this execution cycle & n - number of traj to be used in the next
	vector<double> x_reg, y_reg;                                 

	//Time Measurement in CPU cycles
	double upStart = 0.0, upEnd = 0.0, evStart = 0.0, evEnd = 0.0, lsStart = 0.0, lsEnd = 0.0;
	double upSum = 0.0, evSum = 0.0, lsSum = 0.0;
	unsigned int upCount = 1, evCount = 1, lsCount = 1;
	double cur_up_load = 0.0, cur_ev_load = 0.0, cur_ls_load = 0.0;	//Loads for update, evaluation and load shedding for each timestamp

    //---------------------------------------------------------------------------------------------------------------------------------//
    //------------------ Create grid partiioning --------------------------------------------------------------------------------------//
    //---------------------------------------------------------------------------------------------------------------------------------//
	cout << "Iniializing..." << endl;

    box athensBox = box(450000 , 4175000,  500000, 4225000);     //Parameter box for Athens' synthetic dataset
    box romeBox   = box(1700000, 4500000, 1900000, 4700000);     //Parameter box for Rome's real dataset
    //Define the grid
    RegularGrid *grid;
    //Cases on dataset to be used

    if ( isRealDataset == 1 )
        grid = new RegularGrid( gx, gy, romeBox );
    else
        grid = new RegularGrid( gx, gy, athensBox );
    //Allocate the space
    if ( !grid->Allocate(obj_cnt) ) 
    {
        cout << "Memory Allocation Problem -- Grid partiioning!" << endl; 
        return 1;
    }

	//Create the random box queries & hash them to grid
	queryMap.clear();
	grid->randomQ(qry_cnt, extent, isRealDataset, &queryMap);
	
	//Check the number of queries
	cout << "#queries:\t" << queryMap.size() << endl;
    //Print the queries hashed on the grid
    grid->printGridState(1);

    //Iniialize the U bimap
    for ( unsigned int i = 1 ; i < (obj_cnt + 1) ; i++ )
	    U.insert( {i, 1} );

	//Time, Loads
    if ( op == 1 )
    {
        cout << "t" << ";" << "L_actual" << ";" << "L_pred" << ";" << "C" << ";" << "C-θ" << ";" << "C+θ" << ";" <<
        //newby, profiable, active_traj, total_traj 
        "#active" << ";" << "#newbies" << ";" << "#profiable" << ";" << "#total" << ";" << 
        //traj used, delta, trajs to be used
        "n" << ";" << "s" <<";" << "n'" << ";" <<
        //#tuple_processed, #tuples total
        "#tup_processed" << ";" << "#tup_total" << ";" << "#tup_sofar" << ";" << "#tup_deleted" << ";" << 
        //#UNLS, #OVLS, #totalLS, #cyclesTotal 
        "cycles_under" << ";" << "cycles_over" << ";" << "#cycles_ls" << ";" << "#cyclesTotal" << ";" <<  "LS:_si_o_no?" << ";" << 
        //#queries_answered, case_of_queries 
        "#queries_ans" << ";" <<
        //Times in CPU cycles
        "up_time" << ";" << "ev_time" <<  ";" << "ls_time" << endl;
    }	

    //STREAM INPUT - Keep processing data file until i gets exhausted
    while ( scanStream->exhausted == false )
    {
    	//Read streaming data -- Batch Reading until timestamp t
        inTuples = scanStream->consumeInput(t);          

        //----------------------------------------------------------------------------------------------//
        //UPDATE PHASE: Create new tuples for the current timestamp value ------------------------------//
        //----------------------------------------------------------------------------------------------//
     	cout << "Update Processing..." << endl;
        upStart = get_cpu_time();

        cycles_total++;
        newbies       = 0;
        tup_total     = 0;
        tup_processed = 0;
        tup_deleted   = 0;

        for ( auto i = inTuples.begin(); i != inTuples.end(); i++ )
        {
        	//Current Tuple
        	MyTuple *tup = (*i);

        	//Total number of location updates in this execution cycle
            tup_total++;    

            auto j = O.find( tup->oid );
            if ( j == O.end() )
            {
                traj_sofar++;
                newbies++;
                tup_processed++;

                //Create a new entry about this object -- assign velocity 0 for the first tuple
                TObject tmp_obj = new TObject(t, j->first, true, j->second->x, j->second->y, 0.0);
                //Insert the object's location update to map
                O.insert( {tup->oid, &tmp_obj} );
                //Create an entry in inverted index
                invIndex.insert( {tup->oid, 0} );
            }
            else
            {
                if ( U[ tup->oid ] == 1 )	//If it is an important one
                {
                    tup_processed++;

                    if ( findVelociy( tup->loc, lastTuple[tup->oid]->loc ) == -1 )
                    {
                        //Even if we do NOT include this location update into trajectory -- We update the timeLastUpdate
                        if ( O.find( tup->oid ) != O.end() )
                        	//Make the last location update non-fresh
                        	O[ tup->oid ]->fresh = false;
			            
			            //Maintaint the last tuple for each object
                        lastTuple[ tup->oid ] = tup;

                        continue;  
                    }

                    //Update the velocity of the object
                    j->second->t 	 = tup->loc->t;
                    j->second->x 	 = tup->loc->x;
                    j->second->y 	 = tup->loc->y;
                    j->second->v 	 = tup->loc->v;
                    j->second->fresh = true;
                }
            }

            //Maintaint the last tuple for each object
            lastTuple[ tup->oid ] = tup;
        }

		//Update the object location on grid
		for ( auto i = O.begin() ; i != O.end() ; i++ )
		{
			if ( i->second->fresh )				  //Only fresh objects should be processed!
			{
				active++;        
				grid->UpdateObject( *i->second );
				i->second->fresh = false;		  //Invalidate 'fresh' indicator for the next execution cycle
			}
		}

        upEnd = get_cpu_time();

        //Times per Execution Cycle
        cur_up_load = upEnd - upStart;

        upSum += cur_up_load;                     //Update Time of this cycle in CPU cycles
        upCount++;                                //Increase update counter

        //---------------------------------------------------------------------------------------------------------------//
        //EVALUATION PHASE ----------------------------------------------------------------------------------------------//
        //---------------------------------------------------------------------------------------------------------------//
     	cout << "Query Processing..." << endl;
    
        evStart = get_cpu_time();

        //Query Evaluation -- O(M)
        for ( auto i = queryMap.begin() ; i != queryMap.end() ; i++ )
        {
        	//This vector keeps the ids of the objects that reside inside the i->first query
        	vector<long> qryResults;
        	qryResults.clear();
        	grid->rangeQuery( i->second->qBox, qryResults);

        	//Revise the object to be included in the Query Results
        	for ( auto j = qryResults.begin() ; j != qryResults.end() ; j++ )
        	{
        		if ( U[ *j ] != 1 ) 				//If this object is UNimportant
        			qryResults.erase(j);			//erase it from qryResults
        	}

        	//Temp values
        	double 	 	   vel_sum = 0.0;
    		unsigned int vel_count = 0;

        	//Report Results -- Average instanteneous speed of objects inside query
        	for ( auto i = qryResults.begin() ; i != qryResults.end() ; i++ )
        	{
        		//Update the sum & count
				vel_sum += O[ *i ].v;
				vel_count++;         		
        	}

        	//Print out the results for this query
        	if      ( op == 1 )
                ; //Do nothing with speed measurements
            else if ( op == 2 )
                if ( qryResults.size() == 0 )
                    cout << t << "," << i->first << ",-1000," << qryResults.size() << endl;
                else  
                    if ( vel_sum > 0.00000000001 ) //If sum <> 0
                        cout << t << "," << i->first << "," << vel_sum / (double) vel_count << "," << qryResults.size() << endl;
                    else
                        cout << t << "," << i->first << ",-1000," << qryResults.size() << endl;

        	//Statistics
			cout << "QUERY contains " << qryResults.size() << "important objects:";
	        for ( auto i = qryResults.begin() ; i != qryResults.end() ; i++ )
	            cout << ' ' << *i;
	        cout << endl;

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

        evSum += cur_ev_load;         //Total Sum of Evaluation Time in CPU cycles
        evCount++;                    //Increase evaluation counter

        // //------------------------------------------------------------------------------------------------//
        // //LOAD SHEDDING PHASE ----------------------------------------------------------------------------//
        // //------------------------------------------------------------------------------------------------//
     	cout << "Load Shedding Processing..." << endl;

     	bool happenedLS = false;
        // lsStart = get_cpu_time(); 


        // //LS CODE SEGMENT HERE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        // lsEnd = get_cpu_time();

        // if ( schemaLS != 0 )
        //     cur_ls_load = lsEnd - lsStart;
        // else
        //     cur_ls_load = 0.0;

        // lsSum += cur_ls_load;
        // lsCount++;

        //Resultados -- Time, Loads
        if ( op == 1 )
        {
            cout << std::fixed << t << ";" << setprecision(3) << L << ";" << setprecision(3) << L_pred << ";" << 
            setprecision(3) << C << ";" << setprecision(3) << C - theta << ";" << setprecision(3) << C + theta << ";" <<
            //newby, profitable, active_traj, total_traj 
            active << ";" << newbies << ";" << profitable << ";" << O.size() << ";" <<
            //n: #traj used in this exe cycle, n': #traj TO BE USED in the next exe cycle
            traj_used << ";" << delta << ";" << n << ";" << 
            //#tuple_processed, #tuples total
            tup_processed << ";" << tup_total << ";" << tup_sofar << ";" << tup_deleted << ";" << 
            //#UNLS, #OVLS, #totalLS, #cyclesTotal 
            cycles_under << ";" << cycles_over << ";" << cycles_ls << ";" << cycles_total << ";" << happenedLS << ";" <<
            //#queries_answered, case_of_queries 
            queryMap.size() << ";" <<
            //Times in CPU cycles
            setprecision(3) << cur_up_load << ";" << setprecision(3) <<  cur_ev_load << ";" << setprecision(3) << cur_ls_load << endl;
        }

    	inTuples.clear(); 						  //Clear the tuples you read in the previous step

        tup_sofar += tup_total;					  //Maintain the total number of tuples so far in the experiment

        t = t + delta; 						      //Here δ = delta

    } //End of while loop - input file exhausted

	//-----------------------------------------------------------------------//
	//--------- Memory Deallocation -----------------------------------------//
	//-----------------------------------------------------------------------//
	delete scanStream;
	delete[] ObjData;
	ObjData = NULL;
	delete grid;
	queryMap.clear();
	
	cout << "*****************************************************" << endl;
	cout << "Execution completed successfully!" << endl;
	return 0;
}
