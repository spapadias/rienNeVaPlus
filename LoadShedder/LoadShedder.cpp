//Title: LoadShedder_priority_queue_scoring.cpp
//Description: A class of LoadShedder which applies the two criterions for choosing the significant trajectories
//Author: Serafeim Papadias
//Platform(s): gcc (Ubuntu/Linaro 4.6.3-1ubuntu5) 4.6.3
//Date:     04/04/2015
//Revision: 20/01/18 (previous --> 22/03/2017)

#include "LoadShedder.h"

LoadShedder::LoadShedder(
						 RegularGrid *_grid,
	 					 unordered_map< unsigned int, Trajectory* >  *_T, 
					 	 unordered_map< unsigned int, unsigned int > *_U, 
						 unordered_map< unsigned int, unsigned int > *_invIndex, 
					     unordered_map< unsigned int, queryRes*>     *_results,
	            	     unordered_map< unsigned int, TQuery*>       *_Q,
					     int *_n,
					     unsigned int _obj_cnt,
					     unsigned int _isUnder,
					     unsigned int _tnow,
					     unsigned int _range
						)
{
	grid 	 = _grid;
	T 		 = _T;
	U        = _U;
	invIndex = _invIndex;
	results  = _results;
	Q        = _Q;
	n   	 = _n;
	obj_cnt  = _obj_cnt;
	isUnder  = _isUnder;
	tnow     = _tnow;
	range  	 = _range;
}

LoadShedder::~LoadShedder()
{
}

// SRS for every incoming tuple. Now assume n = 3000 tids to maintain. So, roughly 3*3000=9000 tuples to maintain.
// Also, tupleSoFar is the number of tuples processed so far. Thus, for each incoming tuple we process it with 
// probability 3*n/tupleSoFar.
// int LoadShedder::RandomTupleLS_tuple(unsigned int tuplesArrivedUpTillNow, unsigned int maintained)
// {
// 	// Random process -- Uniformly Distributed Numbers
// 	int min = 0.0;
// 	int max = 1.0;
// 	// srand(unsigned(time(NULL)));
// 	// int output_1 = min + (rand() % static_cast<int>(max - min + 1));	
// 	std::random_device rd;     // only used once to initialise (seed) engine
// 	std::mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)
// 	std::uniform_real_distribution<> dis(min, max);
// 	auto output_2 = dis(rng);

// 	unsigned int random_number = output_2;
// 	double acceptance_probability = 3.0*(*n) / ((double) tuplesArrivedUpTillNow);

// 	if (random_number >= acceptance_probability) // With probability 3n/tupleSoFar --> Process the tuple
// 		return 0; 	 // --> process the tuple
// 	else
// 		return 666;  // --> Do not process the tuple
// }

int LoadShedder::SelectFirstUnorderedLS() 
{
	int counter = *n;

	if (counter < 0)
	{
//		cout << "FTASAME 0" << endl; //CAUTION: Problem if n reaches ZERO
		return 1; 
	}
	else if ((unsigned long)*n >= T->size())
	{
		counter = T->size();
		*n = T->size();

		if (isUnder == 1)
		{
			for (auto i = T->begin() ; i != T->end() ; i++)
				(*U)[ i->first ] = 1;
			return -1; 								  // All tids utilized 																						
		}
	}	

	// for (auto i = T->begin() ; i != T->end() ; i++)
	// 	(*U)[ i->first ] = 0; 		 				  // Invalidate all trajectory ids
	
	for (unsigned int i = 1; i < (obj_cnt + 10); i++) // Invalidate all tids
        (*U)[ i ] = 0;	

	for (auto i = T->begin() ; i != T->end() ; i++)   // Maintain the first n entries of T
	{
		(*U)[ i->first ] = 1;
		counter--;
		if (counter <= 0)
			break;
	}

	return 0; 										  // positive n < T.size
}

int LoadShedder::SelectFirstOrderedLS() 
{
	int counter = *n;

	if (counter < 0)
	{
//		cout << "FTASAME 0" << endl; //CAUTION: Problem if n reaches ZERO
		return 1; 
	}
	else if ((unsigned long)*n >= T->size())
	{
		counter = T->size();
		*n = T->size();

		if (isUnder == 1)
		{
			for (auto i = T->begin() ; i != T->end() ; i++)
				(*U)[ i->first ] = 1;
			return -1; 								  // All tids utilized 																						
		}
	}	

	set<unsigned int> ordered_set_T;

	for (auto i = T->begin() ; i != T->end() ; i++)
		ordered_set_T.insert(i->first);
	// (*U)[ i->first ] = 0; 		 				  // Invalidate all trajectory ids
	
	for (unsigned int i = 1; i < (obj_cnt + 10); i++) // Invalidate all tids
        (*U)[ i ] = 0;	

	for (auto i = ordered_set_T.begin() ; i != ordered_set_T.end() ; i++)   // Maintain the first n entries of T
	{
		(*U)[ *i ] = 1;
		counter--;
		if (counter <= 0)
			break;
	}

	return 0; 										  // positive n < T.size
}

int LoadShedder::PermutationSampleLS()
{
	int counter = *n;
	vector<unsigned int> ids;

	//-----------------
	if (counter <= 0)
	{
//      cout << "FTASAME O" << endl;
		return 1;		
	}
	else if ((unsigned long)counter >= T->size())
	{
		*n       = T->size();
		counter  = T->size(); 																		 

		if ( isUnder == 1 )
		{
			for ( auto i = T->begin() ; i != T->end() ; i++ )
				(*U)[i->first] = 1;
			return -1; 																						
		}
	} //-----------------------------------------------------

	// for ( auto i = T->begin() ; i != T->end() ; i++ )
	// 	(*U)[i->first] = 0;

	for (unsigned int i = 1; i < (obj_cnt + 10); i++) // Invalidate all tids
        (*U)[ i ] = 0;

	//Create a vector with all object/trajectory ids
	for ( auto i = T->begin() ; i != T->end() ; i++ )
		ids.push_back(i->first);

	// srand(unsigned(time(NULL)));
	// random_shuffle( ids.begin(), ids.end() );

	auto rng = std::default_random_engine {};
	shuffle(begin(ids), end(ids), rng);

	//Randomly select #ofTrajectoriesLastlyKept trajectories & update U
	for( auto i = ids.begin() ; i != ids.end() ; i++ )
	{
		(*U)[*i] = 1;
		//--------
		counter--;
		if ( counter <= 0 )
			break;
	}

	return 0;
}

int LoadShedder::ReservoirSampleLS_tuple(unsigned int tid_new, vector<unsigned int> *active_ids)
//CAUTION : Create also a procedure to define only the number n to keep each time --> in the same place you have the other methods --> online algorithm
{
// cout << "mesa1" << endl;
	if (T->size() < (unsigned long)*n) // Not more than n trajectory ids
		(*U)[ tid_new ] = 1;
	else 				// More than n trajectory ids
	{
		// Random process -- Biased Non-uniformly Distributed Numbers
		int min = 0;
		int max = T->size() - 1;
		// srand(unsigned(time(NULL)));
		// int output_1 = min + (rand() % static_cast<int>(max - min + 1));
// cout << "mesa4" << endl;
	
		//C++11 variant -- Uniformly Distributed Numbers
		std::random_device rd;     // only used once to initialise (seed) engine
		std::mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)
		std::uniform_int_distribution<int> uni(min,max); // guaranteed unbiased
		auto output_2 = uni(rng);
// cout << "mesa5" << endl;
		// Here maintain the qids randomly
    	// auto rng = std::default_random_engine {};
    	// shuffle(begin(active_ids), end(active_ids), rng);

		unsigned int random_number = output_2;
		if (random_number < active_ids->size()) //*n)
		{
			(*U)[ (*active_ids)[random_number] ] = 0;
			(*U)[ tid_new ] = 1; 
		}
// cout << "mesa6" << endl;		
	}

	return 0;
}

int LoadShedder::RandomTupleLS_over() 
{
	int counter = *n;
	
	if (counter <= 0)
	{
//      cout << "FTASAME O" << endl;
		return 1;		
	}
	else if ((unsigned long)counter >= T->size())
	{
		*n       = T->size();
		counter  = T->size(); 																		 

		if ( isUnder == 1 )
		{
			for ( auto i = T->begin() ; i != T->end() ; i++ )
				(*U)[i->first] = 1;
			return -1; 																						
		}
	} // ---

	return 0;
}

int LoadShedder::RandomTupleLS_under()
{
	int counter = *n;

	if (counter <= 0)
	{
//      cout << "FTASAME O" << endl;
		return 1;		
	}
	else if ((unsigned long)counter >= T->size())
	{
		*n       = T->size();
		counter  = T->size(); 																		 

		if ( isUnder == 1 )
		{
			for ( auto i = T->begin() ; i != T->end() ; i++ )
				(*U)[i->first] = 1;
			return -1; 																						
		}
	} // ---

	return 0;
}

int LoadShedder::ReservoirSampleLS_over()
{
    vector<unsigned int> temp_active_ids;
	int counter = *n;

	if (counter <= 0)
	{
//      cout << "FTASAME O" << endl;
		return 1;		
	}
	else if ((unsigned long)counter >= T->size())
	{
		*n       = T->size();
		counter  = T->size(); 																		 

		if ( isUnder == 1 )
		{
			for ( auto i = T->begin() ; i != T->end() ; i++ )
				(*U)[i->first] = 1;
			return -1; 																						
		}
	} // ---

    for (auto i = U->begin() ; i != U->end() ; i++)
        if (i->second == 1)
            temp_active_ids.push_back(i->first); // Gather the trajectory ids with U[id]=1
    
    if (temp_active_ids.size() > (unsigned long)n)
    {
        int temp_diff = temp_active_ids.size() - counter; //*n;
        auto rng = std::default_random_engine {};
        shuffle(begin(temp_active_ids), end(temp_active_ids), rng);

        for(auto i = temp_active_ids.begin() ; i != temp_active_ids.end() ; i++)
        {
            (*U)[*i] = 0; // Invalidate the exceeding active trajectories
            //--------
            temp_diff--;
            if (temp_diff <= 0)
                break;
        }
    }

    return 0;
}

int LoadShedder::ReservoirSampleLS_under()
{
    //   --- At each round we must ensure that the reservoir has no more than n active trajectories
    vector<unsigned int> temp_active_ids;
    vector<unsigned int> all_ids;
	int counter = *n;

	if (counter <= 0)
	{
//      cout << "FTASAME O" << endl;
		return 1;		
	}
	else if ((unsigned long)counter >= T->size())
	{
		*n       = T->size();
		counter  = T->size(); 																		 

		if ( isUnder == 1 )
		{
			for ( auto i = T->begin() ; i != T->end() ; i++ )
				(*U)[i->first] = 1;
			return -1; 																						
		}
	} // ---

    // for (auto i = T->begin() ; i != T->end() ; i++)
    // {
    //     all_ids.push_back(i->first);
    //     if ((*U)[i->first] == 1)
    //         temp_active_ids.push_back(i->first);
    // }
    
    // // If we have exceeding number of active trajectories based on the previous LS decision --> flush some : (i) randomly, (ii) semantically
    // if (temp_active_ids.size() < (unsigned long)n)
    // {
    //     vector<unsigned int> temp_unimportant_ids;
    //     // Find the set difference
    //     for (auto i = all_ids.begin() ; i != all_ids.end() ; i++)
    //     {
    //         auto j = find(temp_active_ids.begin(), temp_active_ids.end(), *i);
    //         if (j == temp_active_ids.end()) // If the id *i belongs to an UNimportant trajectory
    //             temp_unimportant_ids.push_back(*i);
    //     }

    //     int temp_diff = counter - temp_active_ids.size(); // *n instead of counter, also
    //     // Randomly permute the active ids
    //     auto rng = std::default_random_engine {};
    //     shuffle(begin(temp_unimportant_ids), end(temp_unimportant_ids), rng);

    //     //Randomly select #ofTrajectoriesLastlyKept trajectories & update U
    //     for( auto i = temp_unimportant_ids.begin() ; i != temp_unimportant_ids.end() ; i++ )
    //     {
    //         (*U)[*i] = 1; // Validate more trajectories randomly
    //         //--------
    //         temp_diff--;
    //         if ( temp_diff <= 0 )
    //             break;
    //     }
    // }

    return 0;
}

//Function that keeps the trajectories which belong to some queries and they form a small set in contrast with the other set in other queries
//Usage
//Parameter func: 1 - mean : If the crowd ( #ofTrajectories ) is below the mean crowd then keep all these trajectory belonging to this query
//Parameter func: 2 - median:If the crowd ( #ofTrajectories ) is below the median crowd then keep all these trajectory belonging to theis query
//Parameter func: 3 - minumum thres: If the crowd  ( #ofTrajectories ) is below the func threashold then keep all these trajectory belonging to this query
//Returns the mean or the median of the the crowd sets
unsigned int LoadShedder::FewObjectsFunc(unsigned int func)
{
	//func is mean
	if (func == 1)
	{
		double mean = 0.0;
		double sum = 0.0;
		unsigned int countOfSets = 1;

		for (auto i = results->begin() ; i != results->end() ; i++)
		{
			sum += i->second->count;
			countOfSets++;
		} 

		if ((countOfSets-1) != 0) 
			mean = sum / ((double)(countOfSets - 1));
		else
			mean = 0.0;

		return (unsigned int) mean;
	}
	else if (func >= 2)
	{
		//Define a priority queue to store the pairs (qID , crowdset's size)
		priority_queue< pair<unsigned int, unsigned int>, vector< pair< unsigned int, unsigned int> >, CompareDistQuantile > priorityQueueQuantile;
		for (auto i = results->begin() ; i != results->end() ; i++)
			//Push the pair (qID , sum_crowdset)
			priorityQueueQuantile.push(make_pair(i->first, i->second->count));

		unsigned int numOfQueriesRunning     = priorityQueueQuantile.size(); //The num of queries that are running at this timestamp
		unsigned int numOfElementsInEachTile =  (unsigned int)( numOfQueriesRunning / ((double) func) ); //Define the crowdset sum that we will use as a threshold
		unsigned int numOfElementsToDrop     = numOfElementsInEachTile; //numOfQueriesRunning - numOfElementsInEachTile;

		while (numOfElementsToDrop > 0)
		{
			priorityQueueQuantile.pop();
			numOfElementsToDrop--;
		}

		return (unsigned int) priorityQueueQuantile.top().second ;
	}
	else 
		return -1;
}

int LoadShedder::SemanticLS(unsigned int funcRes)
{
	unsigned int counter = *n;
	unordered_set<unsigned int> tids_firstc;

	if (counter <= 0)
	{
//      cout << "FTASAME O" << endl;
		return 1;		
	}
	else if ((unsigned long)counter >= T->size())
	{
		*n       = T->size();
		counter  = T->size(); 																		 

		if ( isUnder == 1 )
		{
			for ( auto i = T->begin() ; i != T->end() ; i++ )
				(*U)[i->first] = 1;
			return -1; 																						
		}

		return -1;
	} // ---

	// for (auto i = T->begin() ; i != T->end() ; i++)
	// 	(*U)[ i->first ] = 0;

	for (unsigned int i = 1; i < (obj_cnt + 10); i++) // Invalidate all tids
        (*U)[ i ] = 0;

	// //First Criterion
	// for (auto i = results->begin() ; i != results->end() ; i++)
	// {
	// 	if ((i->second->count <= funcRes) && (funcRes != 0))
	// 		for (auto j = i->second->oids.begin() ; j != i->second->oids.end() ; j++)
	// 			if (counter > 0)
	// 			{
	// 				(*U)[ *j ] = 1;
	// 				tids_firstc.insert(*j); // Mark this tid
	// 				counter--;
	// 			}
	// }

	// TF-scores
	unordered_map<unsigned int, double> tf_score;
	tf_score.clear();
	vector<unsigned int> leftovers, tf_selected;
	priority_queue< pair< unsigned int, double>, vector< pair< unsigned int, double > >, CompareDist> priorityScores;

	//Initialization of the scores -- O(N)
	for (auto i = T->begin() ; i != T->end() ; i++)
		if (tids_firstc.find(i->first) == tids_firstc.end()) // If we did not include this trajectory with 1rst criterion
			tf_score.insert( {i->first, 0.0} );

	//Based on synopsis derive the TF scores as if they were all contributed to queries CAUTION: Really slow!! -- O(M x N) 
	for (auto i = results->begin() ; i != results->end() ; i++)
	{	
		if (i->second->count > 0) // not necessary
		{
			//cout << "OIDS set size:\t" << (i->second->oids.size() == i->second->count) << endl;

			for (auto j = i->second->oids.begin() ; j != i->second->oids.end() ; j++)
				if (tids_firstc.find(*j) == tids_firstc.end()) // If this object should get ranked
					tf_score[ *j ] += 1.0 / (double) i->second->oids.size();
		}
	}

	//Find the tf_selected & leftovers sets -- O(N)
	for (auto i = tf_score.begin() ; i != tf_score.end() ; i++)
		if ( tf_score[ i->first ] > 0.0000000001 )
			tf_selected.push_back(i->first);
		else
			leftovers.push_back(i->first);

	//Distinguish cases -- Total Time Complexity O(N)
	if (tf_selected.size() < counter)
	{
		for (auto i = tf_selected.begin() ; i != tf_selected.end() ; i++)
		{
			(*U)[ *i ] = 1;
			counter--;
		}

		auto rng = std::default_random_engine {};
        shuffle(begin(leftovers), end(leftovers), rng);

		//Randomly pick the leftovers that received a zero score
		for(auto i = leftovers.begin() ; i != leftovers.end() ; i++)
		{
			(*U)[ *i ] = 1;
			//Check how many trajectories have we kept
			counter--;
			if ( counter <= 0 )
				break;
		}
	}
	else
	{
		for (auto i = tf_selected.begin() ; i != tf_selected.end() ; i++)
		 	priorityScores.push( make_pair( *i, tf_score[ *i ] ) );

		while (counter > 0)
		{
			(*U)[ priorityScores.top().first ] = 1; 
			priorityScores.pop();
			counter--;
		}
	}

	// // Second Criterion
	// priority_queue< pair<unsigned int, double>, vector< pair< unsigned int, double> >, CompareDist> priorityQueueScoring;
	// double tmp_score = 0.0;
	// for (auto i = T->begin() ; i != T->end() ; i++)
	// { 
	// 	tmp_score = i->second->scoringOfTrajectory(U, timeBreak);
	// 	if (abs( tmp_score - (-2) ) != 0)
	// 		priorityQueueScoring.push( make_pair( i->first , tmp_score ) ); 
	// } 

	// //Apply 2nd Criterion
	// while ( trajCounter > 0 )
	// {
	// 	(*U)[ priorityQueueScoring.top().first ] = 1;
	// 	priorityQueueScoring.pop();
	// 	trajCounter--;
	// 	(*numSecond)++;
	// }

	return 0;
}











// --- --- ---



int LoadShedder::StratifiedLS()
{
	int counter = *n;
	vector<unsigned int> ids;

	//-----------------
	if (counter <= 0)
	{
		counter = -66666666; 
		cout << "SKATAAAAA" << endl;
		//Problem with COUNTER = 0 ====> TRY it on big data set
		// for ( auto i = T->begin() ; i != T->end() ; i++ )
		// 		(*U)[i->first] = 0;	
		//return -2; //Keep zero location updates in the next cycle 			
	}	
	else if ((unsigned long)counter >= T->size())
	{
		*n       = T->size();
		counter  = T->size(); 																		 

		if ( isUnder == 1 )
		{

			for ( auto i = T->begin() ; i != T->end() ; i++ )
				(*U)[i->first] = 1;
			
			return -1; 																						
		}
	} //-----------------------------------------------------

	//Erase U --------------------------------------	
	for ( auto i = T->begin() ; i != T->end() ; i++ )
		(*U)[i->first] = 0;

	//Find the number of grid cells C' <= C cells
	unsigned int denominator = 0;
	set<unsigned int> effectiveCells;
	this->grid->objectsInQueryCells(&effectiveCells, &denominator);

	// cout << "-------" << endl;
	// cout << "Number of cells to sample n trajectories from:\t" << effectiveCells.size() << endl;
	// cout << "Number of objects to sample from:\t" <<  denominator << "\tand n = " << counter << endl;

	//Set random seed ----------
	srand(unsigned(time(NULL)));
	//--------------------------

	//For every cell sample proportionally and uniformly at random
	//...
	unsigned int summation = 0;
	for ( auto i = effectiveCells.begin() ; i != effectiveCells.end() ; i++ )
	{
		vector<unsigned int> ids;
		//---
		// unsigned int samplingRate = (unsigned int) round( (counter * this->grid->cell[*i].ObjList.size()) / (double) denominator );
		//unsigned int samplingRate = (unsigned int) trunc( (counter * this->grid->cell[*i].ObjList.size()) / (double) denominator );
		unsigned int samplingRate = (unsigned int) ceil( (counter * this->grid->cell[*i].ObjList.size()) / (double) denominator );

		summation += samplingRate;
		//---
		//--- Put every trajectory id in each effective cell in a vector/bucket
		for ( auto j = this->grid->cell[*i].ObjList.begin() ; j != this->grid->cell[*i].ObjList.end() ; j++ )
			ids.push_back(j->first);
		//Permute -------------------------------
		random_shuffle( ids.begin(), ids.end() );
		//---
		//Keep as many as the samplingRate ------
		for( auto k = ids.begin() ; k != ids.end() ; k++ )
		{
			(*U)[*k] = 1;
			//-----------
			samplingRate--;
			if ( samplingRate <= 1 ) //Keep one less from each cell
				break;
		}
		ids.clear();
		//---
	}
	//---

	//cout << "Trajectories selected:\t";
	unsigned int picked = 0;
	for ( auto i = T->begin() ; i != T->end() ; i++ )
		// if ( i->second == 1 )
		if ( (*U)[i->first] == 1 )
		{
			//cout << i->first << ", ";
			picked++;
		}
	//Plus + The newly incoming ones
	// cout << "#traj selected:\t" << picked << "\tcounter:\t" << counter << "\tsummation:\t" << summation  << endl;
	// cout << "-------" << endl;

	//---
	return 0;
}

int LoadShedder::StratifiedLS_2()
{
	int counter = *n;
	vector<unsigned int> ids;

	//-----------------
	if (counter <= 0)
	{
		counter = -66666666; 
		cout << "SKATAAAAA" << endl;
		//Problem with COUNTER = 0 ====> TRY it on big data set
		// for ( auto i = T->begin() ; i != T->end() ; i++ )
		// 		(*U)[i->first] = 0;	
		//return -2; //Keep zero location updates in the next cycle 			
	}	
	else if ((unsigned int)counter >= T->size())
	{
		*n       = T->size();
		counter  = T->size(); 																		 

		if ( isUnder == 1 )
		{

			for ( auto i = T->begin() ; i != T->end() ; i++ )
				(*U)[i->first] = 1;
			
			return -1;																					
		}
	} //-----------------------------------------------------

	//Erase U --------------------------------------	
	for ( auto i = T->begin() ; i != T->end() ; i++ )
		(*U)[i->first] = 0;

	//Find the number of grid cells C' <= C cells
	unsigned int denominator = 0;
	set<unsigned int> effectiveCells;
	this->grid->objectsInQueryCellsHistory(&effectiveCells, &denominator);

	// cout << "-------" << endl;
	// cout << "Number of cells to sample n trajectories from:\t" << effectiveCells.size() << endl;
	// cout << "Number of objects to sample from:\t" <<  denominator << "\tand n = " << counter << endl;	

	//Set random seed ----------
	srand(unsigned(time(NULL)));
	//--------------------------

	//For every cell sample proportionally and uniformly at random
	//...
	unsigned int summation = 0;
	for ( auto i = effectiveCells.begin() ; i != effectiveCells.end() ; i++ )
	{
		vector<unsigned int> ids;
		//---
		// unsigned int samplingRate = (unsigned int) round( (counter * this->grid->cell[*i].HisObjList.size()) / (double) denominator );
		//unsigned int samplingRate = (unsigned int) trunc( (counter * this->grid->cell[*i].HisObjList.size()) / (double) denominator );
		unsigned int samplingRate = (unsigned int) ceil( (counter * this->grid->cell[*i].HisObjList.size()) / (double) denominator );     //Proportionally BUT not correctly!

		summation += samplingRate;
		//---
		//--- Put every trajectory id in each effective cell in a vector/bucket
		for ( auto j = this->grid->cell[*i].HisObjList.begin() ; j != this->grid->cell[*i].HisObjList.end() ; j++ )
			ids.push_back(*j);
		//Permute -------------------------------
		random_shuffle( ids.begin(), ids.end() );
		//---
		//Keep as many as the samplingRate ------
		for( auto k = ids.begin() ; k != ids.end() ; k++ )
		{
			(*U)[*k] = 1;
			//-----------
			samplingRate--;
			if ( samplingRate <= 0 )
				break;
		}
		ids.clear();
		//---
	}
	//---

	//cout << "Trajectories selected:\t";
	unsigned int picked = 0;
	for ( auto i = T->begin() ; i != T->end() ; i++ )
		// if ( i->second == 1 )
		if ( (*U)[i->first] == 1 )
		{
			//cout << i->first << ", ";
			picked++;
		}
	//Plus + The newly incoming ones
//	cout << "#traj selected:\t" << picked << "\tcounter:\t" << counter << "\tsummation:\t" << summation << endl;
//	cout << "-------" << endl;

	//The HISTORICAL SAMPLING in general keeps LESS than expected ==> counter

	//--- The rest do stratifiedLS on NOW ==> COMBINATION OF NOW + HISTORY

	if ( picked < (unsigned int)counter )
	{
		cout << "ALSO stratified sampling for NOW" << endl;
		//Find the number of grid cells C' <= C cells
		denominator = 0;
		effectiveCells.clear();
		this->grid->objectsInQueryCells(&effectiveCells, &denominator);

		for ( auto i = effectiveCells.begin() ; i != effectiveCells.end() ; i++ )
		{
			vector<unsigned int> ids;
			//---
			// unsigned int samplingRate = (unsigned int) round( (counter * this->grid->cell[*i].HisObjList.size()) / (double) denominator );
			unsigned int samplingRate = (unsigned int) trunc( (counter * this->grid->cell[*i].HisObjList.size()) / (double) denominator );
			// unsigned int samplingRate = (unsigned int) ceil( ( (counter - picked) * this->grid->cell[*i].ObjList.size()) / (double) denominator );     //Proportionally BUT not correctly!

			//---
			//--- Put every trajectory id in each effective cell in a vector/bucket
			for ( auto j = this->grid->cell[*i].ObjList.begin() ; j != this->grid->cell[*i].ObjList.end() ; j++ )
				ids.push_back(j->first);
			//Permute -------------------------------
			random_shuffle( ids.begin(), ids.end() );
			//---
			//Keep as many as the samplingRate ------
			for( auto k = ids.begin() ; k != ids.end() ; k++ )
			{
				(*U)[*k] = 1;
				//-----------
				samplingRate--;
				if ( samplingRate <= 1 )
					break;
			}
			ids.clear();
			//---
		}


	}
	else
		cout << "ARKOUN" << endl;

	//cout << "Trajectories selected:\t";
	picked = 0;
	for ( auto i = T->begin() ; i != T->end() ; i++ )
		// if ( i->second == 1 )
		if ( (*U)[i->first] == 1 )
		{
			//cout << i->first << ", ";
			picked++;
		}
	//Plus + The newly incoming ones
		cout << "#traj selected:\t" << picked << "\tcounter:\t" << counter << "\tsummation:\t" << summation << endl;
		cout << "-------" << endl;

	//The HISTORICAL SAMPLING in general keeps LESS than expected ==> counter

	//---
	return 0;
}
















// int LoadShedder::TF_noHistory(unsigned int countPIP, double sumPIP)
// {
// 	int counter = *n; 		
// 	priority_queue< pair< unsigned int, double>, vector< pair< unsigned int, double > >, CompareDist> priorityScores;
// 	unordered_map<unsigned int, unsigned int> countQ; 		//< qid, oCount > -- NOW
// 	unordered_map<unsigned int, unsigned int> o2q;    		//< oid, qid >    -- NOW
// 	unordered_map<unsigned int, double> tf_score;     		//TF score        -- NOW
// 	unsigned int nonZero = 0;

// 	//Check the counters values -- for out of bounds cases
// 	//-------------------------------------------------------
// 	if ( counter <= 0 ) //-----------------------------------
// 		counter = 0; 		
// 	else if ( counter >= T->size() )
// 	{
// 		*n 		= T->size(); 
// 		counter = T->size(); 																		 

// 		if ( isUnder == 1 )
// 		{
// 			for ( auto i = T->begin() ; i != T->end() ; i++ )
// 				(*U)[ i->first ] = 1;
// 			return -1; 																						
// 		}
// 	} //-----------------------------------------------------
// 	//-------------------------------------------------------
	
// 	//Initialize countQ unordered map structure -- O(M) -----
// 	for ( auto i = Q->begin() ; i != Q->end() ; i++ )
// 		countQ.insert( {i->first, 0} );
// 	//-------------------------------------------------------

// 	//Assignment of Trajectory Scores ------------------------------ O(N) 
// 	for ( auto i = T->begin() ; i != T->end() ; i++ )
// 	{
// 		if ( i->second->mrut > (double) tnow - (double) range )
// 		{
// 			//Simulate PIP time by waiting the average time
// 			// if ( (*U)[ i->first ] != 1 )
// 			// 	usleep( sumPIP / (countPIP * 1000.0) ); // sumPIP/countPIP in MILliseconds and divided by 1000 ==> microseconds

// 			//Find the last tuple's qid
// 			unsigned int qid = (*i->second->last)[ i->first ]->qid_cheat;
				
// 			//Check if it belongs to a query
// 			if ( i->second->Q->find(qid) == i->second->Q->end() )
// 				qid = 0;

// 			//Built the partial structures
// 			o2q.insert( {i->first, qid} );

// 			if ( qid != 0 )
// 			{
// 				countQ[ qid ]++; //Increase the count
// 			}
// 		}
// 	}
// 	//-------------------------------------------------------------------

// 	//Finalize the TF scores --------------------------------------- O(N)
// 	for ( auto i = T->begin() ; i != T->end() ; i++ )
// 	{
// 		//Consider ONLY active trajectories ----------------------------------------------------------------------------------
// 		if ( i->second->mrut > (double) tnow - (double) range ) //------------------------------------------------------------
// 		{
// 			if ( o2q[ i->first ] == 0 )
// 				tf_score[ i->first ] = 0.0; 							//Zero score if movement in non-query
// 			else
// 			{
// 				tf_score[ i->first ] = 1.0 / (double)countQ[ o2q[ i->first ] ]; //Proportional score to number of objects in the query  
// 				nonZero++;
// 			}
// 		}
// 		else
// 			tf_score[ i->first ] = -10.0;								//Assign negative score
// 		//-------------------------------------------------------------
// 		//Finally, insert the tf score into the trajectory
// 	 	priorityScores.push( make_pair( i->first, tf_score[ i->first ] ) );
// 	}
// 	//------------------------------------------------------------------------------------------------------------------------

// //	cout << "nonZero:\t" << nonZero << "\tperc:\t" << (nonZero * 100.0) / (double)T->size() << endl;

// 	// Erase the U map	
// 	for ( auto i = T->begin() ; i != T->end() ; i++ )
// 		(*U)[ i->first ] = 0;

// 	while ( counter > 0 )
// 	{
// 		(*U)[ priorityScores.top().first ] = 1;
// 		//Print scores -----------------------------
// 		//cout << priorityScores.top().second << endl;
// 		//------------------------------------------
// 		priorityScores.pop();
// 		counter--;
// 	}
// 	//-----------------------------------------------------------------

// 	return 0;
// }

// int LoadShedder::RandomLS_FisherYates()
// {
// 	int counter = *n;		  //Num Trajs for the LS to maintain

// 	//Check the counters values -- for out of bounds cases --
// 	//-------------------------------------------------------
// 	if ( counter <= 0 )
// 		counter = 0; 		
// 	else if ( counter >= T->size() )
// 	{
// 		*n       = T->size();
// 		counter  = T->size(); 																		 

// 		if ( isUnder == 1 )
// 		{

// 			for ( auto i = T->begin() ; i != T->end() ; i++ )
// 				(*U)[ i->first ] = 1;
// 			return -1; 																						
// 		}
// 	} //------------------------------------------------------

// 	vector<unsigned int> elements;

// 	//Create a vector with all object/trajectory ids
// 	for ( auto i = T->begin() ; i != T->end() ; i++ )
// 		elements.push_back( i->first );

// 	srand(time(NULL));

// 	// the counter to be used to generate the random number
// 	int currentIndexCounter  = elements.size();
// 	for ( auto currentIndex  = elements.rbegin() ; currentIndex != elements.rend() - 1 ; currentIndex++, --currentIndexCounter )
//     {
//         int randomIndex = rand() % currentIndexCounter;
//         // if its not pointing to the same index      
//         if ( *currentIndex != elements.at(randomIndex) )
//         {
//             //then swap the elements
//             swap( elements.at(randomIndex), *currentIndex );
//         }
//     }

// 	// Erase the U map	
// 	for ( auto i = T->begin() ; i != T->end() ; i++ )
// 		(*U)[ i->first ] = 0;

// 	//Randomly select #ofTrajectoriesLastlyKept trajectories & update U
// 	for( auto i = elements.begin() ; i != elements.end() ; i++ )
// 	{
// 		(*U)[ *i ] = 1;
// 		//Check how many trajectories have we kept
// 		counter--;
// 		if ( counter <= 0 )
// 			break;
// 	}

// 	return 0;
// }

// int LoadShedder::QueriesPerObject()
// {
// 	int counter = *n;		  //Num Trajs for the LS to maintain
// 	priority_queue< pair< unsigned int, double>, vector< pair< unsigned int, double > >, CompareDist> priorityScores;
// 	unordered_map<unsigned int, double> score;        //Cell based synopsis score for a trajectory of a moving object

// 	//Check the counters values -- for out of bounds cases --------------
// 	//-------------------------------------------------------------------
// 	if ( counter <= 0 )
// 		counter = 0; 		
// 	else if ( counter >= T->size() )
// 	{
// 		*n       = T->size();
// 		counter  = T->size(); 																		 

// 		//Check if we are below or above the band -- If below then return
// 		if ( isUnder == 1 )
// 		{

// 			for ( auto i = T->begin() ; i != T->end() ; i++ )
// 				(*U)[ i->first ] = 1;
// 			return -1; 																						
// 		}
// 	}

// 	//------- Assign Scores Here --------------------
// 	for ( auto i = T->begin() ; i != T->end() ; i++ )
// 	{
// 		if ( i->second->mrut > (double) tnow - (double) range )
// 			//score[i->first] = i->second->synopsis->scoring();
// 			score[i->first] = (*invIndex)[i->first];
// 		else //In case of inactive object
// 			score[i->first] = -10.0;
// 		priorityScores.push( make_pair( i->first, score[i->first] ) );
// 	}

// 	// Erase the U map	
// 	for ( auto i = T->begin() ; i != T->end() ; i++ )
// 		(*U)[ i->first ] = 0;

// 	while ( counter > 0 )
// 	{
// 		(*U)[ priorityScores.top().first ] = 1;
// 		//Print scores -----------------------------
// //		cout << priorityScores.top().second << endl;
// 		//------------------------------------------
// 		priorityScores.pop();
// 		counter--;
// 	}

// 	return 0;
// }

// int LoadShedder::CombinationLS()
// {
// 	int counter = *n;							//Num of trajs for the LS to maintain
// 	priority_queue< pair< unsigned int, double>, vector< pair< unsigned int, double > >, CompareDist> priorityScores;
// 	unordered_map<unsigned int, double> score;  //Cell based synopsis score for a trajectory of a moving object

// 	//-----------------
// 	if ( counter <= 0 )
// 		counter = 0; 		
// 	else if ( counter >= T->size() )
// 	{
// 		*n       = T->size();
// 		counter  = T->size(); 																		 

// 		if ( isUnder == 1 )
// 		{
// 			for ( auto i = T->begin() ; i != T->end() ; i++ )
// 				(*U)[ i->first ] = 1;
// 			return -1; 																						
// 		}
// 	} //-----------------------------------------------------

// //cout << "counter:\t" << counter << endl;

// 	// Erase the U map ------------------------------
// 	for ( auto i = T->begin() ; i != T->end() ; i++ )
// 		(*U)[ i->first ] = 0;
// 	//-----------------------------------------------

// 	//Half of the trajectories semantically
// 	int half_1 = round(counter / 2.0); // cout << "half_1:\t" << half_1 << endl;
// 	///////////////////////////////////////
// 	set<unsigned int> semElements;

// 	//------- Assign Scores Here --------------------
// 	for ( auto i = T->begin() ; i != T->end() ; i++ )
// 	{
// 		if ( i->second->mrut > (double) tnow - (double) range )
// 			//score[i->first] = i->second->synopsis->scoring();
// 			score[ i->first ] = (*invIndex)[ i->first ];
// 		else //In case of inactive object
// 			score[ i->first ] = -10.0;
// 		priorityScores.push( make_pair( i->first , score[ i->first ] ) );
// 	}

// 	// for ( auto i = results->begin() ; i != results->end() ; i++ )
// 	// {	
// 	// 	if ( i->second->count > 0 )
// 	// 	{
// 	// 		for ( auto j = i->second->oids.begin() ; j != i->second->oids.end() ; j++ )
// 	// 			score[*j] = score[*j] * (1.0 / (double) i->second->oids.size());
// 	// 	}
// 	// }

// 	while ( half_1 > 0 )
// 	{
// 		(*U)[ priorityScores.top().first ] = 1;
// 		semElements.insert( priorityScores.top().first ); //Insert to the set of already semantically picked elements
// 		//Print scores -----------------------------
// 		cout << priorityScores.top().second << endl;
// 		//------------------------------------------
// 		priorityScores.pop();
// 		--half_1;
// 	}

// 	//Half of the trajectories randomly
// 	int half_2 = round(counter / 2.0); // cout << "half_2:\t" << half_2 << endl;
// 	///////////////////////////////////
// 	vector<unsigned int> elements;
// 	set<unsigned int> rndElements;

// 	// Insert the to be randomly choosen elements into elements vector
// 	for ( auto i = T->begin() ; i != T->end() ; i++ )
// 	{
// 		if ( (*U)[ i->first ] == 0 )
// 			elements.push_back( i->first );
// 	}

// //cout << "1" << endl;
// //cout << "elements size:\t" << elements.size() << endl;

// 	srand((time(0)));		//Initialize random generator

// 	// using built-in random generator:
//   	random_shuffle( elements.begin() , elements.end() );

// //cout << "2" << endl;

// 	// Randomly select half of the required elements
// 	for( auto i = elements.begin() ; i != elements.end() ; i++ )
// 	{
// 		(*U)[ *i ] = 1;
// 		rndElements.insert( *i ); //Insert to the set of already randomly picked elements
// 		//Check how many trajectories have we kept
// 		half_2--;
// 		if ( half_2 <= 0 )
// 			break;
// 	}

// //cout << "3" << endl;

// 	// //Test it!
// 	// cout << semElements.size() << " Semantically Selected:\t";
// 	// for ( auto i = semElements.begin() ; i != semElements.end() ; i++ )
// 	// 	cout << *i << "\t";
// 	// cout << endl;

// 	// cout << rndElements.size() << " Randomly Selected:\t";
// 	// for ( auto i = rndElements.begin() ; i != rndElements.end() ; i++ )
// 	// 	cout << *i << "\t";
// 	// cout << endl;

// 	return 0;
// }


// OLD BUT NOT FORGOTTEN



// int LoadShedder::TF()
// {
// 	int counter = *alpha; 		
// 	vector<unsigned int> leftovers, tf_selected;
// 	priority_queue< pair< unsigned int, double>, vector< pair< unsigned int, double > >, CompareDist> priorityScores;

// 	//Check the counters values -- for out of bounds cases
// 	//----------------------------------------------------
// 	if ( counter <= 0 )
// 		counter = 0; 		
// 	else if ( counter >= T->size() )
// 	{
// 		*alpha  = T->size(); 
// 		counter = T->size(); 																		 

// 		//Check if we are below or above the band -- If below then return
// 		if ( isUnder == 1 )
// 		{
// 			for ( auto i = T->begin() ; i != T->end() ; i++ )
// 				(*U)[ i->first ] = 1;
// 			return -1; 																						
// 		}
// 	}
// 	//----------------------------------------------------

// 	// tfs
// 	unordered_map<unsigned int, double> tf_score;
// 	tf_score.clear();

// 	//Initialization of the scores -- O(N)
// 	for ( auto i = T->begin() ; i != T->end() ; i++ )
// 		tf_score.insert( {i->first, 0.0} );

// 	//Based on synopsis derive the TF scores as if they were all contributed to queries CAUTION: Really slow!! -- O(M x N) 
// 	for ( auto i = results->begin() ; i != results->end() ; i++ )
// 	{	
// 		if ( i->second->count > 0 )
// 		{
// 			//cout << "OIDS set size:\t" << (i->second->oids.size() == i->second->count) << endl;

// 			for ( auto j = i->second->oids.begin() ; j != i->second->oids.end() ; j++ )
// 				tf_score[ *j ] += 1.0 / (double) i->second->oids.size();
// 		}
// 	}

// 	//Find the tf_selected & leftovers sets -- O(N)
// 	for ( auto i = T->begin() ; i != T->end() ; i++ )
// 		if ( tf_score[ i->first ] > 0.0000000001 )
// 			tf_selected.push_back( i->first );
// 		else
// 			leftovers.push_back( i->first );

// 	// //Statistics
// 	// if ( tf_selected.size() >= *alpha )
// 	// 	cout << "non_zero_scores:\t" << *alpha << " random_leftovers:\t" << 0 << endl;
// 	// else
// 	// 	cout << "non_zero_scores:\t" << tf_selected.size() << " random_leftovers:\t" << *alpha - tf_selected.size() << endl;

// 	//Distinguish cases -- Total Time Complexity O(N)
// 	if ( tf_selected.size() <= *alpha )
// 	{
// 		//-----------------------------------------------------------------
// 		for ( auto i = tf_selected.begin() ; i != tf_selected.end() ; i++ )
// 		{
// 			(*U)[ *i ] = 1;
// 			counter--;
// 		}

// 		//Perform a random permutation of the leftovers
// 		srand(unsigned(time(NULL)));
// 		random_shuffle( leftovers.begin(), leftovers.end() );

// 		//Randomly pick the leftovers that received a zero score
// 		for( auto i = leftovers.begin() ; i != leftovers.end() ; i++ )
// 		{
// 			(*U)[ *i ] = 1;
// 			//Check how many trajectories have we kept
// 			counter--;
// 			if ( counter <= 0 )
// 				break;
// 		}
// 		//------------------------------------------------------------------
// 	}
// 	else
// 	{
// 		//----------------------------------------------------------------------
// 		for ( auto i = tf_selected.begin() ; i != tf_selected.end() ; i++ )
// 		 	priorityScores.push( make_pair( *i, tf_score[ *i ] ) );

// 		while ( counter > 0 )
// 		{
// 			(*U)[ priorityScores.top().first ] = 1; 
// 			priorityScores.pop();
// 			counter--;
// 		}
// 		//----------------------------------------------------------------------
// 	}

// 	return 0;
// }

// //cout << "pip_1" << endl;
// 		//For every moving object find the qid of its last tuple with PIP
// 	//	unsigned int qid = i->second->pipLS();
// //cout << "pip_2" << endl;

// 		//Fetch the last tuple's qid
// 		unsigned int qid = (*i->second->last)[ i->first ]->qid_cheat;
// //cout << "qid:\t" << qid;
// 		if ( i->second->Q->find(qid) == i->second->Q->end() )
// 			qid = 0;
// //cout << "\tqid:\t" << qid << endl;

// 		//Score Assignment
// 		if ( qid == 0 )
// 		{
// //cout << "oid:\t" << i->first << "\tqid:\t" << qid << endl;
// 			tf_score[ i->first ] = 0.0;
// 		}
// 		else if ( (qid != 0) && ((*results)[ qid ]->oids.size() == 0) )
// 		{
// //cout << "oid:\t" << i->first << "\tqid:\t" << qid << "\t|qid|:\t" << (*results)[ qid ]->oids.size();
// 			tf_score[ i->first ] = 1.0;
// //cout << "\tscore:\t" << tf_score[ i->first ] << endl;
// 			nonZero++;
// 		}
// 		else if ( (qid != 0) && ((*results)[ qid ]->oids.size() != 0) )
// 		{
// //cout << "oid:\t" << i->first << "\tqid:\t" << qid << "\t|qid|:\t" << (*results)[ qid ]->oids.size();
// 			tf_score[ i->first ] = 1.0 / (double) (*results)[ qid ]->oids.size();
// //cout << "\tscore:\t" << tf_score[ i->first ] << endl;
// 			nonZero++;
// 		}

// 		//Push the trajectory scores to the priority queue
// 	 	priorityScores.push( make_pair( i->first, tf_score[ i->first ] ) );
// 	}

// int LoadShedder::CellSynopsisLS()
// {
// 	int counter = *n; 		
// 	priority_queue< pair< unsigned int, double>, vector< pair< unsigned int, double > >, CompareDist> priorityScores;
// 	unordered_map<unsigned int, double> score;        //Cell based synopsis score for a trajectory of a moving object

// 	//Check the counters values -- for out of bounds cases --
// 	//-------------------------------------------------------
// 	if ( counter <= 0 ) //-----------------------------------
// 		counter = 0; 		
// 	else if ( counter >= T->size() )
// 	{
// 		*n 		= T->size(); 
// 		counter = T->size(); 																		 

// 		if ( isUnder == 1 )
// 		{
// 			for ( auto i = T->begin() ; i != T->end() ; i++ )
// 				(*U)[ i->first ] = 1;
// 			return -1; 																						
// 		} //-------------------------------------------------
// 	} //-----------------------------------------------------

// 	for ( auto i = T->begin() ; i != T->end() ; i++ )
// 	{
// 		if ( i->second->mrut > (double) tnow - (double) range )
// 			score[i->first] = i->second->synopsis->scoring();
// 		else //In case of inactive object
// 			score[i->first] = -10.0;
// 		//Feed the queue of scores
// 		priorityScores.push( make_pair( i->first, score[i->first] ) );
// 	}

// 	// Erase the U map	
// 	for ( auto i = T->begin() ; i != T->end() ; i++ )
// 		(*U)[ i->first ] = 0;

// 	while ( counter > 0 )
// 	{
// 		(*U)[ priorityScores.top().first ] = 1;
// 		//Print scores -----------------------------
// 		//cout << priorityScores.top().second << endl;
// 		//------------------------------------------
// 		priorityScores.pop();
// 		counter--;
// 	}

// 	return 0;
// }