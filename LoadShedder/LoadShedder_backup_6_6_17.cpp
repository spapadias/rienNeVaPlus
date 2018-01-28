//Title: LoadShedder_priority_queue_scoring.cpp
//Description: A class of LoadShedder which applies the two criterions for choosing the significant trajectories
//Author: Serafeim Papadias
//Platform(s): gcc (Ubuntu/Linaro 4.6.3-1ubuntu5) 4.6.3
//Date:     04/04/2015
//Revision: 22/03/2017

#include "LoadShedder.h"

LoadShedder::LoadShedder(
	 					 unordered_map< unsigned int, Trajectory* >  *_T, 
					 	 unordered_map< unsigned int, unsigned int > *_U, 
						 unordered_map< unsigned int, unsigned int > *_invIndex, 
					     unordered_map< unsigned int, queryRes*>     *_results,
	            	     unordered_map< unsigned int, TQuery*>       *_Q, 
					     int *_n,
					     //unsigned int _boolW, 
					     unsigned int _isLoadUnderBand,
					     unsigned int _tnow,
					     unsigned int _range
						)
{
	T 				= _T;
	U               = _U;
	invIndex 		= _invIndex;
	results         = _results;
	Q               = _Q;
	n   	   		= _n;
	isLoadUnderBand = _isLoadUnderBand;
	tnow            = _tnow;
	range  			= _range;
}

LoadShedder::~LoadShedder()
{
}

int LoadShedder::RandomLS_Shuffle()
{
	int counter = *n;		  //Num Trajs for the LS to maintain
	vector<unsigned int> ids;

	//-----------------
	if ( counter <= 0 )
		counter = 0; 		
	else if ( counter >= T->size() )
	{
		*n       = T->size();
		counter  = T->size(); 																		 

		//Check if we are below or above the band -- If below then return
		if ( isLoadUnderBand == 1 )
		{
			for ( auto i = T->begin() ; i != T->end() ; i++ )
				(*U)[ i->first ] = 1;
			return -1; 																						
		}
	} //-----------------------------------------------------

	//Create a vector with all object/trajectory ids
	for ( auto i = T->begin() ; i != T->end() ; i++ )
		ids.push_back( i->first );

	//Perform a random permutation of the ids of the existing moving obejcts
	srand(time(NULL));
	random_shuffle( ids.begin(), ids.end() );

	//Randomly select #ofTrajectoriesLastlyKept trajectories & update U
	for( auto i = ids.begin() ; i != ids.end() ; i++ )
	{
		(*U)[ *i ] = 1;
		//Check how many trajectories have we kept
		counter--;
		if ( counter <= 0 )
			break;
	}

	// unsigned int picked = 0;
	// for ( auto i = U->begin() ; i != U->end() ; i++ )
	// 	if ( i->second == 1 )
	// 		picked++;
	// cout << "#traj selected:\t" << picked << endl;

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

// 		if ( isLoadUnderBand == 1 )
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

// 		if ( isLoadUnderBand == 1 )
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
// 		if ( isLoadUnderBand == 1 )
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

// 		if ( isLoadUnderBand == 1 )
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
// 		if ( isLoadUnderBand == 1 )
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

// 		if ( isLoadUnderBand == 1 )
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