//Title: Synopsis.h
//Description: This file defines the class Synopsis
//Platform(s): gcc (Ubuntu/Linaro 4.6.3-1ubuntu5) 4.6.3
//Date: 	04/04/2015
//Revision: 19/05/2017

#ifndef SYNOPSIS_H
#define SYNOPSIS_H
//------------------------------------------------------------
#include "./Utilities/Structures.h"
#include "../RegularGrid/RegularGrid.h" 								 
//------------------------------------------------------------

class Trajectory;
//------------------------------------------------------------

class Synopsis
{

public:
	Synopsis(Trajectory *, unsigned int); 											
	~Synopsis(); 					

	double update(Point *, unordered_map< unsigned int, unsigned int > *, int, int, unordered_map< unsigned int, unsigned int > *, unordered_set<unsigned int> *);
    void   crop(int); 	
    double retInterval(unsigned int); 										//Fast retrieval of time interval -- O(1)
    double scoring(unsigned int); 		
    double scoringTF(list< map < unsigned int, long double >* >::iterator); //Calculates the TF score of this trajectory
    int    size();		 													//Size of the synopsis in [t_NOW - maxWinRange , t_NOW]
	void   print(); 	

	typedef struct partCELL
	{
	    unsigned int cid; //cell ID 							
	   	double t1; 										
	   	double t2; 										

	   	double getInt()
	   	{
	   		return (t2 - t1 + 1.0);
	   	}

	   	double getTin()
	   	{
	   		return t1;
	   	}

	   	double getTout()
	   	{
	   		return t2;
	   	}
	} synPartC;

	typedef struct partQUERY
	{
	    unsigned int qid; //query ID				
	   	double t1; 										
	   	double t2; 										

	   	double getInt()
	   	{
	   		return (t2 - t1 + 1.0);
	   	}

	   	double getTin()
	   	{
	   		return t1;
	   	}

	   	double getTout()
	   	{
	   		return t2;
	   	}
	} synPartQ;

    //Variables
    unsigned int type; 	     //(pointer to trajectory, synopsis type (1 or 2)), 1 : cell based & 2 : query based			
    Trajectory *traj;		 //Pointer to the trajectory 	
	list<synPartQ> Qynopsis; //Query ID synopsis
	list<synPartC> Cynopsis; //Cell  ID Synopsis										
	double duration;	 	 //Total duration of the synopsis
	//----------------------							
	unordered_multimap<unsigned int, synopart> synopsisHash;	//Calculates the interval in O(1)

private:

};

#endif /* SYNOPSIS */
