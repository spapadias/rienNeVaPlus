//Title: Synopsis.h
//Description: Synopsis of a given trajectory of a moving object on a given grid.
//Platform(s): gcc (Ubuntu/Linaro 4.6.3-1ubuntu5) 4.6.3
//Date: 4/4/2015
//Revision: 4/3/2017

#ifndef SYNOPSIS_H
#define SYNOPSIS_H //Inclusion guards

#include "RegularGrid.h" 								 

//#include "Trajectory.h" 	//The forward declaration needed here! Because there is a cycle between Trajectory.h and Synopsis.h
class Trajectory; 			/* Forward declaration of Trajectory */

class Synopsis
{

public:
	Synopsis(Trajectory *); 							
	~Synopsis(); 					

	//Functions of the Synopsis class
	double updateSynopsis(point *, unordered_map< unsigned int, unsigned int > *, int, int, unordered_map< unsigned int, unsigned int > *, unordered_set<unsigned int> *);
    void   cropSynopsis(double); 							
    double findTimeIntervalSynopsisFast(unsigned int); 	//Fast retrieval of time interval -- O(1)
    double scoring(unsigned int); 		
    double scoringTF(list< map < unsigned int, long double >* >::iterator); //Calculates the TF score of this trajectory
    int    sizeSynopsis();		 //Size of the synopsis in the interval [t_NOW - maxWinRange , t_NOW]
	void   printSynopsis(); 	

    //Define the atom of a synopsis of a trajectory (queryID, [t1,t2])
	typedef struct PartOfSynospis
	{
	    unsigned int qid; 							
	   	double t1; 										
	   	double t2; 										

	   	double getTimeDif()
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
	} synopart;

    struct classcomp
    {
		bool operator()(const synopart &lhs, const synopart &rhs) const
		{
			return lhs.qid < rhs.qid;
		}
	};

    //Variables of the class Synopsis
	list<synopart> synopsis; 								
	unordered_multimap<unsigned int, synopart> synopsisHash;	//Calculates the interval in O(1)
	Trajectory *traj; 											
	double total_duration; 										

private:

};

#endif /* SYNOPSIS */
