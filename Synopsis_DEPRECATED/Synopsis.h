//Title: Synopsis.h
//Description: This file defines the class Synopsis
//Platform(s): gcc (Ubuntu/Linaro 4.6.3-1ubuntu5) 4.6.3
//Date: 	04/04/2015
//Revision: 19/05/2017

#ifndef SYNOPSIS_H
#define SYNOPSIS_H
//------------------------------------------------------------
#include "../Utilities/Structures.h"
#include "../RegularGrid/RegularGrid.h"
#include "../Point/Point.h" 								 
//------------------------------------------------------------

class Trajectory;
//------------------------------------------------------------

class Synopsis
{

public:

	Synopsis(RegularGrid *, Trajectory *); 											
	~Synopsis(); 					

	void   updateQ(Point *); // Update Qynopsis
	void   updateC(Point *); // Update Cynopsis
    void   cropQ(double);
    void   cropC(double);

    double scoringC(); 	
    double scoringQ();	

    int    sizeQ();
    int    sizeC();

	void   print();

	typedef struct partQID
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

	typedef struct partCID
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

    //Variables
    RegularGrid *grid;	     // Pointer to the grid
    Trajectory  *traj;		 // Pointer to the trajectory 	
	list<synPartQ> Qynopsis; // Query ID Synopsis
	list<synPartC> Cynopsis; // Cell  ID Synopsis										
	double duration;	 	 // Total duration of the synopsis

private:

};

#endif /* SYNOPSIS */
