//Title: RegularGrid.h
//Description: Regular grid partitioning for manipulating point data streams.
//Platform(s): gcc (Ubuntu/Linaro 4.6.3-1ubuntu5) 4.6.3
//Date:     29/12/2007
//Revision: 11/05/2017

#ifndef REGULARGRID_H
#define REGULARGRID_H

#include "../Utilities/Structures.h"
#include "../Point/Point.h"

class RegularGrid 
{

public:
    RegularGrid(unsigned int &, unsigned int &);                        
    RegularGrid(unsigned int &, unsigned int &, box);                   
    ~RegularGrid();	                                                    

    bool Allocate();
    void printGridState();
    void getCellObjects(unsigned int);
    void getCellQueries(unsigned int);                                 
    inline unsigned int HashLocation(double &, double &);
    inline list<unsigned int> HashBox(box &);
    inline box getCellBox(unsigned int &);
    inline void updateObject(Point *);
    int pnpoly(int , polygon *, double, double );                       //Point in polygon algorithm 
    void objectsInQueryCells(set<unsigned int> *, unsigned int *);
    //---
    void objectsInQueryCellsHistory(set<unsigned int> *, unsigned int *);

    unordered_map<unsigned int, unsigned int> ObjAssignments;           //For each object, remember the cell it has been allocated to <oid, cid>

    double ts;                                                          //Current timestamp value (Heartbeat)
    unsigned int cell_cnt;                                              //Total number of cells
    unsigned int GranX, GranY;                                          //Granularity of each dimension for hashing
    unsigned int obj_cnt;                                               //Number of distinct objects indexed in the grid
    double XMIN, YMIN, XMAX, YMAX;                                      //Space bounds (universe of discourse)
    double width, height;                                               //Width, height of the 2D space
    double dx, dy;                                                      //Dimensions of each cell

    struct GridCell
    {
        unordered_map<unsigned int, Point* > ObjList;                   //Object locations assigned into this cell
        unordered_map<unsigned int, TQuery*> QryList;                   //Queries assigned to this cell of the grid
        box cellBox;                                                    //Cell rectangle specified by its coordinates
        //bool processed;
        //---
        unordered_set<unsigned int> HisObjList;                         //Object ids that moved in this cell inside the time window
        //---
    };

    GridCell *cell; 

private:
    
};

#endif /* REGULARGRID_H */