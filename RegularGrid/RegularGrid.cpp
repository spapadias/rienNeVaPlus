//Title: RegularGrid.cpp
//Description: Regular grid partitioning for manipulating point data streams.
//Platform(s): gcc (Ubuntu/Linaro 4.6.3-1ubuntu5) 4.6.3
//Date:     29/12/2007
//Revision: 11/05/2017

#include "RegularGrid.h"

RegularGrid::RegularGrid(unsigned int &_GranX, unsigned int &_GranY)
{
    XMIN   = X_MIN; YMIN = Y_MIN; XMAX = X_MAX; YMAX = Y_MAX;
	width  = XMAX - XMIN;    
	height = YMAX - YMIN;
	GranX  = _GranX; 
	GranY  = _GranY;
	//Compute cell dimensions
	dx     = width/(double)GranX; 
	dy     = height/(double)GranY;
}

RegularGrid::RegularGrid(unsigned int &_GranX, unsigned int &_GranY, box pBox)
{
    XMIN   = pBox.x_min; YMIN = pBox.y_min; XMAX = pBox.x_max; YMAX = pBox.y_max;
	width  = XMAX - XMIN;    
	height = YMAX - YMIN;
	GranX  = _GranX; 
	GranY  = _GranY;
	//Compute cell dimensions
	dx     = width /(double)GranX; 
	dy     = height/(double)GranY;
}

RegularGrid::~RegularGrid()
{
    //No need to maintain information about object allocations anymore
    ObjAssignments.clear();

	//Destroy contents for all cells
	for( unsigned int cid = 0 ; cid < cell_cnt ; cid++ ) 
    {
        //Delete Objects
        for ( auto i = cell[cid].ObjList.begin() ; i != cell[cid].ObjList.end() ; i++ )
            delete i->second;
		cell[cid].ObjList.clear(); 			
		
        //Delete Queries
        for ( auto i = cell[cid].QryList.begin() ; i != cell[cid].QryList.end() ; i++ )
        {
            delete i->second->startV;
            delete i->second->endV;
            delete i->second->qArea;
            delete i->second;
        }
        cell[cid].QryList.clear(); 			
	}

	delete[] cell;
}

bool RegularGrid::Allocate() 
{
	cell_cnt = GranX * GranY;
	try 
	{
		cell = new GridCell[cell_cnt];					//Allocate Cell table

		//Initialize Cell table
		for ( unsigned int i = 0 ; i < cell_cnt ; ++i ) 
        {
           	cell[i].cellBox = getCellBox(i);     		//Store cell coordinates
            //cell[i].processed = false;
        }
	}
	catch(...) { return false; }

	return true;
}

inline unsigned int RegularGrid::HashLocation(double &x, double &y)
{
    return GranX * (unsigned int)( (y-YMIN)/dy ) + (unsigned int)( (x-XMIN)/dx ); 
}

inline list<unsigned int> RegularGrid::HashBox(box &pBox)
{
    list<unsigned int> cellList;               			 //List of cells returned with indicators for partial overlap
   
    unsigned int lc, uc, c, i, gdx;

    //First, hashing box corners
    lc = HashLocation(pBox.x_min, pBox.y_min); 
    uc = HashLocation(pBox.x_max, pBox.y_max);
        
    //Range of cells affected in each row of the grid
    gdx = uc%GranX - lc%GranX;
    
    //Find all cells covered by this box after examining the cells of its corners
    for(c=lc; c<=uc-gdx; c+=GranX)  
        for(i=c; i<=c+gdx; ++i)     
        {   
            if (i >= cell_cnt) 							 //i is the cell currently being examined for overlap
                continue;
                    
            //Insert this cell into the list    
            cellList.push_back( i );
        }

    return cellList; 
}

inline box RegularGrid::getCellBox(unsigned int &cellID)
{
    //Cell matrix indices along x and y axes
    unsigned int i = cellID % GranX;   
    unsigned int j = cellID / GranX;
    
    double x_min = XMIN + i * dx;
    double y_min = YMIN + j * dy;
    double x_max = XMIN + (i+1) * dx;
    double y_max = YMIN + (j+1) * dy;

    return box(x_min, y_min, x_max, y_max);
}

void RegularGrid::getCellObjects(unsigned int cid)
{
    if ( cell[cid].ObjList.size() > 0 )
    {
        cout << "cid:\t" << cid << "\toids:\t";
	    for( auto i = cell[cid].ObjList.begin(); i != cell[cid].ObjList.end(); i++ ) 
    	{
    		Point *o = i->second;
            cout << o->oid << ", ";
    	}
        cout << endl;
    }
}

void RegularGrid::getCellQueries(unsigned int cid)
{
	if ( cell[cid].QryList.size() > 5 ) //The cells with 0 queries as well  
	{
        //Format cid
        cout << cid << "\t" << cell[cid].QryList.size() << endl;
	}
}

void RegularGrid::printGridState()
{
	for ( unsigned int i = 0 ; i < cell_cnt ; ++i ) 
	{
			getCellObjects(i);  					//Currently retained object locations
//			getCellQueries(i);  					//Currently retained query  locations
	} 	
}

int RegularGrid::pnpoly(int nvert, polygon *pPol , double testx, double testy)
{
    double vertx[nvert];
    double verty[nvert];
    
    vector<vertex *>::iterator vertIt;
    unsigned int counter = 0;
    
    for ( vertIt = pPol->vertices.begin() ; vertIt != pPol->vertices.end() ; vertIt++ )
    {
        vertx[counter] = (*vertIt)->x;
        verty[counter] = (*vertIt)->y;
        counter++; 									//Increase the counter
    }

    int i, j, c = 0;
    for (i = 0, j = nvert-1; i < nvert; j = i++) {
        if ( ((verty[i]>testy) != (verty[j]>testy)) && (testx < (vertx[j]-vertx[i]) * (testy-verty[i]) / (verty[j]-verty[i]) + vertx[i]) )
            c = !c;
    }
    return c;
    //CAUTION: This code gives 0 for the point being at the boundary of the Polygon
}

inline void RegularGrid::updateObject(Point *p)
{
    unsigned int i;

    //Assume ascending timestamp ordering: if update ts < current ts, reject this object update
    if (p->t < ts)
        return;
    
    //Remove this object from the cell it was previously allocated to if it is needed.
    if ( ObjAssignments.find(p->oid) != ObjAssignments.end() )
    {
        i = ObjAssignments[p->oid];                 //Cell id
        cell[i].ObjList.erase(p->oid);              //Delete object from that cell
    }
    
    //Hashing object's positional update to find the proper cell
    i = HashLocation(p->x, p->y);

    //Assign new object location into the cell just found
    cell[i].ObjList.insert( {p->oid, p} );
    //cell[i].processed = false;                    //Mark that this cell must be probed by query reevaluation
    ObjAssignments[p->oid] = i;                     //Remember this object assignment for subsequent execution cycles

}

void RegularGrid::objectsInQueryCells(set<unsigned int> *cellsToSample, unsigned int *totalObj)
{
    //---
    *totalObj  = 0; //In cells with queries
    unsigned int cellWithQ = 0;

    for ( unsigned int i = 0 ; i < cell_cnt ; ++i )
    {
        //---
        if ( !cell[i].QryList.empty() )
        {
            if ( !cell[i].ObjList.empty() )
            {
                //Form the total
                *totalObj += cell[i].ObjList.size();
                cellsToSample->insert(i);
            }
            //---    
            cellWithQ++;
        }
        //---
    }
//    cout << "Percentage of cells with queries:\t" << (cellWithQ * 100.0) / cell_cnt << "\%" << endl;
    //---
}

void RegularGrid::objectsInQueryCellsHistory(set<unsigned int> * cellsToSample, unsigned int *totalObj)
{
    //---
    *totalObj = 0; //In cells with queries
    unsigned int cellWithQ = 0;

    for ( unsigned int i = 0 ; i < cell_cnt ; ++i )
    {
        //---
        if ( !cell[i].QryList.empty() )
        {
            if ( !cell[i].HisObjList.empty() )
            {
                //Form the total
                *totalObj += cell[i].HisObjList.size();
                cellsToSample->insert(i);
            }
            //---    
            cellWithQ++;

        }
        //---
    }
    //    cout << "Percentage of cells with queries:\t" << (cellWithQ * 100.0) / cell_cnt << "\%" << endl;
    //---
}

