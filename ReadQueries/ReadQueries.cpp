//Title: queryRetrieve.cpp
//Description: A class which retrieves the queries from a given query input file and stores the needed query information. 
//Author: Serafeim Papadias
//Platform(s): gcc (Ubuntu/Linaro 4.6.3-1ubuntu5) 4.6.3
//Date:     10/02/2017
//Revision: 20/05/2017

#include "ReadQueries.h"

ReadQueries::ReadQueries(unsigned int _isGlobalWin,
                        unsigned int _global_range, 
                        unsigned int _global_slide, 
                        unordered_map<unsigned int, TQuery*> *_Q, 
                        RegularGrid *_grid, 
                        vector<unsigned int> *_randomQ)
{
    isGlobalWin   = _isGlobalWin;
    global_range  = _global_range;
    global_slide  = _global_slide;
    grid          = _grid;
    Q             = _Q;
    randomQvector = _randomQ;
}

ReadQueries::~ReadQueries()
{
    
}

int ReadQueries::queryRead(char *QueryFileName, unsigned int isRealDataset)
{
    //Variables
    fstream fin;
    fin.open( QueryFileName, ios::in );                                    
    string geom, X_Ystart, X_Yend;                                               
    list<unsigned int> ListOfCells;
    vector<string> orientationVerticesStart, orientationVerticesEnd;
    double x1, y1, x2, y2;
    unsigned int range = 0, slide = 0, qid, t_input, cid = 0;

//    cout << "Parsed Started" << endl;

    // Put this->randomQ into an unordered_map
    unordered_set<unsigned int> randomQ;
    for (auto i = randomQvector->begin() ; i != randomQvector->end() ; i++)
        randomQ.insert(*i);

//    cout << "Parsed random Q into a hashmap" << endl;
//    cout << "Number of Queries: " << randomQ.size() << endl;

    //Scan the Athens - synthetic query file
    if ( isRealDataset == 1 )
    {
        unsigned int cid;
        //Scan the Rome -- real dataset
        while(fin >> t_input >> qid >> cid >> X_Ystart >> X_Yend)              
        {
            // if ( (cid == 1) || (cid == 3) || (cid == 5) || (cid == 7) || (cid == 9) )
            // {
            //     cout << qid << endl;
            // }

            //Check if inside the randomQ set of queries to receive answer
            if ( randomQ.find(qid) == randomQ.end() )
            {
                std::getline(fin, geom, '\n'); //Read the rest of the line
                continue;                      //and continue to the next iteration
            }

            //Create a new query
            TQuery *query = new TQuery();
            //Schema: < t qid class_id X_start,Y_start X_end,Y_end buffer_polygon >
            query->t      = t_input; 
            query->qid    = qid; 
            query->cid    = cid; 
            
            //Distinguish cases
            if (isGlobalWin == 1)
            {
                query->range  = global_range; 
                query->slide  = global_slide;
            }
            else
            {
                query->range  = range; //CAUTION : RANGE IS NOT PROVIDED FOR REAL ROME DATASETS
                query->slide  = slide;   
            } //--------------

            //Define the v1 = (x1,y1) and v2 = (x2,y2) vertices
            orientationVerticesStart.clear(); orientationVerticesStart = split(X_Ystart, SEPARATOR);
            orientationVerticesEnd.clear()  ; orientationVerticesEnd   = split(X_Yend  , SEPARATOR);
        
            //Initialize min/max values according to coordinates of the first vertex
            x1 = atof(orientationVerticesStart[0].c_str()); y1 = atof(orientationVerticesStart[1].c_str());
            x2 = atof(orientationVerticesEnd[0].c_str())  ; y2 = atof(  orientationVerticesEnd[1].c_str());
            vertex *v1 = new vertex(x1, y1); query->startV = v1; 
            vertex *v2 = new vertex(x2, y2); query->endV   = v2;

            //Read the rest of the line
            std::getline(fin, geom, '\n');

            std::pair<polygon*, polygon*> res = constructPolygon(geom);         
            polygon *p   = res.first; query->qArea = p;
            polygon *MBB = res.second;                                          
            box *MBBox = new box( MBB->vertices[0]->x , MBB->vertices[0]->y , MBB->vertices[1]->x , MBB->vertices[1]->y );

            //Hash the MBB to see at which grid cell it falls into
            ListOfCells = grid->HashBox( *MBBox );
            for ( auto ListIt = ListOfCells.begin() ; ListIt != ListOfCells.end() ; ListIt++ )
            {
                //Put the qid into the appropriate QueryList of the corresponding cells that the MBB intersects.
                if ( grid->cell[*ListIt].QryList.find(query->qid) == grid->cell[*ListIt].QryList.end() )
                    grid->cell[*ListIt].QryList.insert( {query->qid, query} );
            }

            //Insert the query into the map structure which keeps all the queries
            if ( Q->find(query->qid) == Q->end() )
                Q->insert( {query->qid, query} );     
        }
        //End of real dataset
    }
    else
    {
        //Scan the Athens - synthetic dataset
        while(fin >> t_input >> qid >> range >> slide >> X_Ystart >> X_Yend)    
        {
            //Check if inside the randomQ set of queries to receive answer
            if ( randomQ.find(qid) == randomQ.end() )
            {
                std::getline(fin, geom, '\n'); //Read the rest of the line
                continue;                      //and continue to the next iteration
            }

            //Create a new query
            // TQuery *query = new TQuery();
            // //Schema: < t qid class_id X_start,Y_start X_end,Y_end buffer_polygon >
            // query->t      = t_input; query->qid = qid; query->cid = cid; 
            // //Currently all queries have same range & slide.
            // query->range  = global_range; query->slide = global_slide;

            //Create a new query
            TQuery *query = new TQuery();
            query->t      = t_input; 
            query->qid    = qid; 
            query->cid    = cid;
// /**/        query->range  = range; //CAUTION : We are using the different windows for each query 
// /**/        query->slide  = slide; 

            //Distinguish cases
            if (isGlobalWin == 1)
            {
                query->range  = global_range; 
                query->slide  = global_slide;
            }
            else
            {
                query->range  = range; //CAUTION : RANGE IS NOT PROVIDED FOR REAL ROME DATASETS
                query->slide  = slide;   
            } //--------------

            //Define the v1 = (x1,y1) and v2 = (x2,y2) vertices for calculating the orientation vector of each query
            orientationVerticesEnd.clear(); orientationVerticesStart.clear(); 
            orientationVerticesStart = split(X_Ystart, SEPARATOR);
            orientationVerticesEnd   = split(X_Yend  , SEPARATOR);
        
            //Initialize min/max values according to coordinates of the first vertex
            x1 = atof(orientationVerticesStart[0].c_str()); y1 = atof(orientationVerticesStart[1].c_str());
            x2 =   atof(orientationVerticesEnd[0].c_str()); y2 = atof(orientationVerticesEnd[1].c_str());
            vertex *v1 = new vertex(x1, y1); vertex *v2 = new vertex(x2, y2);

            query->startV = v1; 
            query->endV = v2;

            //Read the rest of the line
            std::getline(fin, geom, '\n');

            std::pair<polygon*, polygon*> res = constructPolygon(geom); 
            polygon *p   = res.first;                                             
            polygon *MBB = res.second;                                          

            query->qArea = p;

            //Create an instance of a box structure (x_min,y_min,x_max,y_max) 
            box *MBBox = new box( MBB->vertices[0]->x, MBB->vertices[0]->y, MBB->vertices[1]->x, MBB->vertices[1]->y );



            //Hash the MBB to see at which grid cell it falls into
            ListOfCells = grid->HashBox( *MBBox );
            for ( auto ListIt = ListOfCells.begin() ; ListIt != ListOfCells.end() ; ListIt++ )
            {
                //Put the qid into the appropriate QueryList of the corresponding cells that the MBB intersects.
                if ( grid->cell[*ListIt].QryList.find(query->qid) == grid->cell[*ListIt].QryList.end() )
                    grid->cell[*ListIt].QryList.insert( {query->qid, query} );
            }

            //Insert the query into the map structure which keeps all the queries
            if ( Q->find(query->qid) == Q->end() )
                Q->insert( {query->qid, query} );  
        }
        //End of synthetic dataset
    }
    
    fin.close();

    return 0;
}