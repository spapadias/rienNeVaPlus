//Title: queryRetrieve.h
//Description: A class which retrieves the queries from a given query input file and stores the needed query information. 
//Author: Serafeim Papadias
//Platform(s): gcc (Ubuntu/Linaro 4.6.3-1ubuntu5) 4.6.3
//Date:     10/02/2017
//Revision: 20/05/2017

#ifndef QUERYRETRIEVER_H
#define QUERYRETRIEVER_H

#include "../Utilities/Structures.h"
#include "../RegularGrid/RegularGrid.h"

class ReadQueries
{

public:
    ReadQueries(unsigned int,
                unsigned int, 
                unsigned int, 
                unordered_map<unsigned int, TQuery*> *, 
                RegularGrid *, 
                vector<unsigned int> *);
    ~ReadQueries();

    //Functions
    int queryRead(char *, unsigned int);

    //Variables
    unsigned int isGlobalWin;
    unsigned int global_range;
    unsigned int global_slide;
    RegularGrid *grid;
    unordered_map<unsigned int, TQuery*> *Q;
    vector<unsigned int> *randomQvector; 

private:    
    
};

#endif /* QUERYRETRIEVER_H */