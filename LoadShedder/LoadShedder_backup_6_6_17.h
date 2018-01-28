//Title: LoadShedder.h
//Description: Class for LoadShedder
//Author: Serafeim Papadias
//Platform(s): gcc (Ubuntu/Linaro 4.6.3-1ubuntu5) 4.6.3
//Date:     04/04/2015
//Revision: 22/03/2017

#ifndef LOADSHEDDER_H
#define LOADSHEDDER_H

#include "../Utilities/Structures.h"
#include "../RegularGrid/RegularGrid.h"
//#include <random>
//using namespace std;

class LoadShedder
{

public:
    LoadShedder(
                unordered_map< unsigned int, Trajectory* > *, 
                unordered_map< unsigned int, unsigned int > *, 
                unordered_map< unsigned int, unsigned int > *, 
                unordered_map<unsigned int, queryRes*> *, 
                unordered_map<unsigned int, TQuery*      > *Q,
                int*, 
                unsigned int,
                unsigned int,
                unsigned int
               );

    ~LoadShedder();

    //Sampling Based Techniques
    int RandomLS_Shuffle();
    int StratifiedLS();




    //int RandomLS_FisherYates();
    //int QueriesPerObject();
    //int CombinationLS();    

    //2 shades of Random
    //2 shades of Trajectory Frequency (TF)
    //int TF_noHistory(unsigned int, double);
    //int TF_history();
    //1 Shade of Synopsis
    //int CellSynopsisLS();
    //Inverted Index

    //Variables
    unordered_map<unsigned int, Trajectory*  > *T;                            
    unordered_map<unsigned int, unsigned int > *U;
    unordered_map<unsigned int, unsigned int > *invIndex;   
    unordered_map<unsigned int, queryRes*    > *results;
    unordered_map<unsigned int, TQuery*      > *Q;
    int *n;                                                   
    //unsigned int boolW;
    unsigned int isLoadUnderBand;
    unsigned int tnow;
    unsigned int range;

    class CompareDist
    {
        public:
            bool operator()( const std::pair<unsigned int, double> &left, const std::pair<unsigned int, double> &right ) 
            {   
                if ( left.second < right.second )
                    return true;
                else
                    return false;
            }    
    };

    // class CompareDistQuantile
    // {
    //     public:
    //         bool operator()( const std::pair<unsigned int, unsigned int> &left, const std::pair<unsigned int, unsigned int> &right ) 
    //         {   
    //             if ( left.second > right.second )
    //                 return true;
    //             else
    //                 return false;
    //         }    
    // };

private:    
    
};

#endif /* LOADSHEDDER_H */