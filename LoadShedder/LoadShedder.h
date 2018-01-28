//Title: LoadShedder.h
//Description: Class for LoadShedder
//Author: Serafeim Papadias
//Platform(s): gcc (Ubuntu/Linaro 4.6.3-1ubuntu5) 4.6.3
//Date:     04/04/2015
//Revision: 20/01/18 (previous --> 22/03/2017)

#ifndef LOADSHEDDER_H
#define LOADSHEDDER_H

#include "../Utilities/Structures.h"
#include "../RegularGrid/RegularGrid.h"
//#include <random>

class LoadShedder
{

public:
    LoadShedder(
    RegularGrid *,
    unordered_map< unsigned int, Trajectory*  > *, 
    unordered_map< unsigned int, unsigned int > *, 
    unordered_map< unsigned int, unsigned int > *, 
    unordered_map< unsigned int, queryRes*    > *, 
    unordered_map< unsigned int, TQuery*      > *,
    int*, 
    unsigned int,
    unsigned int,
    unsigned int,
    unsigned int
    );

    ~LoadShedder();

    // Semantic Selection
    unsigned int FewObjectsFunc(unsigned int);
    int SemanticLS(unsigned int);
    
    // Simple Random Sample (SRS) -- Reservoir Sample
    // 1 -- The principle of SRS is that every object has the same probability of being chosen n/T.size
    // 2 -- IMPORTANT : Every subset of size s has the same n has the same probability to be the sample! (this does not hold for Systematic Sampling)
    // 3 -- For a small sample from a large population, sampling without replacement is approximately the same as sampling with replacement, since the odds of choosing the same individual twice is low
        // Algorithmic Sketch
    // 1 -- Maintain the trajectories with maintaining a reservoir sample --> for every new trajectory include it with a certain probability n/T.size 
    // 2 -- Gather all the ids and do not decide until the next LS phase (CAUTION: this way we maintain all trajectories)
        // If the sample has unknown size --> it still works --> if size of the reservoir changes
            // case 1 -- increase : Insert more trajectory ids into the (updated in terms of size) reservoir
            // case 2 -- decrease : Delete from the reservoir : (i) randomly, (ii) semantically
    int ReservoirSampleLS_tuple(unsigned int, vector<unsigned int> *);
    int ReservoirSampleLS_under();
    int ReservoirSampleLS_over();

    // Systematic Sampling
    // 1 -- Equal probability for each trajectory id to be selected 
    // 2 -- Different that SRS because not every sample of certain size has the same probability of being selected (e.g. samples with at lest two adjacent elements are impossible)
    // ?? It is however, much more efficient (if variance within systematic sample is more than variance of population) ??
    // 3 -- skip = T.size/n -- equiprobability of selection for each element
//    int SystematicSampleLS();

    // Permutation Sampling
    // 1 -- No all trajectories has the same probability of being in the sample of the "stream of trajectory ids" -- non-uniform sample
    // 2 -- Trajectories which are older have greater probability of being selected more times after the random permutatio
    int PermutationSampleLS();

    // Deterministic Selection
    // 1 -- Selects the first n trajectories of the unordered map (hashmap), thus it is an non-uniform sample
    int SelectFirstUnorderedLS();
    int SelectFirstOrderedLS();

    // Stratified Sampling
    // 1 -- Divide the whole trajectory ids (sampling space) into groups based on location (why? reasonable cause you know where the queries are located)
    // 2 -- Use any sample technique discussed above to maintain the samples inside each strata
    // 3 -- Idea : Hybrid Stratified Sampling --> SRS, Reservoir, Permutation, Systematic, Semantic Selection
    // CAUTION : Concern whether this is stratification or plain filtering based on location 
    // CAUTION : Location is not always known for new tuples --> approximation, heuristics, e.g. Use the last calculated location for incoming unimportant elements
    int StratifiedSampleLS();

    // Random Tuples Sampling
    // 1 -- Sample with SRS/Reservoir each location not the trajectory ids
    // 2 -- Tranlate the n proportionally to number of tuples to maintain using a reservoir from the input
    int RandomTupleLS_tuple(unsigned int, unsigned int);
    int RandomTupleLS_over();
    int RandomTupleLS_under();

    //Old ones
    int StratifiedLS();
    int StratifiedLS_2();

    //Variables
    RegularGrid                                *grid;
    unordered_map<unsigned int, Trajectory*  > *T;                            
    unordered_map<unsigned int, unsigned int > *U;
    unordered_map<unsigned int, unsigned int > *invIndex;   
    unordered_map<unsigned int, queryRes*    > *results;
    unordered_map<unsigned int, TQuery*      > *Q;
    int *n;              
    unsigned int obj_cnt;                                     
    //unsigned int boolW;
    unsigned int isUnder;
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

    class CompareDistQuantile
    {
        public:
            bool operator()(const std::pair<unsigned int, unsigned int> &left, const std::pair<unsigned int, unsigned int> &right) 
            {   
                if ( left.second > right.second )
                    return true;
                else
                    return false;
            }    
    };

private:    
    
};

#endif /* LOADSHEDDER_H */