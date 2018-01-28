//Title: subsetRandomItems.cpp
//Description: Create a subset of items, by randomly choosing a given amount of them.
//Author: Kostas Patroumpas
//Platform(s): gcc 4.8.3
//Date: 26/6/2011
//Revision: 15/6/2015

#include <stdio.h>
#include <stdlib.h>
//#include <math.h>
#include <iostream>
#include <set>
//#include <string>
//#include <ctime>

using namespace std;

int main(int argc, char* argv[])
{		
	unsigned int item_cnt, subset_cnt;
	unsigned int i;
	unsigned int N;
	double percent;

	if (argc<3) 
	{
		//EXAMPLE: ./subset 100 20
		cout << "\nUsage: " << argv[0] << " [Total-items-count] [Subset-count]" << endl; 
		return 1;
	}

	//Interpretation of given parameters
	item_cnt = atoi(argv[1]);	//Total count of input items
	subset_cnt = atof(argv[2]);	//Cardinality of the subset that will be preserved in output
	
	if (subset_cnt > item_cnt)
	{ 
	    cout << "INVALID PARAMETERS! Subset count must be less than the amount of input items!" << endl;
		return 1;
	}
	
	percent = subset_cnt / (double) item_cnt; //Percentage of items that will be preserved in output
	
	N = item_cnt + 1;           //Position of last item (CAUTION!!! starting enumeration of items from 1, not 0)  

	//Initialise memory array to hold input data
	int InputData[item_cnt+1];
    for (i=0; i<=N; i++){
        InputData[i] = i;
    }
	
	//The set that will hold output data
	std::set<int> OutputData;                 //CAUTION: Using set in order to avoid duplicate items
	std::set<int>::iterator it;
	
	srand(unsigned(time(NULL)));		//Initialize random generator

	//Randomly choose items that will be retained
	i = 1;   //Start from the first input item   (CAUTION!!! starting enumeration of items from 1, not 0)  
	while (OutputData.size() < subset_cnt)  //Continue as long as the required amount of output items is not reached
	{
		//Probabilistically determine which object to keep
		if (((float) (rand() % 100 + 1)/100)<= percent)    
		{
			OutputData.insert(i);
		}
		i++;
		if (i > item_cnt)    //Start again in case the required amount is not reached
		   i = 1;
	}
	
	//Print retained items
	cout << OutputData.size() << " items RETAINED:";
	for (it=OutputData.begin(); it!=OutputData.end(); ++it)
        cout << " " << *it;
	cout << endl;
		
	cout << "**************************END***************************" << endl;
	return 0;
}
