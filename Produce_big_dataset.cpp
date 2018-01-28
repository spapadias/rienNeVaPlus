//Title: Create a bigger synthetic dataset from an existing one (real or synthetic)
//Author  : Papadias Serafeim
//Date    : 26/03/2015 
//Revision: 25/03/2017

// #include <algorithm>
// #include <stdio.h>
// #include <stdlib.h>
// #include <sstream> 
// #include <string.h>
// #include <sys/time.h>
// #include <sys/timeb.h>
// #include <string>
// #include <set>
// #include <fstream>
// #include <iostream>
// #include <iomanip>
// #include <cmath>
// #include <ctime>
// #include <cstddef>
// #include <map>
// #include <list>
// #include <vector>
// #include <utility>
// #include <unordered_map>
// #include <unordered_set>

#include "./Utilities/Structures.h" 
#include "./Scan/Scan.cpp"
#include "./Point/Point.cpp"

using namespace std;

int main (int argc, char *argv[])
{
	if (argc != 4)  
	{
        //EXAMPLE:  ./big_dataset VALID trajectory_file.txt 
        cout << "Usage: " << argv[0] << " VALID [input-file] 0" << endl;
        exit(0);
	}

	cout.setf(ios::boolalpha);	                      //Print "true" and "false" instead of 1, 0
    char *modeTime;                
 	modeTime = argv[1]; 

    char *fileName;                
 	fileName = argv[2];                               //ASCII file with object tuples of schema < t id x y >    

    bool mode = (strcmp(modeTime, "VALID") == 0);     //True for VALID timestamps
    //mode = 0 for VALID timestamps , else mode = 1 for TRANS timestamps

	Scan * scanStream = new Scan(fileName, mode);
    if (mode)
	    scanStream->setTimeAttribute(atoi(argv[3]));  //...the timestamp attribute in the schema of tuples
    else
	    scanStream->setRate(atoi(argv[3]));           //...the stream arrival rate (in tuples/sec)

	vector<Point*> inTuples;
    int t = 0;

    while (scanStream->exhausted == false)
    {
        //Read streaming data
        inTuples = scanStream->consumeInput(t);          
        
        for ( auto it = inTuples.begin() ; it != inTuples.end() ; it++ )
        {
            //Check if the tuples have been read 
            cout << (*it)->t << " " << (*it)->oid          << " " << setprecision(8) << fixed << (*it)->x << " " << setprecision(8) << fixed << (*it)->y << " " << (*it)->qid << endl;
            cout << (*it)->t << " " << (*it)->oid + 100000 << " " << setprecision(8) << fixed << (*it)->x << " " << setprecision(8) << fixed << (*it)->y << " " << (*it)->qid << endl; //Create the 2 * 60k 120k dataset
            cout << (*it)->t << " " << (*it)->oid + 200000 << " " << setprecision(8) << fixed << (*it)->x << " " << setprecision(8) << fixed << (*it)->y << " " << (*it)->qid << endl; //Create the 180k dataset
            cout << (*it)->t << " " << (*it)->oid + 300000 << " " << setprecision(8) << fixed << (*it)->x << " " << setprecision(8) << fixed << (*it)->y << " " << (*it)->qid << endl; //Create the 240k dataset
            cout << (*it)->t << " " << (*it)->oid + 400000 << " " << setprecision(8) << fixed << (*it)->x << " " << setprecision(8) << fixed << (*it)->y << " " << (*it)->qid << endl; //Create the 300k dataset      

            // cout << (*it)->t << " " << (*it)->oid + 500000 << " " << setprecision(8) << fixed << (*it)->x << " " << setprecision(8) << fixed << (*it)->y << " " << (*it)->qid << endl; //Create the 360k dataset      
            // cout << (*it)->t << " " << (*it)->oid + 600000 << " " << setprecision(8) << fixed << (*it)->x << " " << setprecision(8) << fixed << (*it)->y << " " << (*it)->qid << endl; //Create the 420k dataset      
            // cout << (*it)->t << " " << (*it)->oid + 700000 << " " << setprecision(8) << fixed << (*it)->x << " " << setprecision(8) << fixed << (*it)->y << " " << (*it)->qid << endl; //Create the 480k dataset      
            // cout << (*it)->t << " " << (*it)->oid + 800000 << " " << setprecision(8) << fixed << (*it)->x << " " << setprecision(8) << fixed << (*it)->y << " " << (*it)->qid << endl; //Create the 540k dataset      
            // cout << (*it)->t << " " << (*it)->oid + 900000 << " " << setprecision(8) << fixed << (*it)->x << " " << setprecision(8) << fixed << (*it)->y << " " << (*it)->qid << endl; //Create the 600k dataset      
            // cout << (*it)->t << " " << (*it)->oid + 1000000 << " " << setprecision(8) << fixed << (*it)->x << " " << setprecision(8) << fixed << (*it)->y << " " << (*it)->qid << endl; //Create the 660k dataset      
//            cout << (*it)->t << " " << (*it)->oid + 1100000 << " " << setprecision(8) << fixed << (*it)->x << " " << setprecision(8) << fixed << (*it)->y << " " << (*it)->qid << endl; //Create the 720k dataset      
//            cout << (*it)->t << " " << (*it)->oid + 1200000 << " " << setprecision(8) << fixed << (*it)->x << " " << setprecision(8) << fixed << (*it)->y << " " << (*it)->qid << endl; //Create the 780k dataset      
//            cout << (*it)->t << " " << (*it)->oid + 1300000 << " " << setprecision(8) << fixed << (*it)->x << " " << setprecision(8) << fixed << (*it)->y << " " << (*it)->qid << endl; //Create the 840k dataset      
//            cout << (*it)->t << " " << (*it)->oid + 1400000 << " " << setprecision(8) << fixed << (*it)->x << " " << setprecision(8) << fixed << (*it)->y << " " << (*it)->qid << endl; //Create the 900k dataset      
//            cout << (*it)->t << " " << (*it)->oid + 1500000 << " " << setprecision(8) << fixed << (*it)->x << " " << setprecision(8) << fixed << (*it)->y << " " << (*it)->qid << endl; //Create the 960k dataset      
//            cout << (*it)->t << " " << (*it)->oid + 1600000 << " " << setprecision(8) << fixed << (*it)->x << " " << setprecision(8) << fixed << (*it)->y << " " << (*it)->qid << endl; //Create the 1020k dataset      
        }

        inTuples.clear();

        t = t + 10;
    } //end of while loop - input file exhausted

	return 0;
}
