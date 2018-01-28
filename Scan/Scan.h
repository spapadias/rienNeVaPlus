//Title: Scan.h
//Description: Scan class to read an input data stream.
//Author: Kostas Patroumpas, Papadias Serafeim
//Platform(s): gcc (Ubuntu/Linaro 4.6.3-1ubuntu5) 4.6.3
//Date:     07/10/2009
//Revision: 19/03/2017

#ifndef SCAN_H
#define SCAN_H

#include "../Utilities/Structures.h"
#include "../Point/Point.h"

//class Point;

//Class for maintaining joined items from the windowing constructs
class Scan 
{

public:
    Scan(char*, bool);
	Scan(char*, unsigned int);
	Scan(char*, unsigned int, unsigned int); //CAUTION: t is not used. It differenciates the constructors
	~Scan();
    void setRate(unsigned int);
    void setTimeAttribute(double);
	vector<Point *> consumeInput(unsigned int);
	bool exhausted;									//Set TRUE at EOF
	unsigned int recCount;							//Count incoming tuples
    double curTime;
    bool mode;

private:
	fstream fin;
	unsigned int inRate;
	string inLine;
    Point* inTuple;
	double attrTime;
	void read();
	void read(unsigned int);
    vector<Point*> batchTuples;
    Point* decodeTuple(fstream &);
};

#endif /* SCAN */
