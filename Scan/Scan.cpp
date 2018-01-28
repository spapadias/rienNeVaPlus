//Title: Scan.cpp 
//Description: Consumes input lines from ASCII file (a) according to the specified arrival rate. 
//OR (b) for a specified timestamp value in the dataset. 
//No tuple manipulation or timestamp assignment is done at that stage.
//Author: Kostas Patroumpas, Papadias Serafeim
//Platform(s): gcc (Ubuntu/Linaro 4.6.3-1ubuntu5) 4.6.3
//Date:     07/10/2009
//Revision: 15/03/2017

#include "Scan.h" 

//Constructor #1 for reading tuples according to the specified mode
Scan::Scan(char *fileName, bool mode)
{
	fin.open(fileName, ios::in);
	this->recCount = 0;
	this->exhausted = false;
    this->curTime = 0.0;
    this->mode = mode;
    this->inLine = "";
}

//Constructor #2 for reading tuples with a fixed arrival rate
Scan::Scan(char *fileName, unsigned int inRate)
{
	fin.open(fileName, ios::in);
	this->recCount = 0;
	this->inRate = inRate;
	this->exhausted = false;
    this->curTime = 0;
}

//Constructor #3 for reading tuples according to incoming timestamp values
Scan::Scan(char *fileName, unsigned int attrTime, unsigned int t)
{
	// Dummy variable
	unsigned int temp_t = t;
	fin.open(fileName, ios::in);
	this->recCount = 0;
	this->attrTime = attrTime;
	this->inLine = "";
	this->exhausted = false;
	temp_t = attrTime;
}

//Destructor
Scan::~Scan()
{
	fin.close();
}

void Scan::setRate(unsigned int inRate)
{
	this->inRate = inRate;
}

void Scan::setTimeAttribute(double attrTime)
{
	this->attrTime = attrTime;
}


//Depending on the mode, it calls a specific function to read input tuples
vector<Point*> Scan::consumeInput(unsigned int t)
{
    batchTuples.clear(); 			//Clears the buffer that holds the tuples read at the previous step

    if (this->mode)
        this->read(t); 				//VALID timestamping
    else
        this->read();  				//TRANSACTION timestamping     

    return batchTuples;
}

//Decode tuple attributes from incoming string value
Point* Scan::decodeTuple(fstream &fin)
{
	Point *point = new Point();
	fin >> point->t >> point->oid >> point->x >> point->y >> point->qid; //Schema < t oid x y qid >

	//-- CAUTION -- Cheating Here -----
	point->qid_cheat = point->qid; //--
	//---------------------------------

    return point;
}

//FIRST OPTION --> TRANS: Read a batch of lines from the input ASCII file representing a streaming source.
void Scan::read()
{
	//Handle input source according to the specified arrival rate
	do
	{
        if (!fin.eof())
		{
            inTuple = this->decodeTuple(fin);
			batchTuples.push_back(inTuple);		//Create a batch of incoming tuples 
			recCount++;
		}
		else
		{
			this->exhausted = true;				//EOF
			break;
		}
	} while ( recCount % inRate != 0);

    //Next timestamp
    this->curTime++;

	//Notify progress
	//cerr << recCount << " records processed..." << "\r";
}

//SECOND OPTION --> VALID: Read a batch of lines from the input ASCII file representing a streaming source 
//                  until the specified timestamp value.
void Scan::read(unsigned int t) //This must be unsigned int
{
	//First return the tuple that had been prefetched in the previous cycle
	if (inTuple != NULL)
    {
        if (inTuple->t <= t)  //Only in case it fits within the upper window bound
        {
			batchTuples.push_back(inTuple);
            recCount++;
        }
        else                //No need to consume more tuples, as the upper window bound has not reached the next timestamp value in input
            return;
    }

	//Handle input source according to the specified arrival rate
	do
	{
		if (!fin.eof())
		{
            inTuple = this->decodeTuple(fin);
            //Current timestamp refers to the one from the last accessed tuple
			this->curTime = inTuple->t;
			if (this->curTime > t)	//Exceeded timestamp limit
				break;

			batchTuples.push_back(inTuple);		//Create a batch of incoming tuples 
			inTuple = NULL;
            recCount++;
		}
		else
		{
			this->exhausted = true;				//EOF
			break;
		}
	} while (true);

	//Notify progress
	//	cerr << t << " -> " << recCount << " records processed..." << "\r";
}
