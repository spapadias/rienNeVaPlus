//Title: linearRegression.cpp 
//Description: Function which performs linear regression
//Author: Serafeim Papadias
//Platform(s): gcc (Ubuntu/Linaro 4.6.3-1ubuntu5) 4.6.3
//Date: 08/6/2016
//Revision: 08/6/2016

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <list>
#include <map>
#include <unordered_map>
#include <sys/time.h>
#include <set>
#include <iomanip>
#include <sstream> 
#include <algorithm>
#include <numeric>

//using namespace std;

//Global variables
double slope = 0.0;
double intersept = 0.0;

typedef struct regValues 
{
    double slope;
    double intersept;

    //Constructors for regression values
    regValues(const double x1, const double y1) 
    {
        slope     = x1;
        intersept = y1;
    }

    //Empty constructor
    regValues() {};
} regVals;

regVals linearReg(const vector<double> &x, const vector<double> &y)
{
    //Define the regVals that you will return as output
    regVals ret;

    if ( x.size() != y.size() )
    {
        cout << "vectors of different size" << endl;
        //throw exception("...");
    }

    double n = x.size();

    double avgX = accumulate(x.begin(), x.end(), 0.0) / n; //cout << avgX << endl;
    double avgY = accumulate(y.begin(), y.end(), 0.0) / n; //cout << avgY << endl;

    double numerator = 0.0;
    double denominator = 0.0;

    //cout << "x_i    y_i" << endl;

    for(int i=0; i<n; ++i)
    {
        numerator   += (x[i] - avgX) * (y[i] - avgY);
        denominator += (x[i] - avgX) * (x[i] - avgX);

        //cout << x[i] << "\t" << y[i] << endl;
    }    

    //CAUTION: In our case this will never happen because time is increasing always
    if(denominator == 0)
    {
        ;//cout << "Denominator equal to zero" << endl;
        //throw exception("...");
    }

    ret.slope     = ((double) numerator) / denominator;
    ret.intersept = avgY - ret.slope * avgX; 

    //cout << "Slope: "     << ret.slope << endl;
    //cout << "Intercept: " << ret.intersept << endl;

    return ret;
}

//Gives the estimate of x-coord based on the model learn already
double giveEstimate(regVals regLine, double x)
{
    //y = a * x + b
    return (double)(regLine.slope * x + regLine.intersept);
}

/*int main() 
{
    //static const double x[] = {10,30,50,70,90,110};
    vector<double> v1; // = (x, x + sizeof(x) / sizeof(x[0]) );
    v1.push_back(10);
    v1.push_back(30);
    v1.push_back(50);
    v1.push_back(70);
    v1.push_back(90);
    v1.push_back(110);

    vector<double> v2;
    v2.push_back(5);
    v2.push_back(8);
    v2.push_back(6.5);
    v2.push_back(9);
    v2.push_back(12);
    v2.push_back(10);

    for ( vector<double>::iterator it = v1.begin() ; it != v1.end() ; it++)
    {
        cout << *it << endl;
    }

    //Initialize the linear Regression
    regVals ret;
    ret = linearReg(v1,v2);

    cout << "slope = " << ret.slope << endl;
    cout << "intersept = " << ret.intersept << endl;

    //Give the estimate of point with x = 130
    cout << giveEstimate(ret, 130) << endl;

    return 0;
}*/