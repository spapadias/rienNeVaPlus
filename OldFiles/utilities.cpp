//Title: utilities.cpp
//Description: Useful functions
//Platform(s): gcc (Ubuntu/Linaro 4.6.3-1ubuntu5) 4.6.3
//Author: PAPADIAS Serafeim
//Date:     13/05/2017
//Revision: 13/05/2017

//---------------------------------------------------------//
//--------From Kyriakos motif paper -- EDBT 2017 ----------//
//---------------------------------------------------------//
double toRadians(double degree)
{
	return degree * PI / 180;
}

double haversine(double latitude1, double longitude1, double latitude2, double longitude2)
{
	double R = 6371;
	latitude1 = toRadians(latitude1);
	latitude2 = toRadians(latitude2);
	longitude1 = toRadians(longitude1);
	longitude2 = toRadians(longitude2);
	double deltaLatitude = latitude2 - latitude1;
	double deltaLongitude = longitude2 - longitude1;
	double a = sin(deltaLatitude / 2)*sin(deltaLatitude / 2) + cos(latitude1) * cos(latitude2) * sin(deltaLongitude / 2)*sin(deltaLongitude / 2);
	double c = 2 * asin(sqrt(a));
	return c * R; // return km
}

double dist(double x1, double y1, double x2, double y2)
{
	//return sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2));
	return haversine(x1, y1, x2, y2);
}

double max(double a, double b)
{
	if (a >= b) return a;
	else return b;
}

double min(double a, double b)
{
	if (a <= b) return a;
	else return b;
}

double min3(double f1, double f2, double f3)
{
	double minf = (f1 < f2) ? f1 : f2;
	minf = (minf < f3) ? minf : f3;
	return minf;
}
double dist3(double x1, double y1, double z1, double x2, double y2, double z2)
{
	double R = 6371;
	double c = x1*x2 + y1*y2 + z1*z2;
	c = min(1.0, c);
	c = max(-1.0, c);
	double a = acos(c);
	return a * R;
	//return haversine(x1, y1, x2, y2);
}

void quickSortVector(vector<Point>& candidate, int l, int r)
{
	int i, j;
	double mid;
	Point tmp;
	i = l; j = r; mid = candidate[(i + j) / 2].value;
	while (i <= j)
	{
		while (candidate[i].value < mid) i++;
		while (candidate[j].value > mid) j--;
		if (i <= j)
		{
			tmp = candidate[i]; candidate[i] = candidate[j]; candidate[j] = tmp;
			i++; j--;
		}
	}
	if (i<r) quickSortVector(candidate, i, r);
	if (j>l) quickSortVector(candidate, l, j);
}

void quickSort(double *value, int *x, int *y, int l, int r)
{
	int i, j;
	double mid;
	double tmp;
	int tmpInt;
	i = l; j = r; mid = value[(i + j) / 2];
	while (i <= j)
	{
		while (value[i] < mid) i++;
		while (value[j] > mid) j--;
		if (i <= j)
		{
			tmp = value[i]; value[i] = value[j]; value[j] = tmp;
			tmpInt = x[i]; x[i] = x[j]; x[j] = tmpInt;
			tmpInt = y[i]; y[i] = y[j]; y[j] = tmpInt;
			i++; j--;
		}
	}
	if (i<r) quickSort(value, x, y, i, r);
	if (j>l) quickSort(value, x, y, l, j);
}