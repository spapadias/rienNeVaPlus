//Title: parsingFunctions.cpp
//Description: Functions for parsing the input queries.
//Author: Kostas Patroumpas, Papadias Serafeim
//Platform(s): gcc (Ubuntu/Linaro 4.6.3-1ubuntu5) 4.6.3
//Date:     04/04/2015
//Revision: 10/05/2016

vector<string> split(string s, string sep=":")  
{
    vector<string> resString;
    string::size_type i=0,j=0;

    while( i != string::npos ) 
    {
        i = s.find_first_not_of(sep, j);
        j = s.find_first_of(    sep, i);
        if( (i != string::npos) && (i != j) )
            resString.push_back( s.substr( i, (j-i) ) );
    }

    return resString;
}

std::pair <polygon*, polygon*> constructPolygon(string sGeom)
{
    polygon *p   = new polygon();                           //The polygon itself
    polygon *box = new polygon();                           //Its axis-aligned  bounding box 
    
    double xMin, xMax, yMin, yMax;                          //Coordinates of the axis-aligned bounding box 
    
    //Decode geometry string into an array of coordinates...
    vector<string> coordinates;
    coordinates = split(sGeom, SEPARATOR);                  //Changed the SEPARATOR to comma ","
    
    //Initialize min/max values according to coordinates of the first vertex
    xMin = xMax = atof(coordinates[0].c_str());
    yMin = yMax = atof(coordinates[1].c_str());
        
    //... and finally into a sequence of vertices
    //ASSUMPTION: Even values correspond to x-ordinates, odd values to y-ordinates
    for ( unsigned int i = 0 ; i < coordinates.size() ; i += 2 )
    {
        vertex *v = new vertex( atof(coordinates[i].c_str()), atof(coordinates[i+1].c_str()) );
        p->vertices.push_back(v);
        
        //Update min/max values along each axis with coordinates of the current vertex
        xMin = min(xMin, v->x);
        xMax = max(xMax, v->x);
        yMin = min(yMin, v->y);
        yMax = max(yMax, v->y);
    }

    //Create axis-aligned bounding box for this geometry (counter-clockwise listing of vertices)
    box->vertices.push_back( new vertex(xMin, yMin) );
    box->vertices.push_back( new vertex(xMax, yMax) );
    //box->vertices.push_back(new vertex(xMax, yMin));
    //box->vertices.push_back(new vertex(xMin, yMax));
    
    //Return a pair of geometry objects, consisting of the polygon and its bounding box
    return std::pair<polygon *, polygon *>(p, box);
}

string polygon2string(polygon *p)
{
    std::stringstream sGeom;
    
    //String enclosed by brackets [...] using at most 6 decimal digits for coordinates
    sGeom << "[";
    for (vector<vertex *>::iterator it = p->vertices.begin() ; it != p->vertices.end(); ++it)
        sGeom << setprecision(6) << fixed << (*it)->x << DELIMITER << setprecision(6) << fixed << (*it)->y << DELIMITER;
    
    sGeom << "]";
    return sGeom.str();
}