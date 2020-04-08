/*
 *
 * EnvironmentCell.h
 *
 *  Created on: May 7, 2018
 *      Author: Mike Bithell
 * 
 *  Derived from Origianl C# code by
 * Emergent Global Patterns of Ecosystem Structure and Function from a Mechanistic General Ecosystem Model , Michael B. J. Harfoot , Tim Newbold, Derek P. Tittensor,  Stephen Emmott, Jon Hutton, Vassily Lyutsarev, Matthew J. Smith, Jörn P. W. Scharlemann, Drew W. Purves PLOS Biology 12(4) 2014 e1001841 https://doi.org/10.1371/journal.pbio.1001841
 */
#ifndef ENVIRONMENTCELL
#define ENVIRONMENTCELL
#include "Constants.h"
#include <vector>
class EnvironmentCell {
public:
	EnvironmentCell();
    EnvironmentCell(int,int);
    void SetRealm( );


    std::string _Realm;


	virtual ~EnvironmentCell() {}
	
	void merge(int&);
    double Temperature();
    double Width();
    double Height();
    double Area();
    double Precipitation();
    double Population();
    double cellSize();
    double Latitude();
    double Longitude();
    unsigned _x,_y;

private:
    double GetVariableFromDatasetNamed(std:: string s);

    double _Width;
    double _Height;
    double _Area;


    

};
#endif


