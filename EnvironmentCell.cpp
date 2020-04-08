/*
 *
 * EnvironmentCell.cpp
 *
 *  Created on: May 7, 2018
 *      Author: Mike Bithell
 * 
 */
#include "EnvironmentCell.h"
#include "DataLayerSet.h"
#include "TimeStep.h"
#include "Parameters.h"
#include "UtilityFunctions.h"


EnvironmentCell::EnvironmentCell() {}

//------------------------------------------------------------------------------
EnvironmentCell::EnvironmentCell(int x,int y):_x(x),_y(y){
    UtilityFunctions Utility;
    _Area=   Utility.CalculateGridCellArea(Latitude(),cellSize());//in sqkm
    _Width=  Utility.CalculateLengthOfDegreeLatitude( Latitude()) *cellSize();
    _Height= Utility.CalculateLengthOfDegreeLongitude( Latitude())*cellSize();


    SetRealm( );

}
//------------------------------------------------------------------------------

void EnvironmentCell::SetRealm( ) {
    _Realm="none";
    
    /*if( DataLayerSet::Data( )->GetDataAtLonLatFor( "Realm", Longitude(),  Latitude() ) <= 1.5 ) {
          _Realm="marine";
    } else if( DataLayerSet::Data( )->GetDataAtLonLatFor( "Realm", Longitude(),  Latitude() ) > 1.5) {
          _Realm="terrestrial";
    }*/
    if( DataLayerSet::Data( )->GetDataAtLonLatFor( "Population", Longitude(),  Latitude() ) <= 0 ) {
          _Realm="marine";
    } else {
          _Realm="terrestrial";
    }

}
//----------------------------------------------------------------------------------------------

double EnvironmentCell::Width(){return _Width;}
double EnvironmentCell::Height(){return _Height;}
double EnvironmentCell::Area(){return  _Area;}

//------------------------------------------------------------------------------
double EnvironmentCell::Temperature(){

 double d = 0;
 
 if( _Realm=="marine" ) {
    d = GetVariableFromDatasetNamed("MarineTemp");
 } else if( _Realm=="terrestrial" ) {
    d = GetVariableFromDatasetNamed( "TerrestrialTemp");
 }

 return d;
}
//------------------------------------------------------------------------------
double EnvironmentCell::Population(){

  double d = GetVariableFromDatasetNamed("Population");

  return d;
    
}
//------------------------------------------------------------------------------
double EnvironmentCell::Precipitation(){

 double d = Constants::cMissingValue;//currently no marine precip
 if( _Realm=="terrestrial" ) {
   d = GetVariableFromDatasetNamed("TerrestrialPre");
  }
  return d;
    
}
//------------------------------------------------------------------------------
double EnvironmentCell::GetVariableFromDatasetNamed(std:: string s){

    double d = Constants::cMissingValue;
    d = DataLayerSet::Data( )->GetDataAtLonLatFor( s, Longitude(),  Latitude() );
    if( d == Constants::cMissingValue ) {
      std::cout << "Warning EnvironmentCell::GetVariableFromDatasetNamed- missing values in "<<s<<" field!! "<<Longitude()<<" E "<<Latitude()<<" N "<< std::endl;
    }
    return d;
}
//------------------------------------------------------------------------------
double EnvironmentCell::Latitude(){
    return Parameters::instance()->GetLatitudeAtIndex(_y);
}
//------------------------------------------------------------------------------
double EnvironmentCell::Longitude(){
    return Parameters::instance()->GetLongitudeAtIndex(_x);
}
//------------------------------------------------------------------------------
double EnvironmentCell::cellSize(){
    return Parameters::instance()->GetGridCellSize();
}



