/*
 *
 * Human.cpp
 *
 *  Created on: March 4, 2020
 *      Author: Mike Bithell
 * 
 */


#include "repast_hpc/RepastProcess.h"
#include "repast_hpc/Random.h"

#include "Groups.h"
#include "Human.h"
#include "Parameters.h"
#include "TimeStep.h"
#include "UtilityFunctions.h"
#include "model.h"


using namespace repast;
    //shared constants - these are defaults - Humans::setParameters will get these from model.props
    
    //Pi!
    double Human::_Pi = acos(-1);
    //cell areas are in sq. km.
    double Human::_CellAreaToHectares=100;
    
   
    
    //needs to be class variable so that Humans can have unique numbers
    unsigned Human::_NextID=0;
    //_Accounting is shared by all cohorts - saves *a lot* of memory
    //only works as these are temporaries used only and completely within each call to step() by a cohort
    //NB do not reset within step() (e.g. by having newly reproduced cohorts call ResetAccounts() in setupOffspring() !)
    std::map < std::string, std::map<std::string,double> > Human::_Accounting;
//------------------------------------------------------------------------------------------------------------
void Human::setParameters(repast::Properties* props){
    //shared constants - static function to read these from parameter file

    
     _Pi                                          = acos(-1.);
     _CellAreaToHectares                          = repast::strToDouble(props->getProperty("HumanParameters.CellAreaToHectares"));


}
//------------------------------------------------------------------------------------------------------------
void Human::ResetAccounts( ) {
    // Initialize delta abundance sorted list with appropriate processes

    _Accounting["abundance"]["mortality"] = 1.0;//NB this applies because of a change from the original model - this value is now a multiplier (reduces possibility of negatives)

    // Initialize delta biomass sorted list with appropriate processes
    _Accounting["biomass"]["metabolism"] = 0.0;
    _Accounting["biomass"]["carnivory"] = 0.0;
    _Accounting["biomass"]["herbivory"] = 0.0;
    _Accounting["biomass"]["reproduction"] = 0.0;

    // Initialize delta reproductive biomass vector with appropriate processes

    _Accounting["reproductivebiomass"]["reproduction"] = 0.0;

    // Initialize organic pool delta vector with appropriate processes
    _Accounting["organicpool"]["herbivory"] = 0.0;
    _Accounting["organicpool"]["carnivory"] = 0.0;
    _Accounting["organicpool"]["mortality"] = 0.0;

    // Initialize respiratory CO2 pool delta vector with appropriate processes
    _Accounting["respiratoryCO2pool"]["metabolism"] = 0.0;
}
//used to create an initial set of cohorts at the start of a run
//------------------------------------------------------------------------------------------------------------
void Human::setup(unsigned functionalGroup,unsigned numHumansThisCell,EnvironmentCell* e,randomizer* r){
    ResetAccounts( );
    _sequencer = 0;
    _FunctionalGroupIndex=functionalGroup;

    _alive = true;
    _moved=false;
    _location={0,0};
    _IsMature=true;

    setPropertiesFromCohortDefinitions(_FunctionalGroupIndex);

    _BirthTimeStep=0;
    _MaturityTimeStep=std::numeric_limits<unsigned>::max( );
   
    _sex='f';
    if (repast::Random::instance()->nextDouble()<0.5)_sex='m';


    _IndividualBodyMass=70000; //70kg for the present

}
//------------------------------------------------------------------------------------------------------------
void Human::setPropertiesFromCohortDefinitions(unsigned functionalGroup){

    _Realm      ="terrestrial";

    //_MinimumMass=CohortDefinitions::Get()->Property(functionalGroup   ,"minimum mass");could read indiv mass from a file
    //_MaximumMass=CohortDefinitions::Get()->Property(functionalGroup   ,"maximum mass");
}
//------------------------------------------------------------------------------------------------------------
//Required by RHPC for cross-core copy - NB "Accounts" do not need to be included as they are instantaneous within a timestep
void Human::PullThingsOutofPackage( const AgentPackage& package ) {

    _sequencer                   = package._contents._sequencer;
    _FunctionalGroupIndex        = package._contents._FunctionalGroupIndex;

    _IndividualBodyMass          = package._contents._IndividualBodyMass;
    _BirthTimeStep               = package._contents._BirthTimeStep;
    _MaturityTimeStep            = package._contents._MaturityTimeStep;

    _alive                       = package._contents._alive;

    setPropertiesFromCohortDefinitions(_FunctionalGroupIndex);
  
    _IsMature=package._contents._IsMature;
    
    _moved=package._contents._moved;
    _location=package._contents._location;
    _destination=package._contents._destination;
    _sex=package._contents._sex;
    _diseases=package._contents._diseases;
}
//------------------------------------------------------------------------------------------------------------
//Required by RHPC for cross-core copy
void Human::PushThingsIntoPackage( AgentPackage& package ) {
    package._contents._sequencer                   =  _sequencer;
    package._contents._FunctionalGroupIndex        =  _FunctionalGroupIndex;

    package._contents._IndividualBodyMass          =  _IndividualBodyMass;
    package._contents._BirthTimeStep               =  _BirthTimeStep;
    package._contents._MaturityTimeStep            =  _MaturityTimeStep;

    package._contents._alive                       =  _alive;
 
   
    package._contents._IsMature= _IsMature;
    
    package._contents._moved=_moved;
    package._contents._location=_location;
    package._contents._destination=_destination;
    package._contents._sex=_sex;
    package._contents._diseases=_diseases;
    
}
//------------------------------------------------------------------------------------------------------------
void Human::setupOffspring( Human* actingHuman, double juvenileBodyMass, double adultBodyMass, double initialBodyMass, double initialAbundance, unsigned birthTimeStep ) {
     
     _newH=NULL;
    _sequencer                   = 0;
    _FunctionalGroupIndex        = actingHuman->_FunctionalGroupIndex;
    _IndividualBodyMass          = initialBodyMass;
    _CurrentTimeStep             = birthTimeStep;
    _BirthTimeStep               = birthTimeStep;
    _MaturityTimeStep            = std::numeric_limits<unsigned>::max( );

    _alive                       = true;

  
    _Realm      =    actingHuman->_Realm;

    _IsMature=false;
       
    _moved=false;
    _location=actingHuman->_location;
    _destination=_location;


}
//------------------------------------------------------------------------------------------------------------
void Human::step(vector<Human*>& others,const unsigned Timestep,MadModel* m) {
    _newH=NULL;//make sure the reproduction pointer has been zeroed out


    _CurrentTimeStep=Timestep;

    ResetAccounts( );
    for (auto A: _Accounting)for (auto B: A.second)if(B.first!="mortality")assert(B.second==0);

    if (m->_interacting )interact(others,m);
    if (m->_metabolism  )metabolize();
    if (m->_reproduction)reproduce();
    if (m->_death       )mort();
    applyEcology();

}
//------------------------------------------------------------------------------------------------------------
void Human::interact(vector<Human*>& others,MadModel* m){
    if (!hasDisease("covid") || recoveredFrom("covid")) return;
    double _DaysInATimeStep=TimeStep::instance()->DaysPerTimeStep();
    
    
    // Loop over potential prey functional groups
    for (auto& agent: others){
        if (inDistance(this,agent,m)){
           for (auto& [name,d]:_diseases){if (d.infectious() && !agent->hasDisease(name) && repast::Random::instance()->nextDouble() < d.infectionProb())agent->infectWith(name);} 
        }
    }
}
//------------------------------------------------------------------------------------------------------------
void Human::infectWith(std::string name){
    _diseases[name].infect();
}
//------------------------------------------------------------------------------------------------------------
bool Human::hasDisease(std::string name){
    return (_diseases.find(name)!=_diseases.end());
}
//------------------------------------------------------------------------------------------------------------
void Human::updateDiseases(){
    for (auto& [name,d]:_diseases)d.update();
}
//------------------------------------------------------------------------------------------------------------
bool Human::recoveredFrom(std::string name){
    return _diseases[name].recovered();
}
//------------------------------------------------------------------------------------------------------------
void Human::markForDeath(){
    if ( _IndividualBodyMass <= 0){ 

      //mark the cohort but don't kill it yet to avoid any problems with movement code in parallel (or there may be disease spread after death)
      _alive=false;
    }
}
//------------------------------------------------------------------------------------------------------------
void Human::moveIt(EnvironmentCell* e,MadModel* m){

        _moved=false;
        if (!_alive)return;
        _destination=_location;
        vector<int> movement={0,0};
       // Calculate the scalar to convert from the time step units used by this implementation of dispersal to the global model time step units
        double DeltaT = TimeStep::instance()->MonthsPerTimeStep();

        double dispersalSpeed=1;// * DeltaT;//move one random cell
        TryToDisperse( dispersalSpeed,e,m );

        //all cohorts need to update their current position
        _location=_destination;
      
}
//------------------------------------------------------------------------------------------------------------
void Human::TryToDisperse(double dispersalSpeed, EnvironmentCell* e,MadModel* m){
    double randomDirection = repast::Random::instance()->nextDouble()* 2 * _Pi;

    // Calculate the u and v components given the dispersal speed
    double uSpeed = dispersalSpeed * cos( randomDirection );
    double vSpeed = dispersalSpeed * sin( randomDirection );
    TryToDisperse(uSpeed, vSpeed,e,m);
 }
  //------------------------------------------------------------------------------------------------------------
void Human::TryToDisperse(double uSpeed, double vSpeed,EnvironmentCell* e, MadModel* m){

      vector<double> uv={0,0};//default to no dispersal
      if (m->_dispersalSelection=="direct"       ) uv=dDirect(uSpeed,vSpeed,e);
      
      double signu=uv[0];
      double signv=uv[1]; 

      double x=_destination[0],y=_destination[1];//need to accumulate this over multiple steps or some movements might get missed/not wrap properly

     //only move if realm matches and we haven't exceeded upper and lower latitude bounds (currently no movement across the pole)
     //NB last cell lies in the range _maxX<= x <_mMax+1
     int yw=floor(y+signv);

     if (!m->_noLongitudeWrap){
        if (x+signu < m->_minX)    {x = x + (m->_maxX - m->_minX + 1);}
        if (x+signu >= m->_maxX + 1){x = x - (m->_maxX - m->_minX + 1);}
     }
     int xw=floor(x+signu);
     if (yw >= m->_minY && yw <= m->_maxY){
         if (xw >= m->_minX && xw <= m->_maxX) {
             EnvironmentCell* E=m->_Env[xw][yw];
             if (E->_Realm==_Realm || _Realm=="all"){// no movement if wrong realm at destination
               _destination[0]=x+signu;_destination[1]=y+signv;
               if (xw!=floor(_location[0]) || yw!=floor(_location[1]))_moved=true;//special treatment needed if we have changed cell
           }          
         }
     } 
   }
//------------------------------------------------------------------------------------------------------------
 vector<double> Human::dDirect(double uSpeed, double vSpeed,EnvironmentCell* e){
    //dispersal where cohorts can take on fractional cell co-ordinates
    //this is *required* when cross-cell interaction becomes important
    // Calculate the fraction of the grid cell in the u direction 
    double ufrac = ( uSpeed / e->Width() );

    // Calculate the fraction of the grid cell in the v direction
    double vfrac = ( vSpeed / e->Height() );
    return vector<double>({ufrac,vfrac});
 }
//------------------------------------------------------------------------------------------------------------
bool Human::inDistance(MadAgent* a1, MadAgent* a2,MadModel* m){
    //return 0.5;//kind of a whole-cell default...
    //simply - the number of cell widths, Manhattan style
    double wrappedDistX=abs(a1->_location[0]-a2->_location[0]);
    if (!m->_noLongitudeWrap)wrappedDistX=min(wrappedDistX,(m->_maxX - m->_minX+1)-wrappedDistX);
    return max(wrappedDistX,abs( a1->_location[1] - a2->_location[1]))<1.;
    //alternative using proper spherical distance.
    double dg=Parameters::instance()->GetGridCellSize(); //in degrees
    //get lon lat at the current location allowing for fractions of a cell
    double lon1 = Parameters::instance()->GetLongitudeAtIndex(int(a1->_location[0])) + (a1->_location[0]-int(a1->_location[0]))*dg;
    double lat1 = Parameters::instance()->GetLatitudeAtIndex (int(a1->_location[1])) + (a1->_location[1]-int(a1->_location[1]))*dg;
    double lon2 = Parameters::instance()->GetLongitudeAtIndex(int(a2->_location[0])) + (a2->_location[0]-int(a2->_location[0]))*dg;
    double lat2 = Parameters::instance()->GetLatitudeAtIndex (int(a2->_location[1])) + (a2->_location[1]-int(a2->_location[1]))*dg;
    UtilityFunctions u;
    return u.HaversineDistanceInDegrees(lon1,lat1,lon2,lat2)<=0.25*dg;//one quarter cell at the equator, more at the poles, at least in longitude
    //Interaction distance - how should this be set, given that the distance of a longitude cell varies with latitude?
    //one choice would be the max lat cell size - but then this varies with the latitudinal range (and could be zero!).
    //Alternatively maybe one could use the latitude mid way inbetween equator and max. lat - the danger is that near-cell-edge cohorts get longitude truncated interaction sets
    //nearer to the pole.
    //or one could use the average of lat 1 and lat 2, but then polar cohorts really do get less interaction (at least when near cell edges) - although this is already true given
    //lat-long based grid cells (equatorial cohorts have a greater effective interaction range!).
    //so this version more or less mimics the within-cell-only case
    //preferably the interaction would span the same range at all latitudes, but if this is done near poles then more than one-cell-range may be needed for interaction
    //implying more cross-process copying. This is really true in any case - a switch to equal area cells of some kind is needed really - e.g. buckyball or hex.
    //However. given that we use grid cell size/4 it should be OK to use the full equatorial distance right out to latitude 60 (where cos lat=0.5), 
    //beyond which there will be some artefacts if cohorts don't eat all biomass that they can see...

}
//------------------------------------------------------------------------------------------------------------
void Human::metabolize(){

}
//------------------------------------------------------------------------------------------------------------
void Human::reproduce(){
    _newH=NULL;//memory allocation for new cohorts is dealt with via the model list of agents 

}

//------------------------------------------------------------------------------------------------------------
void Human::mort(){

}
//------------------------------------------------------------------------------------------------------------
void Human::applyEcology(){



}






