#include "disease.h"
#include "TimeStep.h"
void disease::infect(){_infected=true;_timer=0;}
void disease::recover(){_recovered=true;_infectious=false;}
void disease::becomeInfectious(){_infectious=true;}
double disease::infectionProb(){return _infectionProb;}
bool disease::infected(){return _infected;}
bool disease::recovered(){return _recovered;}
bool disease::infectious(){return _infectious;}
disease::disease(){
           bool _infected=false;
           bool _recovered=false;
           bool _infectious=false;
           unsigned _timer=0;
           _infectionProb=0.1;
}
void disease::update(){
    double _DaysInATimeStep=TimeStep::instance()->DaysPerTimeStep();
    _timer++;
    if (float(_timer)*_DaysInATimeStep> 2) becomeInfectious();
    if (float(_timer)*_DaysInATimeStep> 10) recover();
    //if (someconditione) die();
}
