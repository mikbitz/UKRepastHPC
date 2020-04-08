/*
 *
 * Human.h
 *
 *  Created on: March 4, 2020
 *      Author: Mike Bithell
 * 

 */

#ifndef HUMAN_H_
#define HUMAN_H_
#include "agent.h"
#include "repast_hpc/AgentId.h"
#include "repast_hpc/SharedContext.h"
#include "repast_hpc/SharedDiscreteSpace.h"
#include "repast_hpc/Properties.h"
#include "repast_hpc/Utilities.h"

#include "AgentPackage.h"
#include "EnvironmentCell.h"
#include "randomizer.h"
#include "disease.h"

class MadModel;

class Human: public MadAgent  {

public:

    int _sequencer;//used to co-ordinate ordering of cohort behaviour across threads
    unsigned _FunctionalGroupIndex;   
            
    double _IndividualBodyMass;
    
    unsigned _BirthTimeStep;            
    unsigned _MaturityTimeStep;            
    
    std::string _Realm;      
   
    char _sex;
    
    bool _IsMature;
    
    unsigned _CurrentTimeStep;
    
    static void setParameters(repast::Properties*);

    //Pi!
    static double _Pi ;
    //area units conversion
    static double _CellAreaToHectares;
    
    //temporary store for within timestep changes
    static std::map < std::string, std::map<std::string,double> > _Accounting;

public:

    map<string,disease>_diseases;
    static unsigned _NextID;
    Human* _newH;
    Human(repast::AgentId id): MadAgent(id){_NextID++;_newH=NULL;_sequencer=0;}
    //for copy across threads (needs increaseNextID=false) or restore from file (set increaseNextID to true)
	Human(repast::AgentId id, const AgentPackage& package,bool increaseNextID=false): MadAgent(id){PullThingsOutofPackage(package);_newH=NULL;if (increaseNextID)_NextID++;_sequencer=0;}
    void set(int currentRank, const AgentPackage& package){_id.currentRank(currentRank);PullThingsOutofPackage(package);}
	void setup(unsigned,unsigned,EnvironmentCell*,randomizer*);
    void setPropertiesFromCohortDefinitions(unsigned);
	virtual ~Human() {}

	void step(vector<Human*>&,const unsigned,MadModel*);

    void metabolize();
    void reproduce();
    void interact(vector<Human*>&,MadModel*);
    void moveIt(EnvironmentCell*,MadModel*);
    void mort();
    void markForDeath();
    void applyEcology();
    void setupOffspring( Human* , double , double , double , double , unsigned  );
    void TryToDisperse(double,EnvironmentCell*,MadModel* );
    void TryToDisperse(double,double,EnvironmentCell*,MadModel*);
    vector<double> dDirect(double,double,EnvironmentCell*);
    bool inDistance(MadAgent*, MadAgent*,MadModel *);
    void PushThingsIntoPackage( AgentPackage& );
    void PullThingsOutofPackage( const AgentPackage& );
    void ResetAccounts();
    void infectWith(std::string);
    bool hasDisease(std::string);
    bool recoveredFrom(std::string);
    void updateDiseases();
};
#endif /* HUMAN_H_ */
