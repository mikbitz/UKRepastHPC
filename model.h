/* Model.h */

#ifndef MODEL
#define MODEL

#include <boost/mpi.hpp>
#include "repast_hpc/Schedule.h"
#include "repast_hpc/Properties.h"
#include "repast_hpc/AgentRequest.h"
#include "repast_hpc/TDataSource.h"
#include "repast_hpc/SVDataSet.h"
#include "repast_hpc/SVDataSetBuilder.h"
#include "repast_hpc/SharedContext.h"
#include "repast_hpc/SharedDiscreteSpace.h"
#include "repast_hpc/GridComponents.h"
#include "Environment.h"
#include "EnvironmentCell.h"
#include "Human.h"
#include "agent.h"


class MadModel;
class AgentPackage;
//This is an alias for the background discrete grid
typedef repast::SharedDiscreteSpace<MadAgent, repast::WrapAroundBorders, repast::SimpleAdder<MadAgent> > wrappedSpaceType;


//------------------------------------------------------------------------------------------
//Structs for moving agents across threads
//------------------------------------------------------------------------------------------
class MadAgentPackageProvider {
	
private:
    repast::SharedContext<MadAgent>* agents;
	
public:
	
   MadAgentPackageProvider(repast::SharedContext<MadAgent>* agentPtr);

    void providePackage(MadAgent * agent, std::vector<AgentPackage>& out);

    void provideContent(repast::AgentRequest req, std::vector<AgentPackage>& out);
	
};
//------------------------------------------------------------------------------------------
class MadAgentPackageReceiver {
	
private:
    repast::SharedContext<MadAgent>* agents;
	
public:
	
    MadAgentPackageReceiver(repast::SharedContext<MadAgent>* agentPtr);
	
    MadAgent * createAgent(AgentPackage package);
	
    void updateAgent(AgentPackage package);
	
};
//------------------------------------------------------------------------------------------
//Actual model class
//------------------------------------------------------------------------------------------
class MadModel{
    
    unsigned _startingStep;
	int _stopAt;
    
    bool _verbose;
    unsigned _randomSeed;
	repast::Properties* _props;
	repast::SharedContext<MadAgent> _context;
    std::vector<repast::DataSet*> dataSets;
	
	MadAgentPackageProvider* provider;
	MadAgentPackageReceiver* receiver;

    wrappedSpaceType* discreteSpace;
    int _totalSusceptible;
    int _totalInfected;
    int _totalRecovered;
    int _totalDied;
    int _totalPopulation;
    
    //vector<int> _FinalCohortBreakdown;
      
    map< string,vector<double> > outputMaps;
    map<string,string> outputUnits;
    vector<string> outputNames;
    vector<int> _cellSelector;
    bool _crossCell;
    unsigned _restartInterval;
    unsigned _restartStep;
    std::string _restartDirectory;
    std::string _archiveFormat;
    
    std::string _filePrefix, _filePostfix;
    void dataSetClose();
    void addDataSet(repast::DataSet*) ;
    void setupNcOutput();
    void netcdfOutput( unsigned step );
    void setNcGridFile(std::string,std::string );
    void writeNcGridFile(unsigned,vector<double>&,std::string);
    std::vector<AgentPackage>_packages;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & _packages;
    }

public:
    bool _interacting  ;
    bool _metabolism   ;
    bool _reproduction ;
    bool _death        ;
    bool _dispersal    ;
    bool _mergers      ;
    bool _output       ;
    
    int _minX,_minY,_maxX,_maxY,_dimX,_dimY;
    int _xlo,_xhi,_ylo,_yhi;
    int _noLongitudeWrap; //1 if domain does *not* span the global longitude range
    string _dispersalSelection;
    Environment _Env;
    
	MadModel(repast::Properties& ,  boost::mpi::communicator* comm);
	~MadModel();
	void init(unsigned);
	void moveAgents();
	void step();
    void sync();
	void initSchedule(unsigned, repast::ScheduleRunner& runner);
	void recordResults();
    wrappedSpaceType* space(){return discreteSpace;}

    static int  _humanType;
    //outputs
    void read_restart(unsigned);
    void write_restart();
    void setupOutputs();
    void tests();
    void setupHumanTestValues(Human*);
    void checkHumanTestValues(Human*);
    int PopCount() const {
		return _totalPopulation;
	}
    int SusCount() const {
		return _totalSusceptible;
	}
	double InfCount() const {
		return _totalInfected;
	}
	double RecCount() const {
		return _totalRecovered;
	}
	double DeathCount() const {
		return _totalDied;
	}


};
#endif
