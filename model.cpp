/* Model.cpp */

#include <stdio.h>
#include <vector>
#include <sstream>
#include <fstream>
#include <boost/filesystem.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/archive_exception.hpp>
#include <boost/mpi.hpp>
#include "repast_hpc/AgentId.h"
#include "repast_hpc/RepastProcess.h"
#include "repast_hpc/Utilities.h"
#include "repast_hpc/Properties.h"
#include "repast_hpc/SharedBaseGrid.h"                      // VN2D broken without this!
#include "repast_hpc/VN2DGridQuery.h"
#include "repast_hpc/Moore2DGridQuery.h"
#include "repast_hpc/initialize_random.h"
#include "repast_hpc/Point.h"
#include "repast_hpc/AgentRequest.h"
#include "repast_hpc/SVDataSet.h"
#include "repast_hpc/SVDataSetBuilder.h"
#include "repast_hpc/Random.h"

#ifndef _WIN32
#include "repast_hpc/NCDataSetBuilder.h"
#endif

#include "model.h"
#include "agent.h"
#include "EnvironmentCell.h"
#include "Parameters.h"
#include "Groups.h"
#include "Human.h"
#include "FileReader.h"
#include "PopStats.h"
#include "TimeStep.h"
#include "Constants.h"
#include "randomizer.h"
#include "RandomRepast.h"
#include "AgentPackage.h"
#include "UtilityFunctions.h"

#include <netcdf>

//repast shuffleList only works on pointers
//this version is for vectors of int

void shuffleList(std::vector<int>& elementList){

  if(elementList.size() <= 1) return;
  repast::IntUniformGenerator rnd = repast::Random::instance()->createUniIntGenerator(0, elementList.size() - 1);
  int swap;
  for(size_t i = 0, sz = elementList.size(); i < sz; i++){
    int other = rnd.next();
    swap = elementList[i];
    elementList[i] = elementList[other];
    elementList[other] = swap;
  }
}

using namespace std;
using namespace repast;
//arbitrary numbers to distiguish the agents types

int MadModel::_humanType=2;
//------------------------------------------------------------------------------------------------------------
//Constructor and destructor
//------------------------------------------------------------------------------------------------------------
MadModel::MadModel(repast::Properties& props,  boost::mpi::communicator* comm): _context(comm){
    //switch on all model aspects - these might need to be switched off for test purposes.
    _interacting  = props.getProperty("simulation.IncludeInteraction") !="false";
    _metabolism   = props.getProperty("simulation.IncludeMetabolism")  !="false";
    _reproduction = props.getProperty("simulation.IncludeReproduction")!="false";
    _death        = props.getProperty("simulation.IncludeDeath")       !="false";
    _dispersal    = props.getProperty("simulation.IncludeDispersal")   !="false";
    _mergers      = props.getProperty("simulation.IncludeMergers")     !="false";
    _output       = props.getProperty("simulation.IncludeOutput")      !="false";
    _verbose      = props.getProperty("verbose")=="true";

    //-----------------
    //Pull in parameter from the model.props file
    _props = &props;
    //Number of timesteps
	_stopAt = repast::strToInt(_props->getProperty("stop.at"));
    
    std::string rstrt=props.getProperty("simulation.RestartEvery");
    if (rstrt!="")_restartInterval=repast::strToInt(rstrt); else _restartInterval=0;

    rstrt=props.getProperty("simulation.RestartStep");
    if (rstrt!="")_restartStep=repast::strToInt(rstrt); else _restartStep=0;
    
    rstrt=props.getProperty("simulation.RestartDirectory");
    if (rstrt!="")_restartDirectory=rstrt+"/"; else _restartDirectory="./";

    _archiveFormat ="binary";
    if (props.getProperty("simulation.RestartFormat")=="text")_archiveFormat="text";

    //extent of buffer zones in grid units - this many grid cells are shared at the boundary between cores
    int gridBuffer = repast::strToInt(_props->getProperty("grid.buffer"));
    //Grid extent
    _minX=repast::strToInt(_props->getProperty("min.x"));
    _minY=repast::strToInt(_props->getProperty("min.y"));
    _maxX=repast::strToInt(_props->getProperty("max.x"));
    _maxY=repast::strToInt(_props->getProperty("max.y"));
    _dimX=repast::strToInt(_props->getProperty("proc.per.x"));
    _dimY=repast::strToInt(_props->getProperty("proc.per.y"));
    _noLongitudeWrap=repast::strToInt(_props->getProperty("noLongitudeWrap"));
    _dispersalSelection="direct";
    _randomSeed=repast::strToUInt(_props->getProperty("global.random.seed"));
    //-----------------
    //naming convention for output files
    if(repast::RepastProcess::instance()->rank() == 0 && _output){
       _filePrefix=                  _props->getProperty("experiment.output.directory")+
                      "/experiment."+_props->getProperty("experiment.name");
       if (!boost::filesystem::exists(_filePrefix))boost::filesystem::create_directories(_filePrefix);
       std::string runNumber= _props->getProperty("run.number");
       std::string m00="/run_";
       if (runNumber!=""){
           m00=m00+runNumber;
       }else{
         //auto-increment run number if run.number is not set
         int i=0;
         m00="/run_000",runNumber="000";
         std::string pfx="00";
         while(boost::filesystem::exists(_filePrefix+m00)){    // Find a new directory name
           i++;
           std::stringstream ss;
           if (i>9 ) pfx="0";
           if (i>99) pfx="";
           ss<<pfx<<i;
           runNumber=ss.str();
           m00="/run_"+runNumber;
         }
       }
       if (!boost::filesystem::exists(_filePrefix+m00))boost::filesystem::create_directories(_filePrefix+m00);
       _props->putProperty ("run.number",runNumber);
       _filePrefix= _filePrefix+m00+"/";
       _filePostfix="";
       cout<<"Outputfiles will be named "<<_filePrefix<<"<Data Name>"<<_filePostfix<<".<filenameExtension>"<<endl;
    }
    //slightly tricky to get the file prefix to all other threads (needed, for example, for restart names)
    //only thread 0 can know the name since it needs to create a new name based on existing directory names on disk 
    int prefix_size = _filePrefix.size();
    MPI_Bcast(&prefix_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (repast::RepastProcess::instance()->rank() != 0)_filePrefix.resize(prefix_size);
    MPI_Bcast(const_cast<char*>(_filePrefix.data()), prefix_size, MPI_CHAR, 0, MPI_COMM_WORLD);
    //-----------------
	//create the model grid
    repast::Point<double> origin(_minX,_minY);
    repast::Point<double> extent(_maxX-_minX+1, _maxY-_minY+1);
    
    repast::GridDimensions gd(origin, extent);
    
    std::vector<int> processDims;
    processDims.push_back(_dimX);
    processDims.push_back(_dimY);
  
    //because RHPC uses templates for grid wrapping, if you want any wrapping at all , you have to 
    //use a space that is wrapped in both x and y - to use a non-wrapped space implies templating all occurrences of
    //"model" - far too much like hard work; either that or we need a wrapper class for spaces, or different models (with classes per space)
    discreteSpace = new wrappedSpaceType("AgentDiscreteSpace", gd, processDims, gridBuffer, comm);
	
    //The agent container is a context. Add the grid to it.
   	_context.addProjection(discreteSpace);
    
   //-----------------
    //Set up the cross-thread data transferers
	provider = new MadAgentPackageProvider(&_context);
	receiver = new MadAgentPackageReceiver(&_context);
    //---------------
    //variables to hold totals across all cells
    _totalSusceptible=0;
    _totalInfected=0;
    _totalRecovered=0;
    _totalDied=0;
    _totalPopulation=0;
    //----------------
    //local grid extents on this thread
    _xlo=        discreteSpace->dimensions().origin().getX() ;
    _xhi= _xlo + discreteSpace->dimensions().extents().getX();
    _ylo=        discreteSpace->dimensions().origin().getY() ;
    _yhi= _ylo + discreteSpace->dimensions().extents().getY();

}
//------------------------------------------------------------------------------------------------------------
MadModel::~MadModel(){

	delete provider;
	delete receiver;
    for (size_t i = 0; i < dataSets.size(); ++i) {
		delete dataSets[i];
	}

}
//------------------------------------------------------------------------------------------------------------
//Model Initialisation
//------------------------------------------------------------------------------------------------------------
void MadModel::initSchedule(unsigned startingStep,repast::ScheduleRunner& runner){
	runner.scheduleEvent(startingStep, 1, repast::Schedule::FunctorPtr(new repast::MethodFunctor<MadModel> (this, &MadModel::step)));
	runner.scheduleStop(_stopAt);
    runner.scheduleEndEvent(Schedule::FunctorPtr(new MethodFunctor<MadModel> (this, &MadModel::dataSetClose)));

}
//------------------------------------------------------------------------------------------------------------
void MadModel::init(unsigned startingStep){
    _startingStep=startingStep;//starts at 1 to correspond to timestep 0 - confusing...needs a fix
    _totalSusceptible=0;
    _totalInfected=0;
    _totalRecovered=0;
    _totalDied=0;
    _totalPopulation=0;

    //time the initialisation
    repast::Timer initTimer;
    initTimer.start();
    //rank (i.e. the number of this thread) will be needed to make agents unique *between* threads
    int rank = repast::RepastProcess::instance()->rank();

    _crossCell=(_props->getProperty("simulation.CrossCellInteraction")=="true");
    //now set up the environmental cells - note at present this uses the full grid, not just local to this thread
    //so that off-thread environment can be easily queried. Currently some duplication here, but it is not a huge amount of data.
    _Env=Environment(_minX,_maxX,_minY,_maxY);

    //set up the static (i.e. shared) parameters for the Humans
    Human::setParameters(_props);
    //get the definitions of stocks and cohorts

    CohortDefinitions::Initialise(_props->getProperty("input.DataDirectory")+"/"+_props->getProperty("input.CohortDefinitionsFileName"));
  
    unsigned numCohortGroups=CohortDefinitions::Get()->size();
    unsigned humanCount = strToInt(_props->getProperty("human.count"));
    //these will be output maps
    outputNames.push_back("totalSusceptible");
    outputNames.push_back("totalInfected");
    outputNames.push_back("totalRecovered");
    outputNames.push_back("totalDeaths");
    //Values only reduced on thread 0
    if(repast::RepastProcess::instance()->rank() == 0){
      outputUnits["totalSusceptible"]   ="number/sq. km.";
      outputUnits["totalInfected"]      ="number/sq. km.";
      outputUnits["totalRecovered"]     ="number/sq. km.";
      outputUnits["totalDeaths"]        ="number/sq. km.";
      for (auto name: outputNames) outputMaps[name]    =  vector<double> ( (_maxX-_minX+1) * (_maxY-_minY+1),0.0 );

    }
    
    
    randomizer* random=new RandomRepast;
    //seed is set from model.props file - see main.cpp

 
    //explicitly use the local bounds of the grid on this thread to create countOfAgents per cell.
    //Not doing this can lead to problems with agents in distant cells not within the local thread neighbourhood
    //see SharedBaseGrid.h moveTo method
    //although latest updates I have made to RHPC should have fixed this...
    unsigned totalCohorts=0,totalStocks=0;

    unsigned cNum=0,sNum=0;

    unsigned totalStocksThisCell=0;

  
    int s=0,q=0;
    unsigned hF=10000;//functionalgroup ID for humans
    if (_restartStep==0){
        for (int x = _xlo; x < _xhi; x++){
            for (int y = _ylo; y < _yhi; y++){
                EnvironmentCell* E=_Env[x][y];
                double lon=E->Longitude();
                double lat=E->Latitude();
                
                if (E->_Realm=="terrestrial"){
                    humanCount=E->Population();
                    //cout<<lon<<" "<<lat<<" h "<<humanCount<<": ";
                    for (unsigned j=0;j<humanCount;j++){
                        //make sure the agentId is unique on this thread!!
                        // values are int id, int startProc, int agentType, 
                        repast::AgentId id(Human::_NextID, rank, _humanType);
                        //agent also needs id of its current thread
                        id.currentRank(rank);
                        Human* h = new Human(id);
                        h->setup(hF,humanCount, E,random);
                        
    
                        if (q==0 && lon<=-0.29 && lon>-0.31 && lat<=51.51 && lat > 51.49){q++;h->_diseases["covid"].infect();cout<<lon<<" "<<lat<<endl;}

                        _context.addAgent(h);
                        //to get movement right agent needs its own copy of location
                        double xr=0,yr=0;
                        if (_dispersalSelection=="direct"){
                            //allow cohorts to be at locations other than cell centres initially - note using fractional cell co-ordinates
                            xr=(1-2*random->GetUniform())*0.5;
                            yr=(1-2*random->GetUniform())*0.5;
                            if (y+yr <  _minX)    {yr = -yr;}
                            if (y+yr >= _maxX + 1){yr = -yr;}
                            if (!_noLongitudeWrap){
                                if (x+xr <  _minX)    {xr = xr + (_maxX - _minX + 1);}
                                if (x+xr >= _maxX + 1){xr = xr - (_maxX - _minX + 1);}
                            } else {
                                if (x+xr <  _minX)    {xr = -xr;}
                                if (x+xr >= _maxX + 1){xr = -xr;}
                            }
                        }
                        assert(y+yr >= _minY);
                        repast::Point<int> initialLocation(x+xr,y+yr);
                        discreteSpace->moveTo(id, initialLocation);
                        //to get movement right agent needs its own copy of location
                        h->setLocation(x+xr,y+yr);
                        //code needed here for output totals - S/L/I/R/D etc.
                        if (h->hasDisease("covid")){
                            if (h->_alive){
                                if (!h->recoveredFrom("covid"))_totalInfected++;
                                if ( h->recoveredFrom("covid"))_totalRecovered++;
                                _totalPopulation++; 
                            }else{
                                _totalDied++;
                            }
                            
                        }else{
                            _totalSusceptible++;
                            _totalPopulation++;
                        }
                    }
                }
                
                
            }
        }
    }
                     
    if (_restartStep>0){
        read_restart(_restartStep);
        vector<MadAgent*>agents;
        _context.selectAgents(repast::SharedContext<MadAgent>::LOCAL,agents);
        for (auto a:agents){
            if (a->getId().agentType()==_humanType){
                Human* h=(Human *)a;
               //code needed here for output totals - S/L/I/R/D etc.
                if (h->hasDisease("covid")){
                    if (h->_alive){
                        if (!h->recoveredFrom("covid"))_totalInfected++;
                        if ( h->recoveredFrom("covid"))_totalRecovered++;
                        _totalPopulation++; 
                    }else{
                        _totalDied++;
                    }
                    
                }else{
                    _totalSusceptible++;
                    _totalPopulation++;
                }
            }
        }
    }
    if(_verbose)cout<<"rank "<<rank<<" total Infected "<<_totalInfected<<endl;
    if (_output){
     setupOutputs();
     if (rank==0)setupNcOutput();
    }
    long double t = initTimer.stop();
	std::stringstream ss;
	ss << t;
	_props->putProperty("init.time", ss.str());
    sync();

}
//------------------------------------------------------------------------------------------------------------
//Run the model
//------------------------------------------------------------------------------------------------------------
void MadModel::step(){
    _totalSusceptible=0;
    _totalInfected=0;
    _totalRecovered=0;
    _totalDied=0;
    _totalPopulation=0;
    int rank=repast::RepastProcess::instance()->rank();
    unsigned CurrentTimeStep=RepastProcess :: instance ()->getScheduleRunner ().currentTick () - 1;
    //needed to advance the environmental datalayers to the current timestep
   
    _Env.update(CurrentTimeStep);
    
	if(rank == 0) std::cout << " TICK " << CurrentTimeStep << std::endl;

    //vectors length CohortDefinitions::Get()->size() initialized to 0
    vector<int> cohortBreakdown(CohortDefinitions::Get()->size(),0);
    //spatial distributions - for efficiency these should just the the local part, but easier initially to just keep the whole thing on each thread
    //with zeros for non-local, and then reduce onto thread 0
    map< string,vector<double> > localMaps;
    for (auto& name: outputNames) localMaps[name]    =  vector<double> ( (_maxX-_minX+1) * (_maxY-_minY+1),0.0 );

    int buffer=repast::strToInt(_props->getProperty("grid.buffer"));
    

    int range=0;
    //cohorts need an interaction range in the case of interaction to say how many cells should be included (not the range over which cohorts can interact)
    //synchronize the data - a range of at least 1 cell though I think... (although if radius =3 half cells and one central agent 
    //can see all the way to the edge of the Moore neighbourhood). There's a problem here with latitudinal changes in range as expressed in fractions of a cell.
    //so interaction-range might be 0.5 at the equator, but will be 1 by 60 degrees latitude - if the domain stopped there, this would make buffer = 1, range=1, with agents able to see 0.5, but
    //buffer=2 is needed to get all asynchronous interactions to work properly
    //timestep shoudl be short enough to make the interaction range work...to avoid artefacts it seems interaction-range shoudl be not much above 0.25 of a cell.
    if (_crossCell) range=1;
    
    //loop over agents on this thread - using synchronous updates so that update order of agents is not important for disease transfer.
    //For cell interaction, we need to make sure copies are updated appropriately. Since the infection process is one-way, 
    //remote humans in the buffer zone on this thread can update local ones provided their disease data is up-to-date.
    //Do this over cells as the Moore query is expensive
    std::vector<MadAgent*> agents;
bool _cellBased=true;//cell based may be al lot faster if there are significant numbers of humans per cell - fewer Moore lookups
    if (_cellBased){
    for(int x = _xlo - range; x < _xhi + range; x++){
        for(int y = _ylo - range; y < _yhi + range; y++){
            if (x>= _minX && x<=_maxX && y>= _minY && y<=_maxY){
                repast::Point<int> location(x,y);
                
                //query four neighbouring cells, distance 0 (i.e. just the centre cell) - "true" keeps the centre cell.
                repast::VN2DGridQuery<MadAgent> VN2DQuery(space());
                VN2DQuery.query(location, 0, true, agents);
                std::vector<MadAgent*> thingsToInteract;
                
                //things that can be eaten by the agents above - they can be in any of the 
                //eight neighbouring cells (range =1 - requires grid.buffer=1) plus the current cell, or just the current cell (range=0, grid,buffer=0)
                repast::Moore2DGridQuery<MadAgent> Moore2DQuery(space());
                Moore2DQuery.query(location, range, true, thingsToInteract);
                
                
                std::vector<Human*> humansToInteract;
                
                for (auto a:thingsToInteract){
                    if ( a->_alive){//agents must be living but need only be local!
                        //separate out humans
                        if (a->getId().agentType()==_humanType && a->getId().currentRank()==rank) {Human* h=(Human*)a;if (!h->hasDisease("covid"))humansToInteract.push_back(h);}
                    }
                }

                for (auto& a: agents)((Human *)a)->step(humansToInteract, CurrentTimeStep,this);
                agents.clear();
                
            }
        }
    }
}else{
    _context.selectAgents(agents);
     for (auto& a:agents){ 
         //find the right location for this human - could be anywhere across the thread
        std::vector<int> agentLoc;
        discreteSpace->getLocation(a->getId(), agentLoc);
        repast::Point<int> location(agentLoc);
        Human * human=(Human *)a;
        
        std::vector<MadAgent*> thingsToInteract;
        
        //things that can be eaten by the agents above - they can be in any of the 
        //eight neighbouring cells (range =1 - requires grid.buffer=1) plus the current cell, or just the current cell (range=0, grid,buffer=0)
        repast::Moore2DGridQuery<MadAgent> Moore2DQuery(space());
        Moore2DQuery.query(location, range, true, thingsToInteract);
        
        
        std::vector<Human*> humansToInteract;
        
        for (auto a:thingsToInteract){
            if ( a->_alive){//agents must be living but might not be local!
                //separate out humans
                if (a->getId().agentType()==_humanType ) {Human* h=(Human*)a;if (!h->hasDisease("covid"))humansToInteract.push_back(h);}
            }
        }    

        human->step(humansToInteract, CurrentTimeStep,this);
        
    }
}

    //now diseases can be updated, including newly infected agents. Only local need be updated.
    agents.clear();
    _context.selectAgents(repast::SharedContext<MadAgent>::LOCAL,agents);
    for (auto& a:agents){
        std::vector<int> agentLoc;
        discreteSpace->getLocation(a->getId(), agentLoc);
        int cellIndex=agentLoc[0]-_minX+(_maxX-_minX+1)*(agentLoc[1]-_minY);
        //double area=_Env[agentLoc[0]][agentLoc[1]]->Area();
        Human * h=(Human *)a;
        //advance disease states
        h->updateDiseases();

        //acumulate other totals and spatial maps
        if (h->hasDisease("covid")){
            if (h->_alive){
                if (!h->recoveredFrom("covid")){_totalInfected++;localMaps["totalInfected"][cellIndex]+= 1.;}
                if ( h->recoveredFrom("covid")){_totalRecovered++;localMaps["totalRecovered"][cellIndex]+= 1.;}
                _totalPopulation++; 
            }else{
                _totalDied++;localMaps["totalDeaths"][cellIndex]+= 1.;
            }
            
        }else{
            _totalSusceptible++;localMaps["totalSusceptible"][cellIndex]+= 1.;
            _totalPopulation++;
        }
    }
    //updates of offspring and mergers/death happen after all cells have updated
    //need to keep this separate if there is cross-cell interaction
 
    for(int y = _ylo; y < _yhi; y++){  
     for(int x = _xlo; x < _xhi; x++){
      int cellIndex=x-_minX+(_maxX-_minX+1)*(y-_minY);
            
            //store current location in a repast structure for later use
            repast::Point<int> location(x,y);
            
            std::vector<MadAgent*> agentsInCell;
            //query four neighbouring cells, distance 0 (i.e. just the centre cell) - "true" keeps the centre cell.

            repast::VN2DGridQuery<MadAgent> VN2DQuery(space());
            VN2DQuery.query(location, 0, true, agentsInCell);

            std::vector<Human*> humans;

            for (auto a:agentsInCell){
               if (a->getId().currentRank()==rank && a->_alive){//agents must be local and living!
                //separate out humans
                if (a->getId().agentType()==_humanType) humans.push_back( (Human*) a);
               }
            }

            //humans can have one offspring per timestep - add the offspring to the model
            vector<Human*> newHumans;
            for (auto h:humans)if(h->_newH!=NULL && _reproduction){
                _context.addAgent(h->_newH);
                discreteSpace->moveTo(h->_newH->getId(), location);
                newHumans.push_back(h->_newH);
                //_totalReproductions++;
            }
            //make sure new humans can get marked for death/merged - new humans are guaranteed to be local
            for (auto n:newHumans){humans.push_back(n);}
            newHumans.clear();
            //for (auto h:humans){h->markForDeath();if (!h->_alive)_totalDied++;}


            
            //care with sync() here - need to get rid of not-alive agents:currently this is a lazy delete for new/non-local agents (they get removed one timestep late)?
            for (auto a:agentsInCell )if (!a->_alive && _death)_context.removeAgent(a->getId());//does this delete the agent storage? - yes if Boost:shared_ptr works OK

        }
    }
    if (_dispersal){
    //find out which agents need to move
    //_moved has been set to false for new agents
    //NB this has to happen after above updates to individual Humans (otherwise some cells could get mixed before other have updated, so some humans could get updated twice)
    //Note do this per cell to minimise expensive Env[x][y] lookups
    vector<Human*> movers;
    for(int x = _xlo; x < _xhi; x++){
        for(int y = _ylo; y < _yhi; y++){
            repast::Point<int> location(x,y);
            EnvironmentCell* E=_Env[x][y];
            std::vector<MadAgent*> agentsInCell;
            //query four neighbouring cells, distance 0 (i.e. just the centre cell) - "true" keeps the centre cell.
            repast::VN2DGridQuery<MadAgent> VN2DQuery(space());
            VN2DQuery.query(location, 0, true, agentsInCell);

            for (auto a:agentsInCell){
                 //dispersers must be local and alive!
                 if (a->getId().currentRank()==repast::RepastProcess::instance()->rank() && a->_alive){
                     if (a->getId().agentType()==MadModel::_humanType ) {
                         ((Human*) a)->moveIt(E,this);
                         if (a->_moved){ //humans that have changed cell
                             movers.push_back((Human *) a);a->_moved=false;
                        }
                     }
                }
            }
        }
    }
    //_totalMoved=movers.size();

    //agent data may have changed locally - ensure this is synced before anything gets moved, otherwise values do not move across threads correctly when there are buffers.
    if (buffer==1)repast::RepastProcess::instance()->synchronizeAgentStates<AgentPackage, 
            MadAgentPackageProvider, MadAgentPackageReceiver>(*provider, *receiver);

    vector<int>newPlace={0,0};
    for (auto& m:movers){
        newPlace[0]=int(m->_destination[0]);
        newPlace[1]=int(m->_destination[1]);
        //move things - these are then settled
        space()->moveTo(m,newPlace);
    }
    }
    // ***** state of the model will not be fully consistent until sync() *****
 
 sync();
    
    if (_output){

     //also get the maps
     for (auto name:outputNames)MPI_Reduce(localMaps[name].data(), outputMaps[name].data(), (_maxX-_minX+1) * (_maxY-_minY+1), MPI::DOUBLE, MPI::SUM, 0, MPI_COMM_WORLD);

     if(repast::RepastProcess::instance()->rank() == 0){netcdfOutput( CurrentTimeStep - _startingStep + 1);}
    }

    if (_restartInterval>0 && CurrentTimeStep>_restartStep && (CurrentTimeStep+1-_restartStep)%_restartInterval==0)write_restart();

}


//------------------------------------------------------------------------------------------------------------
void MadModel::sync(){
    //These lines synchronize the agents across all threads - if there is more than one...
    //Question - are threads guaranteed to be in sync?? (i.e. are we sure that all threads are on the same timestep?)
    //Possibilites for sync in terms of where to send agents are POLL, USE_CURRENT, USE_LAST_OR_CURRENT, USE_LAST_OR_POLL
    //USE_CURRENT assumes agents do not move beyond the neighbours in the cartesian grid of threads. This fails
    //for some long distance moves - POLL seems the safest (not sure if guaranteed to work though...also maybe slower)
	discreteSpace->balance();
    repast::RepastProcess::instance()->synchronizeAgentStatus<MadAgent, AgentPackage, 
             MadAgentPackageProvider, MadAgentPackageReceiver>(_context, *provider, *receiver, *receiver,RepastProcess::POLL);
    
    repast::RepastProcess::instance()->synchronizeProjectionInfo<MadAgent, AgentPackage, 
             MadAgentPackageProvider, MadAgentPackageReceiver>(_context, *provider, *receiver, *receiver,RepastProcess::POLL);

	repast::RepastProcess::instance()->synchronizeAgentStates<AgentPackage, 
             MadAgentPackageProvider, MadAgentPackageReceiver>(*provider, *receiver);
             
}
//------------------------------------------------------------------------------------------------------------
// Packages for exchanging agents across threads
//------------------------------------------------------------------------------------------------------------


MadAgentPackageProvider::MadAgentPackageProvider(repast::SharedContext<MadAgent>* agentPtr): agents(agentPtr){ }
//------------------------------------------------------------------------------------------------------------

void MadAgentPackageProvider::providePackage(MadAgent* agent, std::vector<AgentPackage>& out){
    repast::AgentId id = agent->getId();
    if (id.agentType() == MadModel::_humanType){
     AgentPackage package(id);
     ((Human*)agent)->PushThingsIntoPackage(package);
     out.push_back(package);
    }
}

//------------------------------------------------------------------------------------------------------------

void MadAgentPackageProvider::provideContent(repast::AgentRequest req, std::vector<AgentPackage>& out){
    std::vector<repast::AgentId> ids = req.requestedAgents();
    for(size_t i = 0; i < ids.size(); i++){
        providePackage(agents->getAgent(ids[i]), out);
    }
}
//------------------------------------------------------------------------------------------------------------


MadAgentPackageReceiver::MadAgentPackageReceiver(repast::SharedContext<MadAgent>* agentPtr): agents(agentPtr){}
//------------------------------------------------------------------------------------------------------------

MadAgent * MadAgentPackageReceiver::createAgent(AgentPackage package){
    repast::AgentId id=package.getId();
    if (id.agentType() == MadModel::_humanType){
        Human* c=new Human(id,package);
        return c;
    } else {
        return 0;
    }
}
//------------------------------------------------------------------------------------------------------------
//This function is needed if buffers are being used so that agents can interact across cells
void MadAgentPackageReceiver::updateAgent(AgentPackage package){
    repast::AgentId id=package.getId();
    if (id.agentType() == MadModel::_humanType){
      Human* agent = (Human*)(agents->getAgent(id));//I think this matches irrespective of the value of currentRank (AgentId== operator doesn't use it)
      agent->PullThingsOutofPackage(package);
    }
}
//------------------------------------------------------------------------------------------------------------
void MadModel::setupOutputs(){
    
	//The things added to the datasetbuilder will be accumulated over cores each timestep and output to file filename

        
    std::string filename = _filePrefix+"global.outputs"+_filePostfix+".csv";        
                
	SVDataSetBuilder svbuilder(filename, ",", repast::RepastProcess::instance()->getScheduleRunner().schedule());

    SusSum* sSum = new SusSum(this);
    svbuilder.addDataSource(repast::createSVDataSource("Total Susceptible", sSum, std::plus<int>()));
    
    InfSum* iSum = new InfSum(this);
    svbuilder.addDataSource(repast::createSVDataSource("Total Infected", iSum, std::plus<int>()));
    
    RecSum* rSum = new RecSum(this);
    svbuilder.addDataSource(repast::createSVDataSource("Total Recovered", rSum, std::plus<int>()));
    
    DeathSum* dSum = new DeathSum(this);
    svbuilder.addDataSource(repast::createSVDataSource("Total Died", dSum, std::plus<int>()));
    
    PopSum* pSum = new PopSum(this);
	svbuilder.addDataSource(repast::createSVDataSource("Total Population", pSum, std::plus<int>()));
     

	addDataSet(svbuilder.createDataSet());

}
//------------------------------------------------------------------------------------------------------------
void MadModel::setupNcOutput(){
        
        //***//
        for (auto name:outputNames)setNcGridFile(name,outputUnits[name]);

}
//------------------------------------------------------------------------------------------------------------
void MadModel::setNcGridFile(std::string GridName, std::string units){

        std::string filePath = _filePrefix+GridName+_filePostfix+".nc";
        
        netCDF::NcFile gridFile( filePath.c_str(), netCDF::NcFile::replace );             // Creates file
        auto times=TimeStep::instance()->TimeStepArray();
        netCDF::NcDim gTimeNcDim = gridFile.addDim( "time", times.size() );                    // Creates dimension
        netCDF::NcVar gTimeNcVar = gridFile.addVar( "time", netCDF::ncUint, gTimeNcDim ); // Creates variable
        gTimeNcVar.putVar( times.data() );
        gTimeNcVar.putAtt( "units", TimeStep::instance()->TimeStepUnits() );
                
        netCDF::NcDim longitudeDim =   gridFile.addDim( "Longitude", Parameters::instance()->GetLengthLongitudeArray( ) );
        netCDF::NcVar longitudeNcVar = gridFile.addVar( "Longitude", netCDF::ncFloat, longitudeDim );
        longitudeNcVar.putVar( Parameters::instance()->GetLongitudeArray( ) );
        longitudeNcVar.putAtt( "units", "degrees" );

        netCDF::NcDim latitudeDim =   gridFile.addDim( "Latitude", Parameters::instance()->GetLengthLatitudeArray( ) );
        netCDF::NcVar latitudeNcVar = gridFile.addVar( "Latitude", netCDF::ncFloat, latitudeDim );
        latitudeNcVar.putVar( Parameters::instance()->GetLatitudeArray( ) );
        latitudeNcVar.putAtt( "units", "degrees" );
                
        std::vector< netCDF::NcDim > gridDimensions={gTimeNcDim,latitudeDim,longitudeDim};

        netCDF::NcVar gridVar = gridFile.addVar(  GridName, netCDF::ncDouble, gridDimensions );
        gridVar.putAtt("units", units );
}
//------------------------------------------------------------------------------------------------------------
void MadModel::netcdfOutput( unsigned step ){

         for (auto name:outputNames)writeNcGridFile(step,outputMaps[name],name);

}
//------------------------------------------------------------------------------------------------------------

void MadModel::writeNcGridFile(unsigned step, vector<double>& GridDoubleVector,std::string GridName){
        
        std::string filePath = _filePrefix+GridName+_filePostfix+".nc";
        try {

            netCDF::NcFile gridFile( filePath.c_str(), netCDF::NcFile::write );
            netCDF::NcVar gridVar=gridFile.getVar( GridName );

            vector<size_t> pos={step,0,0};vector<size_t> num={1,Parameters::instance()->GetLengthLatitudeArray( ),Parameters::instance()->GetLengthLongitudeArray( )};
            gridVar.putVar(pos, num,GridDoubleVector.data() );
                    

        } catch( netCDF::exceptions::NcException& e ) {
                e.what( );
                std::cout << "ERROR> Write to \"" << filePath << "\" failed." << std::endl;
        }
                
}

//------------------------------------------------------------------------------------------------------------

void MadModel::dataSetClose() {
	for (size_t i = 0; i < dataSets.size(); ++i) {
		(dataSets[i])->write();
		(dataSets[i])->close();
	}
}
//---------------------------------------------------------------------------------------------------------------------------
void MadModel::addDataSet(repast::DataSet* dataSet) {
	dataSets.push_back(dataSet);
	ScheduleRunner& runner = RepastProcess::instance()->getScheduleRunner();
	runner.scheduleEvent(0.1+_startingStep-1, 1, Schedule::FunctorPtr(new MethodFunctor<repast::DataSet> (dataSet,&repast::DataSet::record)));
	Schedule::FunctorPtr dsWrite = Schedule::FunctorPtr(new MethodFunctor<repast::DataSet> (dataSet,&repast::DataSet::write));
    //output every 100 steps
	runner.scheduleEvent(100.2, 100, dsWrite);
}
//---------------------------------------------------------------------------------------------------------------------------
// Restarts
//------------------------------------------------------------------------------------------------------------
 void MadModel::write_restart(){
     //each thread writes its own restart file. Read restart is able to deal with this later...
     unsigned step=RepastProcess :: instance ()->getScheduleRunner ().currentTick ();
     std::vector<MadAgent*> agents;
     _context.selectAgents(repast::SharedContext<MadAgent>::LOCAL,agents);
     std::stringstream s;
     s<<step<<"_"<<repast::RepastProcess::instance()->rank();
     std::ofstream ofs(_filePrefix+"Restart_step_rank_"+s.str());
     //archive saves when destructor called - this block should ensure this happens
     {

      for (auto a:agents){
          if (a->_alive){

              AgentPackage package(a->getId());
              if (a->getId().agentType()==_humanType){
                  ((Human*)a)->PushThingsIntoPackage( package );
              }
              _packages.push_back(package);
          }
      }
      
      if (_archiveFormat =="binary") {boost::archive::binary_oarchive oa(ofs);oa<<_packages;}
      if (_archiveFormat =="text"  ) {boost::archive::text_oarchive   oa(ofs);oa<<_packages;}
     
      
      if (_verbose) cout<<"Wrote "<<_packages.size()<<" objects to restart: "<<"Restart_step_rank_"<<s.str()<<endl;
     }
     _packages.clear();

     
 }
//------------------------------------------------------------------------------------------------------------
void MadModel::read_restart(unsigned step){
    
    unsigned r=0,error=0,rank=repast::RepastProcess::instance()->rank();
    unsigned numProcs=repast::RepastProcess::instance()->worldSize();

    std::stringstream s;
    s<<step<<"_"<<r;
    std::string filename=_restartDirectory+"Restart_step_rank_"+s.str();
    if (rank==0){
        if (!boost::filesystem::exists(filename)){
            cout<<"No restarts found for "<<filename<<endl;
            error=1;
        }
    }
    MPI_Bcast(&error, 1, MPI_INT, 0 , MPI_COMM_WORLD);
    if (error==0){
        int nxtID[numProcs];
        for (int i=0;i<numProcs;i++)nxtID[i]=0;
        while(boost::filesystem::exists(filename)){
            
            if (rank == 0){
                
                if (error==0)cout<<"Reading restarts from "<<filename<<endl;
                std::ifstream ifs(filename);
                {// this block ensures the archive gets closed on block exit
                    try {
                        if (_archiveFormat =="binary"){boost::archive::binary_iarchive ia(ifs);ia>>_packages;}
                        if (_archiveFormat =="text"  ){boost::archive::text_iarchive   ia(ifs);ia>>_packages;}
                        
                        cout<<"Read in "<<_packages.size()<<" mad agents..."<<endl;
                        for (auto& p:_packages){
                            repast::AgentId id=p.getId();
                            int n=Human::_NextID;
                            //Human::_NextID should be incremented here so set increase flag to true
                            MadAgent* a=NULL;
                            if(id.agentType()==_humanType){
                                a=new Human( id,p,true );
                            }
                            if (a!=NULL){
                                //moveTo checks if we are local - agents may previously have been on another core, so make sure local now
                                a->getId().currentRank(0);
                                //agent unique Ids depend on id no. type and startingRank - to make sure any new agents are unique
                                //need to find the max id no. associated with a given previous rank and add 1
                                //number of threads may have decreased so check to be sure we don't overflow nxtID
                                if (a->getId().startingRank()<numProcs)
                                    nxtID[a->getId().startingRank()]=max<int>(a->getId().id()+1,nxtID[a->getId().startingRank()]);
                                _context.addAgent(a);
                                repast::Point<int> initialLocation(int(a->getLocation()[0]),int(a->getLocation()[1]));
                                space()->moveTo(id, initialLocation);
                            }else{ cout<<"Warning NULL agents on reading restart file "<<filename<<endl;}
                        }
                        _packages.clear();
                    } catch (exception e){
                        if (r==0)cout <<"************* Error attempting to open archive: wrong file format? *************";
                        error=1;
                    }
                }
            }
            
            //check for errors
            MPI_Bcast(&error, 1, MPI_INT, 0 , MPI_COMM_WORLD);
            //sync so thread 0 doesn't have to carry all the agents from a possibly multi-core previous run
            //remember every thread needs to do the sync, not just thread 0!
            sync();
            //share the nextID data - this may update progressively, as startingRanks of current live agents can be arbitrary
            //if there are more threads than previously, those with rank not previously present can have nextID=0
            //if there are fewer threads, those with startingRank>=numProcs are safe to ignore (as there will be no new agents with these startingRanks)
            MPI_Bcast(&nxtID, numProcs, MPI_INT, 0 , MPI_COMM_WORLD);
            Human::_NextID=nxtID[rank];
            
            r++;
            std::stringstream s;
            s<<step<<"_"<<r;
            filename=_restartDirectory+"Restart_step_rank_"+s.str();
        }
        
    }

    if (error!=0){
            MPI_Finalize();
            exit(error);
    }
}
//---------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------
//***------------------------------------------------TESTING Section----------------------------------------------------***//
//---------------------------------------------------------------------------------------------------------------------------
void MadModel::tests(){
    int rank = repast::RepastProcess::instance()->rank();
    randomizer* random=new RandomRepast;
    //random->SetSeed(100); seed is set from model.props file - see main.cpp
    int nranks=_dimX*_dimY;
    CohortDefinitions::Initialise(_props->getProperty("input.DataDirectory")+"/"+_props->getProperty("input.CohortDefinitionsFileName"));
    //MB main has already got the environmental data - this is stored in the background as a DataLayerSet

    //check values from Parameters have got correctly placed in class grid extents
    assert(Parameters::instance()->GetLengthLongitudeArray( )==_maxX-_minX+1);
    assert(Parameters::instance()->GetLengthLatitudeArray( )==_maxY-_minY+1);
    //now set up the environmental cells - note at present this uses the full grid, not just local to this thread
    //so that off-thread environment can be easily queried. Currently some duplication here, but it is not a huge amount of data.
    _Env=Environment(_minX,_maxX,_minY,_maxY);

    int x=_xlo,y=_ylo;
    repast::Point<int> initialLocation(x,y);
    repast::Point<int> origin(_minX,_minY);
    EnvironmentCell* E=_Env[x][y];
    
    //---------------------------------------------------
    //***-------------------TEST 1-------------------***//
    //---------------------------------------------------
    //initially all agents are just on rank 0
   int localTotals=0,globalTotals=0,n=1000;
    //create some agents on rank 0
    if (rank==0){
     cout<<"Test1: create n="<<n<<" agents on rank 0"<<endl;
     for (int i=0;i<n;i++){
      repast::AgentId id(Human::_NextID, rank, _humanType);
      id.currentRank(rank);
      Human* c = new Human(id);
      //agent id number should be increasing
      assert(c->getId().id()==i);
      c->setup(0,1, E,random);
      _context.addAgent(c);
      discreteSpace->moveTo(id, initialLocation);
      c->setLocation(x,y);
     }
    }

    std::vector<MadAgent*> agents;
    _context.selectAgents(repast::SharedContext<MadAgent>::LOCAL,agents);
    localTotals=agents.size();
    
    //check across all threads to see if agents all still exist
    
    MPI_Allreduce(&localTotals, &globalTotals, 1, MPI::INT, MPI::SUM,MPI_COMM_WORLD);
    assert(globalTotals==n);
    if (rank==0)cout<<"Test1 succeeded"<<endl;
    
    //---------------------------------------------------
    //***-------------------TEST 2-------------------***//
    //---------------------------------------------------
    //remove agents, but again only on rank 0
    agents.clear();
    _context.selectAgents(repast::SharedContext<MadAgent>::LOCAL,agents);
    if (rank==0){
     cout<<"Test2: delete all "<<n<<" agents on rank 0"<<endl;
     for (auto a:agents)_context.removeAgent(a->getId());
    }
    agents.clear();
    _context.selectAgents(repast::SharedContext<MadAgent>::LOCAL,agents);
    localTotals=agents.size();
    MPI_Allreduce(&localTotals, &globalTotals, 1, MPI::INT, MPI::SUM,MPI_COMM_WORLD);
    assert(globalTotals==0);
    if (rank==0)cout<<"Test2 succeeded "<<endl;

    //---------------------------------------------------
    //***-------------------TEST 3-------------------***//
    //---------------------------------------------------
    //add agents back to rank 0
    if (rank==0){
     cout<<"Test3: re-create "<<n<<" agents on rank 0"<<endl;
     for (int i=0;i<n;i++){
      repast::AgentId id(Human::_NextID, rank, _humanType);
      id.currentRank(rank);
      Human* c = new Human(id);
      c->setup(0,1, E,random);
      _context.addAgent(c);
      discreteSpace->moveTo(id, initialLocation);
      c->setLocation(x,y);
     }
    }
    agents.clear();
    _context.selectAgents(repast::SharedContext<MadAgent>::LOCAL,agents);
    localTotals=agents.size();
    MPI_Allreduce(&localTotals, &globalTotals, 1, MPI::INT, MPI::SUM,MPI_COMM_WORLD);
    assert(globalTotals==n);
    if (rank==0)cout<<"Test3 succeeded "<<endl;
    
    //---------------------------------------------------
    //***-------------------TEST 4-------------------***//
    //---------------------------------------------------
    //all agents still exist: allow any thread to move them
    agents.clear();
    _context.selectAgents(repast::SharedContext<MadAgent>::LOCAL,agents);
    if (rank==0)cout<<"Test4: move agents by 1 unit in x "<<endl;
    std::vector<int> displ{1,0};
    for (auto a:agents)space()->moveByDisplacement(a,displ);
    sync();
    agents.clear();
    _context.selectAgents(repast::SharedContext<MadAgent>::LOCAL,agents);
    localTotals=agents.size();
    MPI_Allreduce(&localTotals, &globalTotals, 1, MPI::INT, MPI::SUM,MPI_COMM_WORLD);
    assert(globalTotals==n);
    vector<int> location;
    for (auto a:agents){
      space()->getLocation(a->getId(), location);
      assert(location[0]==_xlo+1 && location[1]==_ylo+0);
    }
    if (rank==0)cout<<"Test4 succeeded "<<endl;
    //---------------------------------------------------
    //***-------------------TEST 5-------------------***//
    //---------------------------------------------------
    //all agents still exist: allow any thread to move them
    agents.clear();
    _context.selectAgents(repast::SharedContext<MadAgent>::LOCAL,agents);
    if (rank==0)cout<<"Test5: move agents by 1 unit in y "<<endl;
    std::vector<int> displ2{0,1};
    for (auto a:agents)space()->moveByDisplacement(a,displ2);
    sync();
    agents.clear();
    _context.selectAgents(repast::SharedContext<MadAgent>::LOCAL,agents);
    localTotals=agents.size();
    MPI_Allreduce(&localTotals, &globalTotals, 1, MPI::INT, MPI::SUM,MPI_COMM_WORLD);
    assert(globalTotals==n);
    
    for (auto a:agents){
      space()->getLocation(a->getId(), location);
      assert(location[0]==_xlo+1 && location[1]==_ylo+1);
    }
    if (rank==0)cout<<"Test5 succeeded "<<endl;
    //---------------------------------------------------
    //***-------------------TEST 6-------------------***//
    //---------------------------------------------------
    //all agents still exist: allow any thread to move them
    //so moves of any distance will work only if grid is wrapped.
    //all agents are still co-located, so should end up on the same thread
    agents.clear();
    _context.selectAgents(repast::SharedContext<MadAgent>::LOCAL,agents);
    int xm1=2,ym1=-5;
    if (rank==0)cout<<"Test6: move agents by "<<xm1<<" units in x and "<<ym1<<" units in y "<<endl;
    std::vector<int> displ3{xm1,ym1};
    for (auto a:agents)space()->moveByDisplacement(a,displ3);
    sync();
    agents.clear();
    _context.selectAgents(repast::SharedContext<MadAgent>::LOCAL,agents);
    localTotals=agents.size();
    MPI_Allreduce(&localTotals, &globalTotals, 1, MPI::INT, MPI::SUM,MPI_COMM_WORLD);
    assert(globalTotals==n);
    for (auto a:agents){
      space()->getLocation(a->getId(), location);
      //calculate the wrapped position including previous move by +1,+1 
      int xpos=xm1+1+_minX,ypos=ym1+1+_minY;
      while (xpos<_minX)xpos+=(_maxX - _minX + 1);
      while (ypos<_minY)ypos+=(_maxY - _minY + 1);
      xpos=xpos % (_maxX - _minX + 1);
      ypos=ypos % (_maxY - _minY + 1);
      assert(location[0]==xpos && location[1]==ypos);
    }
    if (rank==0)cout<<"Test6 succeeded "<<endl;
    //---------------------------------------------------
    //***-------------------TEST 7-------------------***//
    //---------------------------------------------------
    //all agents still exist: allow any thread to move them
    //note grid is wrapped, so moves of any distance should work.
    //all agents are still co-located, so should end up on the same thread
    agents.clear();
    _context.selectAgents(repast::SharedContext<MadAgent>::LOCAL,agents);
    int xm2=-79,ym2=3355;
    if (rank==0)cout<<"Test7: move agents by "<<xm2<<" units in x and "<<ym2<<" units in y "<<endl;
    std::vector<int> displ4{xm2,ym2};
    for (auto a:agents)space()->moveByDisplacement(a,displ4);
    sync();
    agents.clear();
    _context.selectAgents(repast::SharedContext<MadAgent>::LOCAL,agents);
    localTotals=agents.size();
    MPI_Allreduce(&localTotals, &globalTotals, 1, MPI::INT, MPI::SUM,MPI_COMM_WORLD);
    assert(globalTotals==n);
    for (auto a:agents){
      space()->getLocation(a->getId(), location);
      //calculate the wrapped position including previous move by xm1+1,ym1+1 
      int xpos=xm1+xm2+1+_minX,ypos=ym1+ym2+1+_minY;
      while (xpos<_minX)xpos+=(_maxX - _minX + 1);
      while (ypos<_minY)ypos+=(_maxY - _minY + 1);
      xpos=xpos % (_maxX - _minX + 1);
      ypos=ypos % (_maxY - _minY + 1);
      assert(location[0]==xpos && location[1]==ypos);
    }
    //---------------------------------------------------
    //***-------------------TEST 8-------------------***//
    //---------------------------------------------------
    //move all agents back to origin (thread 0)
    agents.clear();
    _context.selectAgents(repast::SharedContext<MadAgent>::LOCAL,agents);

    if (rank==0)cout<<"Test8: move agents back to origin: "<<_minX<<" "<<_minY<<endl;
    for (auto a:agents)space()->moveTo(a,origin);
    sync();
    agents.clear();
    _context.selectAgents(repast::SharedContext<MadAgent>::LOCAL,agents);
    localTotals=agents.size();
    MPI_Allreduce(&localTotals, &globalTotals, 1, MPI::INT, MPI::SUM,MPI_COMM_WORLD);
    assert(globalTotals==n);
    for (auto a:agents){
      space()->getLocation(a->getId(), location);
      assert(location[0]==_minX && location[1]==_minY);
    }
    if (rank==0)cout<<"Test8 succeeded "<<endl;
    //---------------------------------------------------
    //***-------------------TEST 9-------------------***//
    //---------------------------------------------------
    //move all agents to random locations
    agents.clear();
    _context.selectAgents(repast::SharedContext<MadAgent>::LOCAL,agents);

    if (rank==0)cout<<"Test9: move agents to random location: "<<endl;
    
    for (auto a:agents){
        vector<int> loc{int(random->GetUniform()*1000.),int(random->GetUniform()*1000.)};
        space()->moveTo(a,loc);
    }
    sync();
    agents.clear();
    _context.selectAgents(repast::SharedContext<MadAgent>::LOCAL,agents);
    localTotals=agents.size();
    cout<<"Total on rank:"<<rank<<" "<<localTotals<<endl;
    MPI_Allreduce(&localTotals, &globalTotals, 1, MPI::INT, MPI::SUM,MPI_COMM_WORLD);
    assert(globalTotals==n);
    for (auto a:agents){
      space()->getLocation(a->getId(), location);
      assert(location[0]>=_xlo && location[1]>=_ylo && location[0]<_xhi && location[1]<_yhi);
    }
    if (rank==0)cout<<"Test9 succeeded "<<endl;
    //---------------------------------------------------
    //***-------------------TEST 10-------------------***//
    //---------------------------------------------------
    //test add and delete on all threads, not just rank 0
    agents.clear();
    _context.selectAgents(repast::SharedContext<MadAgent>::LOCAL,agents);

    if (rank==0)cout<<"Test10: delete up to 10 agents per thread, add up to 100 new, move all randomly: "<<endl;
    int nr=0,globalAdded,globalRemoved;
    //remove up to 10 agents locally: move other agents around randomly
    for (auto a:agents){
         if (nr<int(random->GetUniform()*10.)){a->_alive=false;_context.removeAgent(a->getId());nr++;}//care here - after the removal, the pointer to a is no longer valid
         else{
             vector<int> loc{int(random->GetUniform()*781.),int(random->GetUniform()*500.-250)};
             space()->moveTo(a,loc);
         }
    }
    //addup to 100 new agents on this thread and then move these at random 
    int nnew=int(random->GetUniform()*100.);
    for (int i=0;i<nnew;i++){
      repast::AgentId id(Human::_NextID, rank, _humanType);
      id.currentRank(rank);
      Human* c = new Human(id);
      c->setup(0,1, E,random);
      _context.addAgent(c);
      vector<int> loc{int(random->GetUniform()*37.-99),int(random->GetUniform()*188.)};
      space()->moveTo(c,loc);
    }

    sync();
    MPI_Allreduce(&nnew, &globalAdded, 1, MPI::INT, MPI::SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&nr, &globalRemoved, 1, MPI::INT, MPI::SUM,MPI_COMM_WORLD);

    agents.clear();
    _context.selectAgents(repast::SharedContext<MadAgent>::LOCAL,agents);
    localTotals=agents.size();
    MPI_Allreduce(&localTotals, &globalTotals, 1, MPI::INT, MPI::SUM,MPI_COMM_WORLD);

    assert(globalTotals==n+globalAdded-globalRemoved);
    for (auto a:agents){
      space()->getLocation(a->getId(), location);
      assert(location[0]>=_xlo && location[1]>=_ylo && location[0]<_xhi && location[1]<_yhi);
      assert(a->getId().currentRank()==rank);
    }
    if (rank==0)cout<<"Test10 succeeded: "<<"removed:"<<globalRemoved<<" added:"<<globalAdded<<endl;
    //---------------------------------------------------
    //***-------------------TEST 11-------------------***//
    //---------------------------------------------------
    //test whether data in cohorts move correctly across threads
    //first remove all existing agents
    agents.clear();
    _context.selectAgents(repast::SharedContext<MadAgent>::LOCAL,agents);
    if (rank==0)cout<<"Test11: check agent data moves correctly across threads"<<endl;

    for (auto a:agents){
        _context.removeAgent(a->getId());
    }
    sync();
    //now add a new agent on each thread and set its properties to known values
    for (int i=0;i<n;i++){
      repast::AgentId id(Human::_NextID, rank, _humanType);
      id.currentRank(rank);
      Human* c = new Human(id);
      c->setup(0,1, E,random);
      _context.addAgent(c);
      discreteSpace->moveTo(id, initialLocation);
      c->setLocation(x,y);
      setupHumanTestValues(c);
      checkHumanTestValues(c);
     }
    //move agents randomly 100 times
    for (int i=0;i<100;i++){
    agents.clear();
    _context.selectAgents(repast::SharedContext<MadAgent>::LOCAL,agents);
    for (auto a:agents){
            vector<int> loc{int(random->GetUniform()*237.-125),int(random->GetUniform()*981.-250)};
            space()->moveTo(a,loc);
    }
    sync();}
    agents.clear();
    _context.selectAgents(repast::SharedContext<MadAgent>::LOCAL,agents);
    for (auto a:agents){
      checkHumanTestValues((Human*) a);
    }
    //change some values and check they are still OK
    for (auto a:agents){
      ((Human*)a)->setLocation(112,-98);
    }

   //if there's a buffer, and you change local agent data, seems you need a sync of states *before* you do anything like a move...
   //sync after a move leads to errors
   repast::RepastProcess::instance()->synchronizeAgentStates<AgentPackage,MadAgentPackageProvider, MadAgentPackageReceiver>(*provider, *receiver);
   for (auto a:agents){
           vector<int> loc{int(random->GetUniform()*237.-125),int(random->GetUniform()*981.-250)};
           space()->moveTo(a,loc);
   }
   sync();
   agents.clear();
   _context.selectAgents(repast::SharedContext<MadAgent>::LOCAL,agents);
    for (auto a:agents){
      assert(((Human*)a)->_location[0]==112);
      assert(((Human*)a)->_location[1]==-98);
    }
    if (rank==0)cout<<"Test11 succeeded: values copied across threads unchanged"<<endl;

    cout.flush();
    sync();
    //---------------------------------------------------
    //***-------------------TEST 12-------------------***//
    //---------------------------------------------------

    //test whether data in cohorts can be saved out to a boost archive file (intended for model restarts from a saved state -although randoms will make this non-reproducible)
    agents.clear();
    _context.selectAgents(repast::SharedContext<MadAgent>::LOCAL,agents);
    if (rank==0)cout<<"Test12: check agent data saves correctly to a file and can be restored"<<endl;

    for (auto a:agents){
        _context.removeAgent(a->getId());
    }
    sync();
    //now add a new agent on each thread and set its properties to known values
    for (int i=0;i<n;i++){
      repast::AgentId id(Human::_NextID, rank, _humanType);
      id.currentRank(rank);
      Human* c = new Human(id);
      c->setup(0,1, E,random);
      _context.addAgent(c);
      discreteSpace->moveTo(id, initialLocation);
      c->setLocation(x,y);
      setupHumanTestValues(c);
      checkHumanTestValues(c);
     }
     sync();
     agents.clear();
     _context.selectAgents(repast::SharedContext<MadAgent>::LOCAL,agents);
     std::stringstream s;
     s<<rank;
     std::ofstream ofs("TestAgentSerialization_rank_"+s.str());
     {
      boost::archive::text_oarchive oa(ofs);
      for (auto a:agents){
         AgentPackage package(a->getId());
         ((Human*)a)->PushThingsIntoPackage( package );
         _packages.push_back(package);
      }
      oa<<_packages;
     }
     cout<<"Wrote "<<_packages.size()<<" cohorts"<<endl;
     _packages.clear();
     std::ifstream ifs("TestAgentSerialization_rank_"+s.str());

     {
      boost::archive::text_iarchive ia(ifs);
      ia>>_packages;
      cout<<"Read "<<_packages.size()<<" cohorts"<<endl;
      for (auto& p:_packages){
         repast::AgentId id=p.getId();
         int n=Human::_NextID;
         //Human::_NextID should be incremented here so set increase flag to true
         MadAgent* a=new Human( id,p,true );
         assert(Human::_NextID==n+1);
         assert(a->getId()==p._id);
         assert(a->getId().currentRank()==rank);
         assert(a->getId().agentType()==_humanType);
         checkHumanTestValues((Human*)a);
      }
     }
     if (rank==0)cout<<"Test 12 : succeeded"<<endl;

    //---------------------------------------------------
    //***-------------------Finished TESTS-------------------***//
    //---------------------------------------------------
 
}
//---------------------------------------------------------------------------------------------------------------------------
//define some data values for the Cohort to check whether they are preserved on moving across threads
//see Cohort::setup
void MadModel::setupHumanTestValues(Human* c){

    c->_diseases["covid"].infect();
    c->_diseases["flu"];
/*
    c->_Merged                      = false;
    c->_alive                       = true;
    

    c->_moved=false;
    c->_location={-12,75};
    */
}
//---------------------------------------------------------------------------------------------------------------------------
//check that the values set in the above function are still maintained
void MadModel::checkHumanTestValues(Human* c){

    assert(c->_diseases["covid"].infected());
    assert(!c->_diseases["flu"].infected());
    assert(c->_Realm      =="terrestrial");
/*

    assert(c->_alive                       == true);


    assert(c->_moved==false);
    assert(c->_location[0]==-12);
    assert(c->_location[1]==75);
    */
}
//---------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------
//***------------------------------------------------END TESTING Section----------------------------------------------------***//
//---------------------------------------------------------------------------------------------------------------------------
