/*
 *  disease.h
 *  Created on: April 5, 2020
 *      Author: Mike Bithell
 * 
 */

#ifndef DISEASE_H
#define DISEASE_H
#include <boost/serialization/access.hpp>
//serialize method is included so that agents can serialize diseases in AgentPackage.h
//NB to get boost serialize to a file to work data here has to be initialized - otherwise it crashes with an error on archive input.
class disease {
public:
    disease();
    void infect();
    bool infected();
    bool infectious();
    void becomeInfectious();
    double infectionProb();
    void recover();
    bool recovered();
    void update();
private:
    bool _infected=false;
    bool _recovered=false;
    bool _infectious=false;
    double _infectionProb=0.5;
    unsigned _timer=0;
    friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive& ar, const unsigned int version) {
        ar & _infected;
        ar & _recovered;
        ar & _infectious;
        ar & _timer;
        ar & _infectionProb;

    }
};

#endif /* DISEASE_H */
