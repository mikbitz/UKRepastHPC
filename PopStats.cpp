/*
 * PopStats.cpp
 * Created on: April 5 2020
 * Author: Mike Bithell
 *
 *   Based on Infectionsum.cpp
 * Repast for High Performance Computing (Repast HPC)
 *
 *   Copyright (c) 2010 Argonne National Laboratory
 *   All rights reserved.
 *  
 *   Redistribution and use in source and binary forms, with 
 *   or without modification, are permitted provided that the following 
 *   conditions are met:
 *  
 *  	 Redistributions of source code must retain the above copyright notice,
 *  	 this list of conditions and the following disclaimer.
 *  
 *  	 Redistributions in binary form must reproduce the above copyright notice,
 *  	 this list of conditions and the following disclaimer in the documentation
 *  	 and/or other materials provided with the distribution.
 *  
 *  	 Neither the name of the Argonne National Laboratory nor the names of its
 *     contributors may be used to endorse or promote products derived from
 *     this software without specific prior written permission.
 *  
 *   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *   ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 *   PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE TRUSTEES OR
 *   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 *   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 *   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 *   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
 *   EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *  Created on: Sep 2, 2010
 *      Author: nick
 */

#include "PopStats.h"
#include "model.h"

PopSum::PopSum(MadModel* Mobs) : obs(Mobs){}

PopSum::~PopSum() {}

int PopSum::getData() {
	return obs->PopCount();
}
//----------------------------------------------------------
SusSum::SusSum(MadModel* Mobs) : obs(Mobs){}

SusSum::~SusSum() {}

int SusSum::getData() {
	return obs->SusCount();
}
//----------------------------------------------------------
InfSum::InfSum(MadModel* Mobs) : obs(Mobs){}

InfSum::~InfSum() {}

int InfSum::getData() {
	return obs->InfCount();
}
//----------------------------------------------------------
RecSum::RecSum(MadModel* Mobs) : obs(Mobs){}

RecSum::~RecSum() {}

int RecSum::getData() {
	return obs->RecCount();
}
//----------------------------------------------------------
DeathSum::DeathSum(MadModel* Mobs) : obs(Mobs){}

DeathSum::~DeathSum() {}

int DeathSum::getData() {
	return obs->DeathCount();
}
