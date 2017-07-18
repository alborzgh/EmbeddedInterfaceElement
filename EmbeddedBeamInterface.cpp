/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */
                                                                              
// Written: Alborz Ghofrani, Pedro Arduino, U.Washington
// Created: Oct 2014
// Description: This file contains the class definition for EmbeddedBeamInterface.

#include <EmbeddedBeamInterface.h>
#include <Node.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <Renderer.h>
#include <Domain.h>
#include <string.h>
#include <Information.h>
#include <Parameter.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <ElementResponse.h>
#include <elementAPI.h>
#include <cmath>

static int num_EmbeddedBeamInterface = 0;
Matrix       EmbeddedBeamInterface::m_ContactStiffness(QBEC_NUM_DOF,QBEC_NUM_DOF);
Vector       EmbeddedBeamInterface::m_ContactForces(QBEC_NUM_DOF);
const double EmbeddedBeamInterface::m_Pi = 3.14159265359;

void *
OPS_EmbeddedBeamInterface(void)  
{
	if (num_EmbeddedBeamInterface == 0) {
        num_EmbeddedBeamInterface++;
	opserr<<"EmbeddedBeamInterface element - Written: A.Ghofrani, D.Turello, P.Arduino, U.Washington\n";
	}
	
	Element *theElement = 0;

	int numArgs = OPS_GetNumRemainingInputArgs();
	if (numArgs < 1) {
		opserr << "Want: EmbeddedBeamInterface tag? \n";
		return 0;
	}
    
	int iData[1];
	int eleTag = 0;
	int numData = 1;
	if (OPS_GetIntInput(&numData, iData) != 0) {
		opserr << "WARNING invalid integer data: element EmbeddedBeamInterface" << endln;
		return 0;
	}
		
	eleTag = iData[0];
	
	theElement = new EmbeddedBeamInterface(iData[0]);

	if (theElement == 0) {
		opserr << "WARNING could not create element of type EmbeddedBeamInterface\n";
		return 0;
	}

	return theElement;
}


EmbeddedBeamInterface::EmbeddedBeamInterface(int tag): Element(tag, ELE_TAG_EmbeddedBeamInterface),
		 externalNodes(6)
{
	for(int i = 0; i < 6; i++) 
		theNodes[i] = 0;
}

EmbeddedBeamInterface::EmbeddedBeamInterface()
		: Element(0, ELE_TAG_EmbeddedBeamInterface),
		 externalNodes(6)
{
	for(int i = 0; i < 6; i++) 
		theNodes[i] = 0;
}

EmbeddedBeamInterface::~EmbeddedBeamInterface()
{

}

int 
EmbeddedBeamInterface::getNumExternalNodes(void) const
{
	return QBEC_NUM_NODE;
}

const ID&
EmbeddedBeamInterface::getExternalNodes(void)
{
	return externalNodes;
}
    
Node **
EmbeddedBeamInterface::getNodePtrs(void)
{
	return theNodes;
}

int 
EmbeddedBeamInterface::getNumDOF(void)
{
	return QBEC_NUM_DOF;
}
    
void 
EmbeddedBeamInterface::setDomain(Domain *theDomain)
{
	this->DomainComponent::setDomain(theDomain);
	return;
}

int 
EmbeddedBeamInterface::commitState(void)
{
	int retVal = 0;

	// call element commitState to do any base class stuff
	if ((retVal = this->Element::commitState()) != 0) {
	    opserr << "EmbeddedBeamInterface::commitState() - failed in base class";
	}
	
	return retVal;
}

int 
EmbeddedBeamInterface::revertToLastCommit(void)
{
	return 0;
}

int 
EmbeddedBeamInterface::revertToStart(void)
{
	return 0;
}


const Matrix&
EmbeddedBeamInterface::getTangentStiff(void)
{
	m_ContactStiffness.Zero();
	return m_ContactStiffness;
}

const Matrix&
EmbeddedBeamInterface::getInitialStiff(void)
{
	return this->getTangentStiff();
}

const Vector&
EmbeddedBeamInterface::getResistingForce(void)
{
	m_ContactForces.Zero();
	return m_ContactForces;
}

int 
EmbeddedBeamInterface::sendSelf(int commitTag, Channel &theChannel)
{
	return 0;
}

int 
EmbeddedBeamInterface::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker 
		  &theBroker)
{
	return 0;
}

int 
EmbeddedBeamInterface::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
{
	return 0;
}
 
void 
EmbeddedBeamInterface::Print(OPS_Stream &s, int flag)
{
	return;
}

Response*
EmbeddedBeamInterface::setResponse(const char **argv, int argc, 
			  OPS_Stream &s)
{
	return 0;
}

int 
EmbeddedBeamInterface::getResponse(int responseID, Information &eleInformation)
{	
	return 0;
}

int 
EmbeddedBeamInterface::setParameter(const char **argv, int argc, Parameter &param)
{
	return 0;
}

int
EmbeddedBeamInterface::updateParameter(int parameterID, Information &info)
{
	return 0;
}

int
EmbeddedBeamInterface::update(void)
{
	return 0;
}