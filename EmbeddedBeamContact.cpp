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

// Written: Alborz Ghofrani, Diego Turello, Pedro Arduino, U.Washington 
// Created: May 2017
// Description: This file contains the class definition for EmbeddedBeamContact.

#include <EmbeddedBeamContact.h>
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
#include <CrdTransf.h>
#include <elementAPI.h>
#include <cmath>

static int num_EmbeddedBeamContact = 0;
Matrix       EmbeddedBeamContact::m_InterfaceStiffness(EBC_NUM_DOF, EBC_NUM_DOF);
Vector       EmbeddedBeamContact::m_InterfaceForces(EBC_NUM_DOF);
const double EmbeddedBeamContact::m_Pi = 3.14159265359;

void *
OPS_EmbeddedBeamContact(void)
{
    if (num_EmbeddedBeamContact == 0) {
        num_EmbeddedBeamContact++;
        opserr << "EmbeddedBeamContact element - Written: A.Ghofrani, D.Turello, P.Arduino, U.Washington\n";
    }

    Element *theElement = 0;

    int numArgs = OPS_GetNumRemainingInputArgs();
    if (numArgs < 1) {
        opserr << "Want: EmbeddedBeamContact tag? \n";
        return 0;
    }

    int iData[1];
    int eleTag = 0;
    int numData = 1;
    if (OPS_GetIntInput(&numData, iData) != 0) {
        opserr << "WARNING invalid integer data: element EmbeddedBeamContact" << endln;
        return 0;
    }

    eleTag = iData[0];

    theElement = new EmbeddedBeamContact(iData[0]);

    if (theElement == 0) {
        opserr << "WARNING could not create element of type EmbeddedBeamContact\n";
        return 0;
    }

    return theElement;
}


EmbeddedBeamContact::EmbeddedBeamContact(int tag) : Element(tag, ELE_TAG_EmbeddedBeamContact),
externalNodes(10)
{

}

EmbeddedBeamContact::EmbeddedBeamContact(int tag, int beamTag, int solidTag, int crdTransfTag,
    double beamRho, double beamTheta, double solidXi, double solidEta, double solidZeta,
    double radius, double area) : Element(tag, ELE_TAG_EmbeddedBeamContact),
    externalNodes(10),
    theSolidTag(solidTag), theBeamTag(beamTag),
    m_solid_xi(solidXi), m_solid_eta(solidEta), m_solid_zeta(solidZeta),
    m_beam_rho(beamRho), m_beam_theta(beamTheta), m_beam_radius(radius),
    m_area(area), m_ep(1.0e15),
    m_Ba_rot_n(3), m_Bb_rot_n(3),
    m_Ba_disp_n(3), m_Bb_disp_n(3),
    m_Ba1(3), m_Bb1(3),
    m_Bcl_pos(3), m_Bcl_pos_n(3), m_pos(3),
    m_B_loc(3), m_S_disp(3),
    mQa(3, 3), mQb(3, 3), mQc(3, 3), mc1(3),
    mBphi(3, 12), mBu(3, 12), m_Ns(8), m_n(3)
{
    // get domain to access element tags and their nodes
#ifdef _PARALLEL_PROCESSING
#include <PartitionedDomain.h>
    extern PartitionedDomain theDomain;
#else
    extern Domain theDomain;
#endif

    // set element nodes using the solid and beam tags
    ID tempNodes = theDomain.getElement(solidTag)->getExternalNodes();
    for (int ii = 0; ii < 8; ii++)
        externalNodes(ii) = tempNodes(ii);
    tempNodes = theDomain.getElement(beamTag)->getExternalNodes();
    externalNodes(8) = tempNodes(0);
    externalNodes(9) = tempNodes(1);

    for (int i = 0; i < 10; i++)
        theNodes[i] = 0;

    // get the coordinate transformation object
    crdTransf = OPS_GetCrdTransf(crdTransfTag)->getCopy3d();

    //opserr << "Beam tag : " << theBeamTag << ", Solid tag : " << theSolidTag << ", Beam Iso Coordinates = (" << m_beam_rho << "," << m_beam_theta << 
    //  "), Solid Iso Coordinates = (" << m_solid_xi << "," << m_solid_eta << "," << m_solid_zeta << ")" << endln;
}

EmbeddedBeamContact::EmbeddedBeamContact()
    : Element(0, ELE_TAG_EmbeddedBeamContact),
    externalNodes(10)
{
    for (int i = 0; i < 10; i++)
        theNodes[i] = 0;
}

EmbeddedBeamContact::~EmbeddedBeamContact()
{

}

int
EmbeddedBeamContact::getNumExternalNodes(void) const
{
    return EBC_NUM_NODE;
}

const ID&
EmbeddedBeamContact::getExternalNodes(void)
{
    return externalNodes;
}

Node **
EmbeddedBeamContact::getNodePtrs(void)
{
    return theNodes;
}

int
EmbeddedBeamContact::getNumDOF(void)
{
    return EBC_NUM_DOF;
}

int
EmbeddedBeamContact::revertToLastCommit(void)
{
    return 0;
}

int
EmbeddedBeamContact::revertToStart(void)
{
    return 0;
}


const Matrix&
EmbeddedBeamContact::getTangentStiff(void)
{
    m_InterfaceStiffness.Zero();
    
    if (inContact)
    {
        Vector c2(3), c3(3);
        for (int ii = 0; ii < 3; ii++)
        {
            c2(ii) = mQc(ii, 1);
            c3(ii) = mQc(ii, 2);
        }

        Vector temp(24), temp1(12);

        for (int ii = 0; ii < 8; ii++)
            for (int jj = 0; jj < 3; jj++)
                temp(3 * ii + jj) = m_Ns(ii) * m_n(jj);

        Matrix temp2(3, 12), temp2T(12,3);
        temp2 = mBu - (m_beam_radius*(cos(m_beam_theta)*ComputeSkew(c2) + sin(m_beam_theta)*ComputeSkew(c3))) * mBphi;
        temp2T = Transpose(3, 12, temp2);
        temp1 = temp2T * m_n;

        for (int ii = 0; ii < 24; ii++)
            for (int jj = 0; jj < 24; jj++)
                m_InterfaceStiffness(ii, jj) = temp(ii) * temp(jj);

        for (int ii = 0; ii < 24; ii++)
            for (int jj = 0; jj < 12; jj++)
                m_InterfaceStiffness(ii, 24 + jj) = m_InterfaceStiffness(24 + jj, ii) = -1.0  * temp(ii) * temp1(jj);
        
        for (int ii = 0; ii < 12; ii++)
            for (int jj = 0; jj < 12; jj++)
                m_InterfaceStiffness(24 + ii, 24 + jj) = temp1(ii) * temp1(jj);


        m_InterfaceStiffness *= m_ep * m_area;
    }

    return m_InterfaceStiffness;
}

const Matrix&
EmbeddedBeamContact::getInitialStiff(void)
{
    return this->getTangentStiff();
}

const Vector&
EmbeddedBeamContact::getResistingForce(void)
{
    m_InterfaceForces.Zero();

    if (inContact)
    {
        double temp;
        temp = (m_pos + m_S_disp - m_B_loc) ^ m_n;
        temp *= m_ep * m_area;

        Vector c2(3), c3(3);
        for (int ii = 0; ii < 3; ii++)
        {
            c2(ii) = mQc(ii, 1);
            c3(ii) = mQc(ii, 2);
        }

        Matrix temp2(3, 12);
        temp2 = mBu - (m_beam_radius*(cos(m_beam_theta)*ComputeSkew(c2) + sin(m_beam_theta)*ComputeSkew(c3))) * mBphi;


        for (int ii = 0; ii < 8; ii++)
            for (int jj = 0; jj < 3; jj++)
                m_InterfaceForces(3 * ii + jj) = temp * m_Ns(ii) * m_n(jj);

        for (int ii = 0; ii < 12; ii++)
            for (int jj = 0; jj < 3; jj++)
                m_InterfaceForces(24 + ii) += -1.0 * temp * temp2(jj, ii) * m_n(jj);
    }

    return m_InterfaceForces;
}

int
EmbeddedBeamContact::sendSelf(int commitTag, Channel &theChannel)
{
    return 0;
}

int
EmbeddedBeamContact::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker
    &theBroker)
{
    return 0;
}

int
EmbeddedBeamContact::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
{
    return 0;
}

void
EmbeddedBeamContact::Print(OPS_Stream &s, int flag)
{
    return;
}

Response*
EmbeddedBeamContact::setResponse(const char **argv, int argc,
    OPS_Stream &s)
{
    if (strcmp(argv[0], "locationBeam") == 0 || strcmp(argv[0], "locBeam") == 0) {
        return new ElementResponse(this, 1, Vector(3));

    }
    else if (strcmp(argv[0], "locationSolid") == 0 || strcmp(argv[0], "locSolid") == 0) {
        return new ElementResponse(this, 2, Vector(3));

    }
    else if (strcmp(argv[0], "beamCL") == 0 || strcmp(argv[0], "beamCenterLine") == 0) {
        return new ElementResponse(this, 3, Vector(3));

    }
    else if (strcmp(argv[0], "c1") == 0 || strcmp(argv[0], "tangent") == 0) {
        return new ElementResponse(this, 4, Vector(3));

    }
    else if (strcmp(argv[0], "c2") == 0 || strcmp(argv[0], "perp2") == 0) {
        return new ElementResponse(this, 5, Vector(3));

    }
    else if (strcmp(argv[0], "c3") == 0 || strcmp(argv[0], "perp3") == 0) {
        return new ElementResponse(this, 6, Vector(3));

    }
    else if (strcmp(argv[0], "force") == 0 || strcmp(argv[0], "globalForce") == 0) {
        return new ElementResponse(this, 7, Vector(3));

    }
    else {
        opserr << "EmbeddedBeamContact Recorder, " << argv[0] << "is an unknown recorder request"
            << "  Element tag : " << this->getTag() << endln;
        return 0;
    }
}

int
EmbeddedBeamContact::getResponse(int responseID, Information &eleInformation)
{
    if (responseID == 1) { // location

        return eleInformation.setVector(m_B_loc);

    }
    else if (responseID == 2) { // location

        return eleInformation.setVector(m_pos + m_S_disp);
    }
    else if (responseID == 3) { // centerline

        return eleInformation.setVector(m_Bcl_pos);
    }
    else if (responseID == 4) { // c1
        Vector temp(3);
        for (int ii = 0; ii < 3; ii++)
            temp(ii) = mQc(ii, 0);
        return eleInformation.setVector(temp);
    }
    else if (responseID == 5) { // c2
        Vector temp(3);
        for (int ii = 0; ii < 3; ii++)
            temp(ii) = mQc(ii, 1);
        return eleInformation.setVector(temp);
    }
    else if (responseID == 6) { // c3
        Vector temp(3);
        for (int ii = 0; ii < 3; ii++)
            temp(ii) = mQc(ii, 2);
        return eleInformation.setVector(temp);
    }
    else if (responseID == 7) { // contact force
        Vector temp(3);
        temp = m_ep * m_area * (m_pos + m_S_disp - m_B_loc);
        return eleInformation.setVector(temp);
    }
    else {
        opserr << "EmbeddedBeamContact, tag = " << this->getTag()
            << " -- unknown request" << endln;
        return -1;
    }
}

int
EmbeddedBeamContact::setParameter(const char **argv, int argc, Parameter &param)
{
    return 0;
}

int
EmbeddedBeamContact::updateParameter(int parameterID, Information &info)
{
    return 0;
}

void
EmbeddedBeamContact::setDomain(Domain *theDomain)
{
    for (int ii = 0; ii < 10; ii++)
    {
        theNodes[ii] = theDomain->getNode(externalNodes(ii));
        if (theNodes[ii] == 0)
        {
            opserr << "Could not find node " << externalNodes(ii) << "." << endln;
            return;
        }
        if ((theNodes[ii]->getNumberDOF() != 3) && (ii < 8))
        {
            opserr << "Solid node " << externalNodes(ii) << " has to have 3 degrees of freedom." << endln;
            return;
        }
        if ((theNodes[ii]->getNumberDOF() != 6) && (ii > 7))
        {
            opserr << "Beam node " << externalNodes(ii) << " has to have 6 degrees of freedom." << endln;
            return;
        }
    }

    // initialize the transformation
    if (crdTransf->initialize(theNodes[8], theNodes[9]))
    {
        opserr << "EmbeddedBeamContact::setDomain(): Error initializing coordinate transformation";
        return;
    }

    m_beam_length = crdTransf->getInitialLength();
    if (m_beam_length < 1.0e-12) {
        opserr << "FATAL ERROR EmbeddedBeamContact (tag: " << this->getTag() << ") : "
            << "Beam element has zero length." << endln;
        return;
    }

    Vector initXAxis(3);
    Vector initYAxis(3);
    Vector initZAxis(3);
    crdTransf->getLocalAxes(initXAxis, initYAxis, initZAxis);
    // fill mQa
    for (int i = 0; i < 3; i++) {
        mQa(i, 0) = initXAxis(i);
        mQa(i, 1) = initYAxis(i);
        mQa(i, 2) = initZAxis(i);
    }
    // set mQb = mQa : beam column element requires zero initial twist
    // if mQa = mQb -> mchi = 0
    mQc = mQb = mQa;
    mchi = 0;

    updateShapeFuncs();


    // set interface point position
    Vector	Q1 = theNodes[0]->getCrds();		// quad node 1 coordinates
    Vector	Q2 = theNodes[1]->getCrds();		// quad node 2 coordinates
    Vector	Q3 = theNodes[2]->getCrds();		// quad node 3 coordinates
    Vector	Q4 = theNodes[3]->getCrds();		// quad node 4 coordinates
    Vector	Q5 = theNodes[4]->getCrds();		// quad node 5 coordinates
    Vector	Q6 = theNodes[5]->getCrds();		// quad node 6 coordinates
    Vector	Q7 = theNodes[6]->getCrds();		// quad node 7 coordinates
    Vector	Q8 = theNodes[7]->getCrds();		// quad node 8 coordinates

    // find interface point coordinates
    m_pos = m_Ns(0) * Q1 + m_Ns(1) * Q2 + m_Ns(2) * Q3 + m_Ns(3) * Q4 + m_Ns(4) * Q5 + m_Ns(5) * Q6 + m_Ns(6) * Q7 + m_Ns(7) * Q8;

    // alternatively can use beam nodes to calculate the location
    /*
    Vector	m_B1 = theNodes[8]->getCrds();		// beam node 1 coordinates
    Vector	m_B2 = theNodes[9]->getCrds();		// beam node 2 coordinates
    Vector	m_Ba;		// beam tangent at node 1
    Vector	m_Bb;		// beam tangent at node 2
    m_Bb = m_Ba = (m_B2 - m_B1) / m_beam_length;
    m_pos = m_Hb1 * m_B1 + m_Hb2 * m_Ba + m_Hb3 * m_B2 + m_Hb4 * m_Bb + m_beam_radius * (cos(m_beam_theta)*loc_y + sin(m_beam_theta)*loc_z);
    */

    Vector	B1 = theNodes[8]->getCrds();		// beam node 1 coordinates
    Vector	B2 = theNodes[9]->getCrds();		// beam node 2 coordinates
    Vector	Ba = Geta1();		// beam tangent at node 1
    Vector	Bb = Getb1();		// beam tangent at node 2
    m_Bcl_pos_n = m_Hb1 * B1 + m_Hb2 * Ba + m_Hb3 * B2 + m_Hb4 * Bb;

    // opserr << "element " << this->getTag() << " location = " << m_pos;

    ComputeBphiAndBu(mBphi, mBu);

    Vector c2(3), c3(3);
    for (int i = 0; i<3; i++)
    {
        c2(i) = mQc(i, 1);
        c3(i) = mQc(i, 2);
    }

    m_n = cos(m_beam_theta) * c2 + sin(m_beam_theta) * c3;
    inContact = true;

    this->DomainComponent::setDomain(theDomain);
    return;
}

int
EmbeddedBeamContact::update(void)
{
    Vector solidDisp(3), beamDisp(3);
    Vector Q1_disp(3), Q2_disp(3), Q3_disp(3), Q4_disp(3), Q5_disp(3), Q6_disp(3), Q7_disp(3), Q8_disp(3);
    Vector B1_pos(3), B2_pos(3);
    Vector B1_rot(3), B2_rot(3);
    Vector tempB1(6), tempB2(6);
    Vector Ba1(3), Bb1(3);
    Vector solidPtDisp(3), beamPtDisp(3), beamPtRot(3);

    Q1_disp = theNodes[0]->getTrialDisp();
    Q2_disp = theNodes[1]->getTrialDisp();
    Q3_disp = theNodes[2]->getTrialDisp();
    Q4_disp = theNodes[3]->getTrialDisp();
    Q5_disp = theNodes[4]->getTrialDisp();
    Q6_disp = theNodes[5]->getTrialDisp();
    Q7_disp = theNodes[6]->getTrialDisp();
    Q8_disp = theNodes[7]->getTrialDisp();

    // calculate location of the interface point in solid and beam
    solidPtDisp = m_Ns(0) * Q1_disp + m_Ns(1) * Q2_disp + m_Ns(2) * Q3_disp + m_Ns(3) * Q4_disp + m_Ns(4) * Q5_disp + m_Ns(5) * Q6_disp + m_Ns(6) * Q7_disp + m_Ns(7) * Q8_disp;

    Vector phi_c(3), u_c(3);
    Vector c2n(3);
    Vector c3n(3);
    Vector c2n1(3);
    Vector c3n1(3);
    Vector incDisp_ab(12);


    tempB1 = theNodes[8]->getTrialDisp() - theNodes[8]->getDisp();
    tempB2 = theNodes[9]->getTrialDisp() - theNodes[9]->getDisp();
    
    for (int i = 0; i<3; i++) {
        c2n(i) = mQc(i, 1);
        c3n(i) = mQc(i, 2);

        incDisp_ab(i)     = tempB1(i);
        incDisp_ab(i + 3) = tempB1(i+3);
        incDisp_ab(i + 6) = tempB2(i);
        incDisp_ab(i + 9) = tempB2(i+3);
    }

    phi_c     = mBphi * incDisp_ab;
    u_c       = mBu   * incDisp_ab;

    m_Bcl_pos = m_Bcl_pos_n + u_c;
    c2n1 = c2n + CrossProduct(phi_c, c2n);
    c3n1 = c3n + CrossProduct(phi_c, c3n);

    beamPtDisp = m_Bcl_pos + m_beam_radius * (cos(m_beam_theta) * c2n1 + sin(m_beam_theta) * c3n1);

    m_B_loc  = beamPtDisp;
    m_S_disp = solidPtDisp;

    double gap = (m_pos + m_S_disp - m_B_loc) ^ m_n;
    inContact = (gap <= 0.0);
    
    return 0;
}

int
EmbeddedBeamContact::commitState(void)
{
    int retVal = 0;

    // update Qa and Qb
    UpdateTransforms();

    m_Ba1 = Geta1();
    m_Bb1 = Getb1();

    //compute c1
    Vector temp_a_new = theNodes[8]->getCrds();
    Vector temp_b_new = theNodes[9]->getCrds();
    Vector temp_disp_a = theNodes[8]->getTrialDisp();
    Vector temp_disp_b = theNodes[9]->getTrialDisp();
    Vector temp_c1(3);

    for (int ii = 0; ii < 3; ii++)
    {
        temp_a_new(ii) += temp_disp_a(ii);
        temp_b_new(ii) += temp_disp_b(ii);
    }

    temp_c1 = m_dH1 * temp_a_new + m_dH2 * m_Ba1 + m_dH3 * temp_b_new + m_dH4 * m_Bb1;
    temp_c1 = temp_c1 / (temp_c1.Norm());
    Setc1(temp_c1);

    // update Qc
    ComputeQc();

    // compute Bphi
    ComputeBphiAndBu(mBphi, mBu);

    // update stuff for next step
    temp_a_new = theNodes[8]->getTrialDisp();
    temp_b_new = theNodes[9]->getTrialDisp();
    for (int ii = 0; ii < 3; ii++)
    {
        m_Ba_disp_n(ii) = temp_a_new(ii);
        m_Bb_disp_n(ii) = temp_b_new(ii);
        m_Ba_rot_n(ii)  = temp_a_new(ii+3);
        m_Bb_rot_n(ii)  = temp_b_new(ii+3);
    }
    
    // update the normal vector
    Vector c2(3), c3(3);
    for (int i = 0; i<3; i++)
    {
        c2(i) = mQc(i, 1);
        c3(i) = mQc(i, 2);
    }

    m_n = cos(m_beam_theta) * c2 + sin(m_beam_theta) * c3;
    
    m_Bcl_pos_n = m_Bcl_pos;

    // call element commitState to do any base class stuff
    if ((retVal = this->Element::commitState()) != 0) {
        opserr << "EmbeddedBeamContact::commitState() - failed in base class";
    }

    return retVal;
}

int EmbeddedBeamContact::updateShapeFuncs()
{
    if ((m_solid_xi < -1.0) || (m_solid_xi > 1.0) || (m_solid_eta < -1.0) || (m_solid_eta > 1.0) || (m_solid_zeta < -1.0) || (m_solid_zeta > 1.0))
    {
        opserr << "Error in shape function." << endln;
        return -1;
    }

    if ((m_beam_rho < -1.0) || (m_beam_rho > 1.0))
    {
        opserr << "Error in shape function." << endln;
        return -1;
    }

    double rho  = m_beam_rho;
    double rho2 = rho * rho;
    double rho3 = rho * rho2;

    m_Ns(0) = -0.125 * (m_solid_xi - 1) * (m_solid_eta - 1) * (m_solid_zeta - 1);
    m_Ns(1) =  0.125 * (m_solid_xi + 1) * (m_solid_eta - 1) * (m_solid_zeta - 1);
    m_Ns(2) = -0.125 * (m_solid_xi + 1) * (m_solid_eta + 1) * (m_solid_zeta - 1);
    m_Ns(3) =  0.125 * (m_solid_xi - 1) * (m_solid_eta + 1) * (m_solid_zeta - 1);
    m_Ns(4) =  0.125 * (m_solid_xi - 1) * (m_solid_eta - 1) * (m_solid_zeta + 1);
    m_Ns(5) = -0.125 * (m_solid_xi + 1) * (m_solid_eta - 1) * (m_solid_zeta + 1);
    m_Ns(6) =  0.125 * (m_solid_xi + 1) * (m_solid_eta + 1) * (m_solid_zeta + 1);
    m_Ns(7) = -0.125 * (m_solid_xi - 1) * (m_solid_eta + 1) * (m_solid_zeta + 1);

    m_Hb1 = 0.125 * (4.0 - 6.0 * rho + 2.0 * rho3);
    m_Hb3 = 0.125 * (4.0 + 6.0 * rho - 2.0 * rho3);
    m_Hb2 = 0.125 * m_beam_length * (1.0 - rho - rho2 + rho3);
    m_Hb4 = 0.125 * m_beam_length * (-1.0 - rho + rho2 + rho3);

    m_Nb1 = 0.5 * (1 - rho);
    m_Nb2 = 0.5 * (1 + rho);

    m_dH1 = 0.75 * (-1.0 + rho2);
    m_dH3 = 0.75 * (1.0 - rho2);
    m_dH2 = 0.125 * m_beam_length * (-1.0 - 2.0 * rho + 3.0 * rho2);
    m_dH4 = 0.125 * m_beam_length * (-1.0 + 2.0 * rho + 3.0 * rho2);

    return 0;
}

Vector
EmbeddedBeamContact::CrossProduct(const Vector &V1, const Vector &V2)
{
    Vector V3(3);

    V3(0) = V1(1)*V2(2) - V1(2)*V2(1);
    V3(1) = V1(2)*V2(0) - V1(0)*V2(2);
    V3(2) = V1(0)*V2(1) - V1(1)*V2(0);

    return V3;
}

Vector
EmbeddedBeamContact::Geta1(void) 
{
    Vector a1(3);
    int i;

    for (i = 0; i<3; i++) {
        a1(i) = mQa(i, 0);
    }

    return a1;
}

Vector
EmbeddedBeamContact::Getb1(void) 
{
    Vector b1(3);
    int i;

    for (i = 0; i<3; i++) {
        b1(i) = mQb(i, 0);
    }

    return b1;
}

void
EmbeddedBeamContact::Setc1(Vector c1_vec) 
{
    mc1 = c1_vec;

    return;
}

Vector
EmbeddedBeamContact::Getc1(void) {
    return mc1;
}


Matrix
EmbeddedBeamContact::Transpose(int dim1, int dim2, const Matrix &M)
{
    // copied from transpose function in Brick.cpp

    Matrix Mtran(dim2, dim1);

    for (int i = 0; i < dim1; i++)
        for (int j = 0; j < dim2; j++)
            Mtran(j, i) = M(i, j);

    return Mtran;
}


Matrix
EmbeddedBeamContact::ComputeSkew(Vector th)
{
    Matrix skew_th(3, 3);

    skew_th(0, 0) = 0.0;
    skew_th(0, 1) = -th(2);
    skew_th(0, 2) = th(1);
    skew_th(1, 0) = th(2);
    skew_th(1, 1) = 0.0;
    skew_th(1, 2) = -th(0);
    skew_th(2, 0) = -th(1);
    skew_th(2, 1) = th(0);
    skew_th(2, 2) = 0.0;

    return skew_th;
}

void
EmbeddedBeamContact::ComputeBphiAndBu(Matrix &Bphi, Matrix &Bu)
{
    int i, j;
    Matrix dummy1(3, 3);
    Matrix dummy2(3, 3);
    Matrix dummy3(3, 3);
    Matrix dummy4(3, 3);
    double L = m_beam_length / 2.0;

    Bphi.Zero();
    Bu.Zero();
    
    // Compute Bphi(0:2, 3:5)
    dummy1.Zero();
    dummy2.Zero();
    dummy3.Zero();
    dummy4.Zero();

    // dummy1 = N1 * Qc*(E1 dyadic E1)
    dummy1(0, 0) = m_Nb1*mQc(0, 0);
    dummy1(1, 0) = m_Nb1*mQc(1, 0);
    dummy1(2, 0) = m_Nb1*mQc(2, 0);
    // dummy1 += dH2 * Qc*P1
    dummy1(0, 1) = m_dH2*mQc(0, 1)/L;  // dH2 * mQc(0:2,1:2)
    dummy1(1, 1) = m_dH2*mQc(1, 1)/L;
    dummy1(2, 1) = m_dH2*mQc(2, 1)/L;
    dummy1(0, 2) = m_dH2*mQc(0, 2)/L;
    dummy1(1, 2) = m_dH2*mQc(1, 2)/L;
    dummy1(2, 2) = m_dH2*mQc(2, 2)/L;
    // dummy2 = Qa^T
    dummy2 = Transpose(3, 3, mQa);
    // dummy3 = Qc * (N1*(E1 dyadic E1)+ dH2 * P1) * Qa^T
    dummy3 = dummy1*dummy2;
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            Bphi(i, 3 + j) = dummy3(i, j);

    // Reuse parts of dummy1 and dummy2 to calculate Bu(0:2,0:2)
    // dummy1 += H1 * Qc*P1
    dummy1(0, 1) = m_Hb1*mQc(0, 1);  // H1 * mQc(0:2,1:2)
    dummy1(1, 1) = m_Hb1*mQc(1, 1);
    dummy1(2, 1) = m_Hb1*mQc(2, 1);
    dummy1(0, 2) = m_Hb1*mQc(0, 2);
    dummy1(1, 2) = m_Hb1*mQc(1, 2);
    dummy1(2, 2) = m_Hb1*mQc(2, 2);
    // dummy3 = Qc * (N1*(E1 dyadic E1)+ H1 * P1) * Qa^T
    dummy3 = dummy1*dummy2;
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            Bu(i, j) = dummy3(i, j);

    // Reuse dummy2 and Compute Bphi(0:2, 0:2) and Bu(0:2, 3:5)
    dummy1.Zero();
    dummy3.Zero();
    // dummy1 = Qc*E^R*P1 (E^R is the skew symmetric meatrix for E1 cross product => E1 x a = [E^R].a)
    dummy1(0, 1) = mQc(0, 2);
    dummy1(0, 2) = -mQc(0, 1);
    dummy1(1, 1) = mQc(1, 2);
    dummy1(1, 2) = -mQc(1, 1);
    dummy1(2, 1) = mQc(2, 2);
    dummy1(2, 2) = -mQc(2, 1);
    // dummy3 = Qc*E^R*P1*Qa^T
    dummy3 = dummy1*dummy2;
    // Compute Bphi(0:2, 0:2)
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            Bphi(i, j) = m_dH1 / L * dummy3(i, j);
    // Compute Bu(0:2, 3:5)
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            Bu(i, 3 + j) = - m_Hb2 * dummy3(i, j);

    // Reuse dummy1 and Compute Bphi(0:2, 6:8) and Bu(0:2, 9:11)
    dummy2.Zero();
    dummy3.Zero();

    // dummy2 = Qb^T
    dummy2 = Transpose(3, 3, mQb);
    // dummy3 = Qc*E^R*P1*Qb^T
    dummy3 = dummy1*dummy2;
    // Compute Bphi(0:2, 6 : 8)
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            Bphi(i, 6 + j) = m_dH3 / L * dummy3(i, j);
    // Compute Bu(0:2, 9:11)
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            Bu(i, 9 + j) = -m_Hb4 * dummy3(i, j);

    // Reuse dummy2 and Compute Bphi(0:2, 9:11)
    dummy1.Zero();
    dummy3.Zero();
    // dummy1 = N2 * Qc*(E1 dyadic E1)
    dummy1(0, 0) = m_Nb2*mQc(0, 0);  // N2 * mQc(0:2,0)
    dummy1(1, 0) = m_Nb2*mQc(1, 0);
    dummy1(2, 0) = m_Nb2*mQc(2, 0);
    // dummy1 += dH4 * Qc*P1
    dummy1(0, 1) = m_dH4*mQc(0, 1)/L;     // dH4 * mQc(0:2,1:2)
    dummy1(1, 1) = m_dH4*mQc(1, 1)/L;
    dummy1(2, 1) = m_dH4*mQc(2, 1)/L;
    dummy1(0, 2) = m_dH4*mQc(0, 2)/L;
    dummy1(1, 2) = m_dH4*mQc(1, 2)/L;
    dummy1(2, 2) = m_dH4*mQc(2, 2)/L;
    // dummy3 = Qc * (N2*(E1 dyadic E1)+ dH4 * P1) * Qb^T
    dummy3 = dummy1*dummy2;
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++) 
            Bphi(i, 9 + j) = dummy3(i, j);

    // Reuse parts of dummy1 and dummy2 to calculate Bu(0:2,6:8)
    // dummy1 += dH4 * Qc*P1
    dummy1(0, 1) = m_Hb3*mQc(0, 1);     // H3 * mQc(0:2,1:2)
    dummy1(1, 1) = m_Hb3*mQc(1, 1);
    dummy1(2, 1) = m_Hb3*mQc(2, 1);
    dummy1(0, 2) = m_Hb3*mQc(0, 2);
    dummy1(1, 2) = m_Hb3*mQc(1, 2);
    dummy1(2, 2) = m_Hb3*mQc(2, 2);
    // dummy3 = Qc * (N2*(E1 dyadic E1)+ H3 * P1) * Qb^T
    dummy3 = dummy1*dummy2;
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            Bu(i, 6 + j) = dummy3(i, j);

    return;
}

void
EmbeddedBeamContact::UpdateTransforms(void)
{
    Vector temp_a(6);           // trial disp/rot vector at a(total disp/rot)
    Vector temp_b(6);           // trial disp/rot vector at a(total disp/rot) 
    Vector rot_a(3);            // incr. rot vector at a (from n->n+1)
    Vector rot_b(3);            // incr. rot vector at a (from n->n+1)
    Matrix Omega(3, 3);         // Matrix used for Exponential Map

                                               // Recalculate incremental rotations from step n to n+1
    temp_a = theNodes[8]->getTrialDisp();
    temp_b = theNodes[9]->getTrialDisp();
    
    for (int ii = 0; ii < 3; ii++) {
        rot_a(ii)  = temp_a(ii + 3) - m_Ba_rot_n(ii);
        rot_b(ii)  = temp_b(ii + 3) - m_Bb_rot_n(ii);
    }

    // Perform exponential update of Qa
    //   calculate exponential map of current incremental rotations
    Omega = ExponentialMap(rot_a);
    //   calculate new Qa
    mQa = Omega*mQa;    

    // Perform exponential update of Qb
    //   calculate exponential map of current incremental rotations
    Omega = ExponentialMap(rot_b);
    //   calculate new Qb
    mQb = Omega*mQb;


    return;
}

void
EmbeddedBeamContact::ComputeQc()
{
    Vector c1(3);         // tangent vector at projection point, c
    Vector a1(3);         // tangent vector at a
    Vector b1(3);         // tangent vector at b
    Vector temp(3);       // dummy vector for use in calcs
    Matrix Qc_df(3, 3);   // Drill free transformation matrix for c
    Matrix Qc_chi(3, 3);  // Twist transformation matrix for c
    Matrix Qb_df(3, 3);   // Drill free transf. matrix from a to b

                                                // Fill tangent vectors
    a1 = Geta1();
    b1 = Getb1();
    c1 = Getc1();
    temp.Zero();

    // Calculate the drill free transformation from a to c, Qc_df
    temp = CrossProduct(a1, c1);
    Qc_df = ExponentialMap(temp);

    // Calculate the drill free transformation from a to b, Qb_df
    // for determination of twist angle mchi
    temp = CrossProduct(a1, b1);
    Qb_df = ExponentialMap(temp);
    Qb_df = Qb_df * mQa;

    // mchi = arcsin( b3 dot b2_df) = arcsin(mQb(:,2) dot Qb_df(:,1))
    // WATCH SIGN!!!!
    mchi = mQb(0, 2)*Qb_df(0, 1) + mQb(1, 2)*Qb_df(1, 1) + mQb(2, 2)*Qb_df(2, 1);
    mchi = -asin(mchi);

    // Calculate twist transformation from a to c, Qc_df
    // based upon linear scaling of twist angle: mxi * mchi * c1
    temp = m_Nb2*mchi*c1;
    Qc_chi = ExponentialMap(temp);

    mQc = (Qc_chi * Qc_df) * mQa;

    return;
}

Matrix
EmbeddedBeamContact::ExponentialMap(Vector th)
{
    double theta = th.Norm();             // vector norm
    Matrix sk_theta(3,3);    // skew of vector
    if (theta > 1.0e-10)
        sk_theta = ComputeSkew(th / theta);
    else
        sk_theta.Zero();
    Matrix sk_theta2 = sk_theta * sk_theta; // dyadic product of vector
    Matrix Q(3, 3);           // Exonential Map Vector  

    
    Q.Zero();
    
    Matrix meye1(3, 3);
    meye1(0, 0) = 1.0;
    meye1(1, 1) = 1.0;
    meye1(2, 2) = 1.0;

    Q = meye1 + sin(theta) * sk_theta + (1 - cos(theta)) * sk_theta2;


    return Q;
}