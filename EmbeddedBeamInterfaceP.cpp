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
// Description: This file contains the class definition for EmbeddedBeamInterfaceP.

#include <EmbeddedBeamInterfaceP.h>
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
#include <NodeIter.h>

static int num_EmbeddedBeamInterfaceP = 0;
static const double m_Pi = 3.14159265359;

void *
OPS_EmbeddedBeamInterfaceP(void)
{
    if (num_EmbeddedBeamInterfaceP == 0) {
        num_EmbeddedBeamInterfaceP++;
        opserr << "EmbeddedBeamInterfaceP element - Written: A.Ghofrani, D.Turello, P.Arduino, U.Washington\n";
    }

    Element *theElement = 0;

    int numArgs = OPS_GetNumRemainingInputArgs();
    if (numArgs < 1) {
        opserr << "Want: EmbeddedBeamInterfaceP tag? \n";
        return 0;
    }

    int iData[1];
    int eleTag = 0;
    int numData = 1;
    if (OPS_GetIntInput(&numData, iData) != 0) {
        opserr << "WARNING invalid integer data: element EmbeddedBeamInterfaceP" << endln;
        return 0;
    }

    eleTag = iData[0];

    theElement = new EmbeddedBeamInterfaceP(iData[0]);

    if (theElement == 0) {
        opserr << "WARNING could not create element of type EmbeddedBeamInterfaceP\n";
        return 0;
    }

    return theElement;
}


EmbeddedBeamInterfaceP::EmbeddedBeamInterfaceP(int tag) : Element(tag, ELE_TAG_EmbeddedBeamInterfaceP)
{

}

EmbeddedBeamInterfaceP::EmbeddedBeamInterfaceP(int tag, int beamTag, std::vector <int> solidTag, int crdTransfTag, 
    std::vector <double>  beamRho, std::vector <double>  beamTheta, std::vector <double>  solidXi, std::vector <double>  solidEta,
    std::vector <double>  solidZeta, double radius, double area): Element(tag, ELE_TAG_EmbeddedBeamInterfaceP),
    m_beam_radius(radius), m_area(area), m_ep(1.0e12),
    m_Ba_rot_n(3), m_Bb_rot_n(3),
    m_Ba_disp_n(3), m_Bb_disp_n(3),
    m_Ba1(3), m_Bb1(3),
    m_Bcl_pos(3), m_Bcl_pos_n(3), m_pos(3),
    m_B_loc(3), m_S_disp(3),
    mQa(3, 3), mQb(3, 3), mQc(3, 3), mc1(3),
    mBphi(3, 12), mBu(3, 12), mHf(3, 12), m_Ns(8)
{
    // get domain to access element tags and their nodes
#ifdef _PARALLEL_PROCESSING
#include <PartitionedDomain.h>
    extern PartitionedDomain theDomain;
#else
    extern Domain theDomain;
#endif

    m_numEmbeddedPoints = solidTag.size();
    theSolidTag   = new int[m_numEmbeddedPoints];
    solidNodeTags = new int[8 * m_numEmbeddedPoints];
    m_beam_rho = m_beam_theta = m_solid_xi = m_solid_eta = m_solid_zeta = Vector(m_numEmbeddedPoints);

    std::set <int> uniqueNodeTags;
    Element *theElement;
    for (int ii = 0; ii < m_numEmbeddedPoints; ii++)
    {
        theSolidTag[ii]     = solidTag[ii];
        m_solid_xi(ii)      = solidXi[ii];
        m_solid_eta(ii)     = solidEta[ii];
        m_solid_zeta(ii)    = solidZeta[ii];
        m_beam_rho(ii)      = beamRho[ii];
        m_beam_theta(ii)    = beamTheta[ii];

        theElement = theDomain.getElement(solidTag[ii]);
        // opserr << "Point " << ii +1 << " : element " << solidTag[ii] << " at (" << solidXi[ii] << "," << solidEta[ii] << "," << solidZeta[ii] << ") , beam: " << beamTag << " at (" << beamRho[ii] << "," << beamTheta[ii] << ")" << endln;
        for (int jj = 0; jj < 8; jj++)
        {
            uniqueNodeTags.insert(theElement->getNodePtrs()[jj]->getTag());
            solidNodeTags[ii * 8 + jj] = theElement->getNodePtrs()[jj]->getTag();
        }
    }

    // opserr << m_beam_rho << m_beam_theta;
    // for(int ii = 0 ; ii < 8*m_numEmbeddedPoints; ii++)
    //     opserr << *(solidNodeTags+ii) << " ";
    // opserr << endln;

    m_numSolidNodes = (int)uniqueNodeTags.size();
    EBIP_numNodes = m_numSolidNodes     + 2;
    EBIP_numDOF   = m_numSolidNodes * 3 + 12;

    // opserr << m_numSolidNodes << endln;

    externalNodes = ID(EBIP_numNodes);
    theNodes = new Node*[EBIP_numNodes];

    int count = 0;
    for (std::set <int>::iterator it = uniqueNodeTags.begin(); it != uniqueNodeTags.end(); ++it)
    {
        m_nodeMap[*it] = count;
        externalNodes(count) = *it;

        // opserr << count << " = " << "map(" << *it << ")\n";

        theNodes[count] = theDomain.getNode(*it);
        count++;
    }

    theElement = theDomain.getElement(beamTag);
    for (int ii = 0; ii < 2; ii++)
    {
        Node * thisNode = theElement->getNodePtrs()[ii];
        theNodes[count] = thisNode;
        externalNodes(count) = thisNode->getTag();
        count++;
    }

    m_InterfaceForces = Vector(EBIP_numDOF);
    m_InterfaceStiffness = Matrix(EBIP_numDOF, EBIP_numDOF);

    mN = Matrix(3, 3*m_numSolidNodes);
    mH = Matrix(3, 12);

    mA = Matrix(3*m_numSolidNodes, 12);
    mB = Matrix(12, 12);

    mKss = Matrix(3 * m_numSolidNodes, 3 * m_numSolidNodes);
    mKsb = Matrix(3 * m_numSolidNodes, 12);
    mKbb = Matrix(12, 12);
    


    // get the coordinate transformation object
    crdTransf = OPS_GetCrdTransf(crdTransfTag)->getCopy3d();

}

EmbeddedBeamInterfaceP::EmbeddedBeamInterfaceP()
    : Element(0, ELE_TAG_EmbeddedBeamInterfaceP)
{

}

EmbeddedBeamInterfaceP::~EmbeddedBeamInterfaceP()
{

}

int
EmbeddedBeamInterfaceP::getNumExternalNodes(void) const
{
    return EBIP_numNodes;
}

const ID&
EmbeddedBeamInterfaceP::getExternalNodes(void)
{
    return externalNodes;
}

Node **
EmbeddedBeamInterfaceP::getNodePtrs(void)
{
    return theNodes;
}

int
EmbeddedBeamInterfaceP::getNumDOF(void)
{
    return EBIP_numDOF;
}

int
EmbeddedBeamInterfaceP::revertToLastCommit(void)
{
    return 0;
}

int
EmbeddedBeamInterfaceP::revertToStart(void)
{
    return 0;
}


const Matrix&
EmbeddedBeamInterfaceP::getTangentStiff(void)
{
    m_InterfaceStiffness.Zero();


    for (int ii = 0; ii < 3 * m_numSolidNodes; ii++)
    {
        for (int jj = 0; jj < 3 * m_numSolidNodes; jj++)
            m_InterfaceStiffness(ii, jj) = mKss(ii, jj);

        for (int jj = 0; jj < 12; jj++)
        {
            m_InterfaceStiffness(ii, 3 * m_numSolidNodes + jj) = mKsb(ii, jj);
            m_InterfaceStiffness(3 * m_numSolidNodes + jj, ii) = mKsb(ii, jj);
        }
    }

    for (int ii = 0; ii < 12; ii++)
        for (int jj = 0; jj < 12; jj++)
            m_InterfaceStiffness(3 * m_numSolidNodes + ii, 3 * m_numSolidNodes + jj) = mKbb(ii, jj);

    m_InterfaceStiffness *= m_ep;

    return m_InterfaceStiffness;
}

const Matrix&
EmbeddedBeamInterfaceP::getInitialStiff(void)
{
    return this->getTangentStiff();
}

const Vector&
EmbeddedBeamInterfaceP::getResistingForce(void)
{
    m_InterfaceForces.Zero();
    Vector temp1(3 * m_numSolidNodes + 12);

    for (int ii = 0; ii < m_numSolidNodes; ii++)
    {
        temp1(3 * ii)     = theNodes[ii]->getTrialDisp()(0);
        temp1(3 * ii + 1) = theNodes[ii]->getTrialDisp()(1);
        temp1(3 * ii + 2) = theNodes[ii]->getTrialDisp()(2);
    }

    for (int ii = 0; ii < 6; ii++)
    {
        temp1(3 * m_numSolidNodes + ii)     = theNodes[m_numSolidNodes]->getTrialDisp()(ii);
        temp1(3 * m_numSolidNodes + ii + 6) = theNodes[m_numSolidNodes + 1]->getTrialDisp()(ii);
    }

    m_InterfaceForces = this->getTangentStiff() * temp1;

    return m_InterfaceForces;
}

int
EmbeddedBeamInterfaceP::sendSelf(int commitTag, Channel &theChannel)
{
    return 0;
}

int
EmbeddedBeamInterfaceP::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker
    &theBroker)
{
    return 0;
}

int
EmbeddedBeamInterfaceP::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
{
    return 0;
}

void
EmbeddedBeamInterfaceP::Print(OPS_Stream &s, int flag)
{
    return;
}

Response*
EmbeddedBeamInterfaceP::setResponse(const char **argv, int argc,
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
        opserr << "EmbeddedBeamInterfaceP Recorder, " << argv[0] << "is an unknown recorder request"
            << "  Element tag : " << this->getTag() << endln;
        return 0;
    }
}

int
EmbeddedBeamInterfaceP::getResponse(int responseID, Information &eleInformation)
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
        opserr << "EmbeddedBeamInterfaceP, tag = " << this->getTag()
            << " -- unknown request" << endln;
        return -1;
    }
}

int
EmbeddedBeamInterfaceP::setParameter(const char **argv, int argc, Parameter &param)
{
    return 0;
}

int
EmbeddedBeamInterfaceP::updateParameter(int parameterID, Information &info)
{
    return 0;
}

void
EmbeddedBeamInterfaceP::setDomain(Domain *theDomain)
{
    for (int ii = 0; ii < m_numSolidNodes + 2; ii++)
    {
        theNodes[ii] = theDomain->getNode(externalNodes(ii));
        if (theNodes[ii] == 0)
        {
            opserr << "Could not find node " << externalNodes(ii) << "." << endln;
            return;
        }
        if ((theNodes[ii]->getNumberDOF() != 3) && (ii < m_numSolidNodes))
        {
            opserr << "Solid node " << externalNodes(ii) << " has to have 3 degrees of freedom." << endln;
            return;
        }
        if ((theNodes[ii]->getNumberDOF() != 6) && (ii > m_numSolidNodes - 1))
        {
            opserr << "Beam node " << externalNodes(ii) << " has to have 6 degrees of freedom." << endln;
            return;
        }
    }
    
    // initialize the transformation
    if (crdTransf->initialize(theNodes[m_numSolidNodes], theNodes[m_numSolidNodes + 1]))
    {
        opserr << "EmbeddedBeamInterfaceP::setDomain(): Error initializing coordinate transformation";
        return;
    }
    
    m_beam_length = crdTransf->getInitialLength();
    if (m_beam_length < 1.0e-12) {
        opserr << "FATAL ERROR EmbeddedBeamInterfaceP (tag: " << this->getTag() << ") : "
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
    

    // calculate A and B
    mH.Zero();
    mN.Zero();
    mKss.Zero();
    mKsb.Zero();
    mKbb.Zero();
    for (int ii = 0; ii < m_numEmbeddedPoints; ii++)
    {
        updateShapeFuncs(m_solid_xi(ii), m_solid_eta(ii), m_solid_zeta(ii), m_beam_rho(ii));

        Element * theElement = theDomain->getElement(theSolidTag[ii]);
        for (int jj = 0; jj < 8; jj++)
        {
            int nodeInK1 = m_nodeMap[theElement->getNodePtrs()[jj]->getTag()];
            mN(0, 3 * nodeInK1)     += m_Ns(jj);
            mN(1, 3 * nodeInK1 + 1) += m_Ns(jj);
            mN(2, 3 * nodeInK1 + 2) += m_Ns(jj);

        }
                
        //opserr << mA << mC;

        ComputeBphiAndBu(mBphi, mBu);

        // I need to update Qc as well!

        Vector c2(3), c3(3);
        for (int ii = 0; ii < 3; ii++)
        {
            c2(ii) = mQc(ii, 1);
            c3(ii) = mQc(ii, 2);
        }

        Matrix Hb(3, 12);
        Hb  = mBu - (m_beam_radius*(cos(m_beam_theta(ii))*ComputeSkew(c2) + sin(m_beam_theta(ii))*ComputeSkew(c3))) * mBphi;
        mH += Hb;
    }

    mN *= m_area;
    mH *= m_area;

    mKss = Transpose(3, 3 * m_numSolidNodes, mN) * mN;
    mKsb = -1.0 * Transpose(3, 3 * m_numSolidNodes, mN) * mH;
    mKbb = Transpose(3, 12, mH) * mH;

    // opserr << "mA = " << mA;
    // opserr << "mB = " << mB;

    this->DomainComponent::setDomain(theDomain);
    return;
}

int
EmbeddedBeamInterfaceP::update(void)
{

    return 0;
}

int
EmbeddedBeamInterfaceP::commitState(void)
{
    int retVal = 0;

    // call element commitState to do any base class stuff
    if ((retVal = this->Element::commitState()) != 0) {
        opserr << "EmbeddedBeamInterfaceP::commitState() - failed in base class";
    }

    return retVal;
}

int EmbeddedBeamInterfaceP::updateShapeFuncs(double xi, double eta, double zeta, double rho)
{
    if ((xi < -1.0) || (xi > 1.0) || (eta < -1.0) || (eta > 1.0) || (zeta < -1.0) || (zeta > 1.0))
    {
        opserr << "Error in shape function." << endln;
        return -1;
    }

    if ((rho < -1.0) || (rho > 1.0))
    {
        opserr << "Error in shape function." << endln;
        return -1;
    }

    double rho2 = rho * rho;
    double rho3 = rho * rho2;

    m_Ns(0) = -0.125 * (xi - 1) * (eta - 1) * (zeta - 1);
    m_Ns(1) =  0.125 * (xi + 1) * (eta - 1) * (zeta - 1);
    m_Ns(2) = -0.125 * (xi + 1) * (eta + 1) * (zeta - 1);
    m_Ns(3) =  0.125 * (xi - 1) * (eta + 1) * (zeta - 1);
    m_Ns(4) =  0.125 * (xi - 1) * (eta - 1) * (zeta + 1);
    m_Ns(5) = -0.125 * (xi + 1) * (eta - 1) * (zeta + 1);
    m_Ns(6) =  0.125 * (xi + 1) * (eta + 1) * (zeta + 1);
    m_Ns(7) = -0.125 * (xi - 1) * (eta + 1) * (zeta + 1);

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
EmbeddedBeamInterfaceP::CrossProduct(const Vector &V1, const Vector &V2)
{
    Vector V3(3);

    V3(0) = V1(1)*V2(2) - V1(2)*V2(1);
    V3(1) = V1(2)*V2(0) - V1(0)*V2(2);
    V3(2) = V1(0)*V2(1) - V1(1)*V2(0);

    return V3;
}

Vector
EmbeddedBeamInterfaceP::Geta1(void) 
{
    Vector a1(3);
    int i;

    for (i = 0; i<3; i++) {
        a1(i) = mQa(i, 0);
    }

    return a1;
}

Vector
EmbeddedBeamInterfaceP::Getb1(void) 
{
    Vector b1(3);
    int i;

    for (i = 0; i<3; i++) {
        b1(i) = mQb(i, 0);
    }

    return b1;
}

void
EmbeddedBeamInterfaceP::Setc1(Vector c1_vec) 
{
    mc1 = c1_vec;

    return;
}

Vector
EmbeddedBeamInterfaceP::Getc1(void) {
    return mc1;
}


Matrix
EmbeddedBeamInterfaceP::Transpose(int dim1, int dim2, const Matrix &M)
{
    // copied from transpose function in Brick.cpp

    Matrix Mtran(dim2, dim1);

    for (int i = 0; i < dim1; i++)
        for (int j = 0; j < dim2; j++)
            Mtran(j, i) = M(i, j);

    return Mtran;
}


Matrix
EmbeddedBeamInterfaceP::ComputeSkew(Vector th)
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
EmbeddedBeamInterfaceP::ComputeBphiAndBu(Matrix &Bphi, Matrix &Bu)
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
    dummy1(0, 1) =  mQc(0, 2);
    dummy1(0, 2) = -mQc(0, 1);
    dummy1(1, 1) =  mQc(1, 2);
    dummy1(1, 2) = -mQc(1, 1);
    dummy1(2, 1) =  mQc(2, 2);
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
            Bu(i, 3 + j) = -m_Hb2 * dummy3(i, j);

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

void EmbeddedBeamInterfaceP::ComputeHf(Matrix & Hf, double theta)
{
    Hf.Zero();
    double oneOver2PiR = 0.5 / m_Pi / m_beam_radius * m_area;
    double oneOver2PiR2 = oneOver2PiR / m_beam_radius;
    for (int ii = 0; ii < 3; ii++)
    {
        for (int jj = 0; jj < 3; jj++)
        {
            Hf(ii, jj)     = oneOver2PiR * m_Nb1 * mQa(jj, ii);
            Hf(ii, jj + 6) = oneOver2PiR * m_Nb2 * mQb(jj, ii);
        }
        Hf(0, ii + 3)  = 2.0 * oneOver2PiR2 * m_Nb1 * (mQa(ii, 1) * sin(theta) - mQa(ii, 2) * cos(theta));
        Hf(1, ii + 3)  = -oneOver2PiR2 * mQa(ii, 0) * m_Nb1 * sin(theta);
        Hf(2, ii + 3)  =  oneOver2PiR2 * mQa(ii, 0) * m_Nb1 * cos(theta);
        Hf(0, ii + 9)  = 2.0 * oneOver2PiR2 * m_Nb2 * (mQb(ii, 1) * sin(theta) - mQb(ii, 2) * cos(theta));
        Hf(1, ii + 9)  = -oneOver2PiR2 * mQb(ii, 0) * m_Nb2 * sin(theta);
        Hf(2, ii + 9)  =  oneOver2PiR2 * mQb(ii, 0) * m_Nb2 * cos(theta);
    }
    
    Hf = mQc * Hf;
    return;
}

void
EmbeddedBeamInterfaceP::UpdateTransforms(void)
{
    Vector temp_a(6);           // trial disp/rot vector at a(total disp/rot)
    Vector temp_b(6);           // trial disp/rot vector at a(total disp/rot) 
    Vector rot_a(3);            // incr. rot vector at a (from n->n+1)
    Vector rot_b(3);            // incr. rot vector at a (from n->n+1)
    Matrix Omega(3, 3);         // Matrix used for Exponential Map

                                               // Recalculate incremental rotations from step n to n+1
    temp_a = theNodes[m_numSolidNodes]->getTrialDisp();
    temp_b = theNodes[m_numSolidNodes + 1]->getTrialDisp();
    
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
EmbeddedBeamInterfaceP::ComputeQc()
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
EmbeddedBeamInterfaceP::ExponentialMap(Vector th)
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