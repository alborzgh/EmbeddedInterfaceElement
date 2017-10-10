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
// Description: This file contains the class definition for EmbeddedEPBeamInterface.

#include <EmbeddedEPBeamInterface.h>
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

static int num_EmbeddedEPBeamInterface = 0;
static const double m_Pi = 3.14159265359;

void *
OPS_EmbeddedEPBeamInterface(void)
{
    if (num_EmbeddedEPBeamInterface == 0) {
        num_EmbeddedEPBeamInterface++;
        opserr << "EmbeddedEPBeamInterface element - Written: A.Ghofrani, D.Turello, P.Arduino, U.Washington\n";
    }

    Element *theElement = 0;

    int numArgs = OPS_GetNumRemainingInputArgs();
    if (numArgs < 1) {
        opserr << "Want: EmbeddedEPBeamInterface tag? \n";
        return 0;
    }

    int iData[1];
    int eleTag = 0;
    int numData = 1;
    if (OPS_GetIntInput(&numData, iData) != 0) {
        opserr << "WARNING invalid integer data: element EmbeddedEPBeamInterface" << endln;
        return 0;
    }

    eleTag = iData[0];

    theElement = new EmbeddedEPBeamInterface(iData[0]);

    if (theElement == 0) {
        opserr << "WARNING could not create element of type EmbeddedEPBeamInterface\n";
        return 0;
    }

    return theElement;
}


EmbeddedEPBeamInterface::EmbeddedEPBeamInterface(int tag) : Element(tag, ELE_TAG_EmbeddedEPBeamInterface)
{

}

EmbeddedEPBeamInterface::EmbeddedEPBeamInterface(int tag, int beamTag, std::vector <int> solidTag, int crdTransfTag, int matTag,
    std::vector <double>  beamRho, std::vector <double>  beamTheta, std::vector <double>  solidXi, std::vector <double>  solidEta,
    std::vector <double>  solidZeta, double radius, double area, double width) : Element(tag, ELE_TAG_EmbeddedEPBeamInterface),
    m_beam_radius(radius), m_area(area), theMatTag(matTag), m_ep(1.0e15), m_Lambda(12),
    m_Ba_rot_n(3), m_Bb_rot_n(3),
    m_Ba_disp_n(3), m_Bb_disp_n(3),
    m_Ba1(3), m_Bb1(3),
    m_Bcl_pos(3), m_Bcl_pos_n(3), m_pos(3),
    m_B_loc(3), m_S_disp(3),
    mQa(3, 3), mQb(3, 3), mQc(3, 3), mc1(3), m_uRel(12),
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
    EEPBI_numNodes = m_numSolidNodes     + 2;
    EEPBI_numDOF   = m_numSolidNodes * 3 + 12;

    // opserr << m_numSolidNodes << endln;

    externalNodes = ID(EEPBI_numNodes);
    theNodes = new Node*[EEPBI_numNodes];

    int count = 0;
    for (std::set <int>::iterator it = uniqueNodeTags.begin(); it != uniqueNodeTags.end(); ++it)
    {
        m_nodeMap[*it] = count;
        externalNodes(count) = *it;
        count++;
    }

    theElement = theDomain.getElement(beamTag);
    for (int ii = 0; ii < 2; ii++)
    {
        externalNodes(count) = theElement->getNodePtrs()[ii]->getTag();
        count++;
    }

    m_InterfaceForces    = Vector(EEPBI_numDOF);
    m_InterfaceStiffness = Matrix(EEPBI_numDOF, EEPBI_numDOF);

    mA = Matrix(3 * m_numSolidNodes, 12);
    mB = Matrix(12, 12);
    mAt = Matrix(12, 3 * m_numSolidNodes);
    mBt = Matrix(12, 12);
    mAAt = Matrix(3 * m_numSolidNodes, 3 * m_numSolidNodes);
    mBBt = Matrix(12, 12);
    mABt = Matrix(3 * m_numSolidNodes, 12);
    mBinv = Matrix(12, 12);
    mBinvTAT = Matrix(12, 3 * m_numSolidNodes);
    
    theMat = new NDMaterial*[m_numEmbeddedPoints];
    for (int ii = 0; ii < m_numEmbeddedPoints; ii++)
        theMat[ii] = OPS_getNDMaterial(matTag)->getCopy("3D");

    if (width > 0)
        m_intWidth = width;
    else
        m_intWidth = 0.1 * m_beam_radius;

    // get the coordinate transformation object
    crdTransf = OPS_GetCrdTransf(crdTransfTag)->getCopy3d();

}

EmbeddedEPBeamInterface::EmbeddedEPBeamInterface()
    : Element(0, ELE_TAG_EmbeddedEPBeamInterface)
{

}

EmbeddedEPBeamInterface::~EmbeddedEPBeamInterface()
{

}

int
EmbeddedEPBeamInterface::getNumExternalNodes(void) const
{
    return EEPBI_numNodes;
}

const ID&
EmbeddedEPBeamInterface::getExternalNodes(void)
{
    return externalNodes;
}

Node **
EmbeddedEPBeamInterface::getNodePtrs(void)
{
    return theNodes;
}

int
EmbeddedEPBeamInterface::getNumDOF(void)
{
    return EEPBI_numDOF;
}

int
EmbeddedEPBeamInterface::revertToLastCommit(void)
{
    return 0;
}

int
EmbeddedEPBeamInterface::revertToStart(void)
{
    return 0;
}


const Matrix&
EmbeddedEPBeamInterface::getTangentStiff(void)
{
    m_InterfaceStiffness.Zero();
    int Cep_map [3][3];

    Cep_map[0][0] = 0;
    Cep_map[0][1] = 3;
    Cep_map[0][2] = 5;
    Cep_map[1][0] = 3;
    Cep_map[1][1] = 1;
    Cep_map[1][2] = 4;
    Cep_map[2][0] = 5;
    Cep_map[2][1] = 4;
    Cep_map[2][2] = 2;


    Matrix Kd(12, 12);

    for (int ii = 0; ii < m_numEmbeddedPoints; ii++)
    {
        Matrix Cep(6, 6);
        Matrix temp(3, 3);

        updateShapeFuncs(m_solid_xi(ii), m_solid_eta(ii), m_solid_zeta(ii), m_beam_rho(ii));
        ComputeBphiAndBu(mBphi, mBu);

        Cep = theMat[ii]->getTangent();

        // set c1 first
        Vector c1(3), c2(3), c3(3);
        for (int jj = 0; jj < 3; jj++)
        {
            c1(jj) = mQc(jj, 0);
            c2(jj) = mQc(jj, 1);
            c3(jj) = mQc(jj, 2);
        }

        Matrix Hb(3, 12), HbT(12,3);
        Hb = mBu - (m_beam_radius*(cos(m_beam_theta(ii))*ComputeSkew(c2) + sin(m_beam_theta(ii))*ComputeSkew(c3))) * mBphi;
        HbT = Transpose(3, 12, Hb);

        Vector n(3);
        n = cos(m_beam_theta(ii)) *c2 + sin(m_beam_theta(ii))*c3;

        temp.Zero();
        for (int jj = 0; jj < 3; jj++)
            for (int kk = 0; kk < 3; kk++)
                for (int mm = 0; mm < 3; mm++)
                    for (int nn = 0; nn < 3; nn++)
                        temp(jj, kk) += Cep(Cep_map[jj][nn], Cep_map[kk][mm]) * n(mm) * n(nn);

        Kd += HbT * temp * Hb;

    }
    Kd *= m_area * m_area / m_intWidth;

    Matrix Kss(3 * m_numSolidNodes, 3 * m_numSolidNodes);
    Matrix Ksb(3 * m_numSolidNodes, 12);

    Kss = mA * mBinv * Kd * mBinvTAT;
    Ksb = -1.0 * mA * mBinv * Kd;

    for (int ii = 0; ii < 3 * m_numSolidNodes; ii++)
        for (int jj = 0; jj < 3 * m_numSolidNodes; jj++)
            m_InterfaceStiffness(ii, jj) = Kss(ii, jj);

    for (int ii = 0; ii < 3 * m_numSolidNodes; ii++)
        for (int jj = 0; jj < 12; jj++)
        {
            m_InterfaceStiffness(ii, 3 * m_numSolidNodes + jj) = Ksb(ii, jj);
            m_InterfaceStiffness(3 * m_numSolidNodes + jj, ii) = Ksb(ii, jj);
        }

    for (int ii = 0; ii < 12; ii++)
        for (int jj = 0; jj < 12; jj++)
            m_InterfaceStiffness(3 * m_numSolidNodes + ii, 3 * m_numSolidNodes + jj) = Kd(ii, jj);

    return m_InterfaceStiffness;
}

const Matrix&
EmbeddedEPBeamInterface::getInitialStiff(void)
{
    return this->getTangentStiff();
}

const Vector&
EmbeddedEPBeamInterface::getResistingForce(void)
{
    m_InterfaceForces.Zero();
    Vector temp2(12), temp(3 * m_numSolidNodes);

    temp = -1.0 * mA * mBinv * m_Lambda;
    temp2 = m_Lambda;

    for (int ii = 0; ii < 3 * m_numSolidNodes; ii++)
        m_InterfaceForces(ii) = temp(ii);
    for (int ii = 0; ii < 12; ii++)
        m_InterfaceForces(3 * m_numSolidNodes + ii) = temp2(ii);

    return m_InterfaceForces;
}

int
EmbeddedEPBeamInterface::sendSelf(int commitTag, Channel &theChannel)
{
    return 0;
}

int
EmbeddedEPBeamInterface::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker
    &theBroker)
{
    return 0;
}

int
EmbeddedEPBeamInterface::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
{
    return 0;
}

void
EmbeddedEPBeamInterface::Print(OPS_Stream &s, int flag)
{
    return;
}

Response*
EmbeddedEPBeamInterface::setResponse(const char **argv, int argc,
    OPS_Stream &s)
{
    if (strcmp(argv[0], "force") == 0 || strcmp(argv[0], "globalForce") == 0) {
        return new ElementResponse(this, 1, Vector(3 * m_numEmbeddedPoints));

    }
    else if (strcmp(argv[0], "displacement") == 0 || strcmp(argv[0], "disp") == 0)
    {
        return new ElementResponse(this, 2, Vector(3 * m_numEmbeddedPoints));

    }
    else if (strcmp(argv[0], "beamCL") == 0 || strcmp(argv[0], "beamCenterLine") == 0) {
        return new ElementResponse(this, 3, Vector(3 * m_numEmbeddedPoints));

    }
    else if (strcmp(argv[0], "c1") == 0 || strcmp(argv[0], "tangent") == 0) {
        return new ElementResponse(this, 4, Vector(3 * m_numEmbeddedPoints));

    }
    else if (strcmp(argv[0], "c2") == 0 || strcmp(argv[0], "perp2") == 0) {
        return new ElementResponse(this, 5, Vector(3 * m_numEmbeddedPoints));

    }
    else if (strcmp(argv[0], "c3") == 0 || strcmp(argv[0], "perp3") == 0) {
        return new ElementResponse(this, 6, Vector(3 * m_numEmbeddedPoints));

    }
    else {
        opserr << "EmbeddedEPBeamInterface Recorder, " << argv[0] << "is an unknown recorder request"
            << "  Element tag : " << this->getTag() << endln;
        return 0;
    }
}

int
EmbeddedEPBeamInterface::getResponse(int responseID, Information &eleInformation)
{
    // TODO: Needs to be updated
    if (responseID == 1) // force
        return eleInformation.setVector(GetInteractionPtForce());
    else if (responseID == 2) // displacement
        return eleInformation.setVector(GetInteractionPtDisp());
    else if (responseID == 3) // centerline
        return eleInformation.setVector(m_Bcl_pos);
    else if (responseID == 4) // c1
    { 
        Vector temp(3);
        for (int ii = 0; ii < 3; ii++)
            temp(ii) = mQc(ii, 0);
        return eleInformation.setVector(temp);
    }
    else if (responseID == 5) // c2
    { 
        Vector temp(3);
        for (int ii = 0; ii < 3; ii++)
            temp(ii) = mQc(ii, 1);
        return eleInformation.setVector(temp);
    }
    else if (responseID == 6) // c3
    { 
        Vector temp(3);
        for (int ii = 0; ii < 3; ii++)
            temp(ii) = mQc(ii, 2);
        return eleInformation.setVector(temp);
    }
    else {
        opserr << "EmbeddedEPBeamInterface, tag = " << this->getTag()
            << " -- unknown request" << endln;
        return -1;
    }
}

int
EmbeddedEPBeamInterface::setParameter(const char **argv, int argc, Parameter &param)
{
    return 0;
}

int
EmbeddedEPBeamInterface::updateParameter(int parameterID, Information &info)
{
    return 0;
}

void
EmbeddedEPBeamInterface::setDomain(Domain *theDomain)
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
        opserr << "EmbeddedBeamInterfaceAL2::setDomain(): Error initializing coordinate transformation";
        return;
    }

    m_beam_length = crdTransf->getInitialLength();
    if (m_beam_length < 1.0e-12)
    {
        opserr << "FATAL ERROR EmbeddedBeamInterfaceAL2 (tag: " << this->getTag() << ") : "
            << "Beam element has zero length." << endln;
        return;
    }

    Vector initXAxis(3);
    Vector initYAxis(3);
    Vector initZAxis(3);
    crdTransf->getLocalAxes(initXAxis, initYAxis, initZAxis);
    // fill mQa
    for (int i = 0; i < 3; i++)
    {
        mQa(i, 0) = initXAxis(i);
        mQa(i, 1) = initYAxis(i);
        mQa(i, 2) = initZAxis(i);
    }
    // set mQb = mQa : beam column element requires zero initial twist
    // if mQa = mQb -> mchi = 0
    mQc = mQb = mQa;
    mchi = 0;


    // calculate A and B
    mA.Zero();
    mB.Zero();
    for (int ii = 0; ii < m_numEmbeddedPoints; ii++)
    {
        updateShapeFuncs(m_solid_xi(ii), m_solid_eta(ii), m_solid_zeta(ii), m_beam_rho(ii));
        mHf.Zero();
        ComputeHf(mHf, m_beam_theta(ii));

        Element * theElement = theDomain->getElement(theSolidTag[ii]);
        double oneOver2PiR = 0.5 / m_Pi / m_beam_radius * m_area;
        double oneOver2PiR2 = oneOver2PiR / m_beam_radius;

        for (int jj = 0; jj < 8; jj++)
        {
            int nodeInA = m_nodeMap[theElement->getNodePtrs()[jj]->getTag()];

            // opserr << "Element " << theSolidTag[ii] << " - Node " << theElement->getNodePtrs()[jj]->getTag() << " which is " << nodeInA << " in the local Mat." << endln;

            for (int kk = 0; kk < 12; kk++)
            {
                mA(3 * nodeInA, kk) += m_Ns(jj) * mHf(0, kk);
                mA(3 * nodeInA + 1, kk) += m_Ns(jj) * mHf(1, kk);
                mA(3 * nodeInA + 2, kk) += m_Ns(jj) * mHf(2, kk);
            }
        }

        ComputeBphiAndBu(mBphi, mBu);

        // I need to update Qc as well!

        Vector c2(3), c3(3);
        for (int ii = 0; ii < 3; ii++)
        {
            c2(ii) = mQc(ii, 1);
            c3(ii) = mQc(ii, 2);
        }

        Matrix Hb(3, 12);
        Hb = mBu - (m_beam_radius*(cos(m_beam_theta(ii))*ComputeSkew(c2) + sin(m_beam_theta(ii))*ComputeSkew(c3))) * mBphi;

        Matrix HbT = Transpose(3, 12, Hb);
        mB += HbT * mHf;
    }
    mA *= m_area;
    mB *= m_area;
    mAt = Transpose(3 * m_numSolidNodes, 12, mA);
    mBt = Transpose(12, 12, mB);
    mAAt = mA * mAt;
    mBBt = mB * mBt;
    mABt = mA * mBt;
    mB.Invert(mBinv);
    mBinvTAT = Transpose(12, 12, mBinv) * mAt;

    this->DomainComponent::setDomain(theDomain);
    return;
}

int
EmbeddedEPBeamInterface::update(void)
{
    Vector sDisp(3 * m_numSolidNodes), bDisp(12);
    for (int ii = 0; ii < m_numSolidNodes; ii++)
    {
        sDisp(3 * ii) = theNodes[ii]->getTrialDisp()(0);
        sDisp(3 * ii + 1) = theNodes[ii]->getTrialDisp()(1);
        sDisp(3 * ii + 2) = theNodes[ii]->getTrialDisp()(2);
    }

    for (int jj = 0; jj < 6; jj++)
    {
        bDisp(jj) = theNodes[m_numSolidNodes]->getTrialDisp()(jj);
        bDisp(jj + 6) = theNodes[m_numSolidNodes + 1]->getTrialDisp()(jj);
    }

    m_Lambda.Zero();
 
    m_uRel = bDisp - mBinvTAT * sDisp;

    for (int ii = 0; ii < m_numEmbeddedPoints; ii++) 
    {

        updateShapeFuncs(m_solid_xi(ii), m_solid_eta(ii), m_solid_zeta(ii), m_beam_rho(ii));
        ComputeBphiAndBu(mBphi, mBu);

        // set c1 first
        Vector c1(3), c2(3), c3(3);
        for (int jj = 0; jj < 3; jj++)
        {
            c1(jj) = mQc(jj, 0);
            c2(jj) = mQc(jj, 1);
            c3(jj) = mQc(jj, 2);
        }
        
        Matrix Hb(3, 12), HbT(12, 3);
        Hb  = mBu - (m_beam_radius*(cos(m_beam_theta(ii))*ComputeSkew(c2) + sin(m_beam_theta(ii))*ComputeSkew(c3))) * mBphi;
        HbT = Transpose(3, 12, Hb);

        Vector loc_relDisp(3);
        loc_relDisp = Hb * m_uRel;

        Vector n(3);
        n = cos(m_beam_theta(ii)) *c2 + sin(m_beam_theta(ii))*c3;

        Vector eps(6), sig(6);
        eps(0) = loc_relDisp(0) * n(0);
        eps(1) = loc_relDisp(1) * n(1);
        eps(2) = loc_relDisp(2) * n(2);
        eps(3) = loc_relDisp(0) * n(1) + loc_relDisp(1) * n(0);
        eps(4) = loc_relDisp(1) * n(2) + loc_relDisp(2) * n(1);
        eps(5) = loc_relDisp(0) * n(2) + loc_relDisp(2) * n(0);

        eps *= 1/m_intWidth;

        theMat[ii]->setTrialStrain(eps);
        sig = theMat[ii]->getStress();

        Vector sig_n(3);
        sig_n(0) = sig(0) * n(0) + sig(3) * n(1) + sig(5) * n(2);
        sig_n(1) = sig(3) * n(0) + sig(1) * n(1) + sig(4) * n(2);
        sig_n(2) = sig(5) * n(0) + sig(4) * n(1) + sig(2) * n(2);

        m_Lambda += HbT * sig_n;

    }
    m_Lambda *= m_area * m_area;


    return 0;
}

int
EmbeddedEPBeamInterface::commitState(void)
{
    int retVal = 0;

    // call element commitState to do any base class stuff
    if ((retVal = this->Element::commitState()) != 0) {
        opserr << "EmbeddedEPBeamInterface::commitState() - failed in base class";
    }

    return retVal;
}

int EmbeddedEPBeamInterface::updateShapeFuncs(double xi, double eta, double zeta, double rho)
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
EmbeddedEPBeamInterface::CrossProduct(const Vector &V1, const Vector &V2)
{
    Vector V3(3);

    V3(0) = V1(1)*V2(2) - V1(2)*V2(1);
    V3(1) = V1(2)*V2(0) - V1(0)*V2(2);
    V3(2) = V1(0)*V2(1) - V1(1)*V2(0);

    return V3;
}

Vector
EmbeddedEPBeamInterface::Geta1(void) 
{
    Vector a1(3);
    int i;

    for (i = 0; i<3; i++) {
        a1(i) = mQa(i, 0);
    }

    return a1;
}

Vector
EmbeddedEPBeamInterface::Getb1(void) 
{
    Vector b1(3);
    int i;

    for (i = 0; i<3; i++) {
        b1(i) = mQb(i, 0);
    }

    return b1;
}

void
EmbeddedEPBeamInterface::Setc1(Vector c1_vec) 
{
    mc1 = c1_vec;

    return;
}

Vector
EmbeddedEPBeamInterface::Getc1(void) {
    return mc1;
}


Matrix
EmbeddedEPBeamInterface::Transpose(int dim1, int dim2, const Matrix &M)
{
    // copied from transpose function in Brick.cpp

    Matrix Mtran(dim2, dim1);

    for (int i = 0; i < dim1; i++)
        for (int j = 0; j < dim2; j++)
            Mtran(j, i) = M(i, j);

    return Mtran;
}


Matrix
EmbeddedEPBeamInterface::ComputeSkew(Vector th)
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
EmbeddedEPBeamInterface::ComputeBphiAndBu(Matrix &Bphi, Matrix &Bu)
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

void EmbeddedEPBeamInterface::ComputeHf(Matrix & Hf, double theta)
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
EmbeddedEPBeamInterface::UpdateTransforms(void)
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
EmbeddedEPBeamInterface::ComputeQc()
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
EmbeddedEPBeamInterface::ExponentialMap(Vector th)
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

Vector EmbeddedEPBeamInterface::GetInteractionPtDisp()
{
    Vector res(3 * m_numEmbeddedPoints);
    Vector c2(3), c3(3);
    Vector bDisp(12);
    Matrix Hb(3, 12);
    Vector ptDisp(3);

    // update local coordinate system
    for (int ii = 0; ii < 3; ii++)
    {
        c2(ii) = mQc(ii, 1);
        c3(ii) = mQc(ii, 2);
    }

    for (int ii = 0; ii < 6; ii++)
    {
        bDisp(ii) = theNodes[m_numSolidNodes]->getTrialDisp()(ii);
        bDisp(ii + 6) = theNodes[m_numSolidNodes + 1]->getTrialDisp()(ii);
    }


    for (int ii = 0; ii < m_numEmbeddedPoints; ii++)
    {
        // update the interpolation functions for displacements and interaction forces
        updateShapeFuncs(m_solid_xi(ii), m_solid_eta(ii), m_solid_zeta(ii), m_beam_rho(ii));

        // update the matrices that define kinematics of the points on the beam surface
        ComputeBphiAndBu(mBphi, mBu);

        Hb = mBu - (m_beam_radius*(cos(m_beam_theta(ii))*ComputeSkew(c2) + sin(m_beam_theta(ii))*ComputeSkew(c3))) * mBphi;
        ptDisp = Hb * bDisp;

        for (int jj = 0; jj < 3; jj++)
            res(3 * ii + jj) = ptDisp(jj);
    }

    return res;
}

Vector EmbeddedEPBeamInterface::GetInteractionPtForce()
{
    Vector res(3 * m_numEmbeddedPoints);
    Vector ptForces(3);
    Matrix Hf(3, 12);


    for (int ii = 0; ii < m_numEmbeddedPoints; ii++)
    {
        // update the interpolation functions for displacements and interaction forces
        updateShapeFuncs(m_solid_xi(ii), m_solid_eta(ii), m_solid_zeta(ii), m_beam_rho(ii));

        // update the matrices that define kinematics of the points on the beam surface
        ComputeHf(Hf, m_beam_theta(ii));

        ptForces = Hf * m_Lambda;

        for (int jj = 0; jj < 3; jj++)
            res(3 * ii + jj) = ptForces(jj);
    }

    return res;
}