/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "v2ftpFinal.H"
#include "addToRunTimeSelectionTable.H"
#include "backwardsCompatibilityWallFunctions.H"
#include "components.H"
#include "fvCFD.H"
#include "volFields.H"
#include "symmTensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(v2ftpFinal, 0);
addToRunTimeSelectionTable(RASModel, v2ftpFinal, dictionary);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

tmp<volScalarField> v2ftpFinal::Ts() const
{ 
	if(tslimiter_ == "true")
	{
        return max(k_/(epsilon_ + epsilonSmall_), 6.0*sqrt(nu()/(epsilon_ + epsilonSmall_)));
	}
	
    return ((k_+k0_)/(epsilon_ + epsilonSmall_));
}


tmp<volScalarField> v2ftpFinal::Ls() const
{
	
	volScalarField trueL = pow(k_+k0_, 1.5)/(epsilon_ + epsilonSmall_);
	
	if(lslimiter_ == "true")
	{
		return cL1_*max(pow(k_+k0_, 1.5)/(epsilon_ + epsilonSmall_),cL2_*pow(pow3(nu())/(epsilon_ + epsilonSmall_),0.25));
	}
	
	//Info << "Max trueL: " << gMax(trueL) << " Min trueL: " << gMin(trueL) << endl;
	//Info << "Dims: " << trueL.dimensions() << endl;
	
	return pow(k_+k0_, 1.5)/(epsilon_ + epsilonSmall_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

v2ftpFinal::v2ftpFinal 
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& lamTransportModel
)
:
    RASModel(typeName, U, phi, lamTransportModel),


    cEp1_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cEp1",
            coeffDict_,
            1.44
        )
    ),
    cEp2con_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cEp2con",
            coeffDict_,
            1.83
        )
    ),
    cEp3_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cEp3",
            coeffDict_,
            0.15
        )
    ),
    cEp4_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cEp4",
            coeffDict_,
            0.04
        )
    ),
    cD1_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cD1",
       	    coeffDict_,
            0.0
        )
    ),
    cD2_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cD2",
       	    coeffDict_,
            1.0
        )
    ),
    cG_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cG",
       	    coeffDict_,
            0.67
        )
    ),
    cGw_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cGw",
       	    coeffDict_,
            0.0
        )
    ),
    cP1_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cP1",
            coeffDict_,
            2.0
        )
    ),
    cP2_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cP2",
            coeffDict_,
            0.3
        )
    ),
    cP3_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cP3",
            coeffDict_,
            0.42
        )
    ),
    cP4_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cP4",
            coeffDict_,
            0.02
        )
    ),
    cP5_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cP5",
            coeffDict_,
            0.0
        )
    ),
    cL1_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cL1",
            coeffDict_,
            0.36
        )
    ),
    cL2_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cL2",
            coeffDict_,
            85.0
        )
    ),
    cN1_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cN1",
            coeffDict_,
            0.6
        )
    ),
    cN2_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cN2",
            coeffDict_,
            0.6
        )
    ),
    cND1_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cND1",
            coeffDict_,
            0.8
        )
    ),
    cND2_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cND2",
            coeffDict_,
            0.6
        )
    ),
    cMu_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cMu",
            coeffDict_,
            0.21
        )
    ),
    betaK_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "betaK",
            coeffDict_,
            0.09
        )
    ),
    cT_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cT",
            coeffDict_,
            0.0033
        )
    ),
    cA_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cA",
            coeffDict_,
            1.0
        )
    ),
    cEhmM_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cEhmM",
            coeffDict_,
            10.0
        )
    ),
	cNF_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cNF",
            coeffDict_,
            1.0
        )
    ),
	cNL_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cNF",
            coeffDict_,
            0.00001
        )
    ),
	pMix_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "pMix",
            coeffDict_,
            0.33
        )
    ),
	cPrK_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cPrK",
            coeffDict_,
            0.67
        )
    ),
	cPrP_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cPrP",
            coeffDict_,
            0.83
        )
    ),
	nutRatMax_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "nutRatMax",
            coeffDict_,
            1.0e5
        )
    ),
    sigmaK_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "sigmaK",
            coeffDict_,
            1.0
        )
    ),
    sigmaEps_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "sigmaEps",
            coeffDict_,
            0.833
        )
    ),
    sigmaPhi_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "sigmaPhi",
            coeffDict_,
            0.33
        )
    ),
    sigmaPsi_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "sigmaPsi",
            coeffDict_,
            1.0
        )
    ),
	prodType_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "prodType",
            coeffDict_,
            1.0
        )
    ),
	nutType_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "nutType",
            coeffDict_,
            3.0
        )
    ),
	cp1Type_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cp1Type",
            coeffDict_,
            1.0
        )
    ),
	psExtraType_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "psExtraType",
            coeffDict_,
            1.0
        )
    ),
   solveK_
   (
       coeffDict_.lookup("solveK")
   ),

   solveEps_
   (
       coeffDict_.lookup("solveEps")
   ),

   solvePsi_
   (
       coeffDict_.lookup("solvePsi")
   ),

   solvePhi_
   (
       coeffDict_.lookup("solvePhi")
   ),

   solveNut_
   (
       coeffDict_.lookup("solveNut")
   ),

   eqncEp1_
   (
       coeffDict_.lookup("eqncEp1")
   ),
   
   eqncEp2_
   (
       coeffDict_.lookup("eqncEp2")
   ),

   eqnEpsHat_
   (
       coeffDict_.lookup("eqnEpsHat")
   ),

   debugWrite_
   (
       coeffDict_.lookup("debugWrite")
   ),
   tslimiter_
   (
       coeffDict_.lookup("tslimiter")
   ),
   lslimiter_
   (
       coeffDict_.lookup("lslimiter")
   ),
   eqncMu_
   (
       coeffDict_.lookup("eqncMu")
   ),
   y_
   (
   mesh_
   ),
    
	k_
    (
        IOobject
        (
            "k",
            runTime_.timeName(),
            U_.db(),
            IOobject::MUST_READ,  
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    
	gradk_
    (
        IOobject
        (
            "gradk",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (fvc::grad(k_))
    ),
    
	epsilon_
    (
        IOobject
        (
            "epsilon",
            runTime_.timeName(),
            U_.db(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    
	nut_
    (
        IOobject
        (
            "nut",
            runTime_.timeName(),
            U_.db(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    
	nutNorm_
    (
        IOobject
        (
            "nutNorm",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (nut_/max(nut_))
    ),
    
	tpphi_
    (
        IOobject
        (
            "tpphi",
            runTime_.timeName(),
            U_.db(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    
	tpphiSqrt_
    (
        IOobject
        (
            "tpphiSqrt",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (sqrt(tpphi_))
    ),
    
	vorticity_
    (
        IOobject
        (
            "vorticity",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        (fvc::curl(U_))
    ),
    
	tppsi_
    (
        IOobject
        (
            "tppsi",
            runTime_.timeName(),
            U_.db(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    
	uGrad_
    (
        IOobject
        (
            "uGrad",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        (fvc::grad(U_))
    ),

	S_
    (
        IOobject
        (
            "uGrad",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        symm(fvc::grad(U_))
    ),
	
	G_
    (
        IOobject
        (
            "G",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        nut_*mag(S_)*mag(S_)
    ),
	
	ndsPhi_
    (
        IOobject
        (
            "ndsPhi",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        tpphi_
    ),
	
	ndsPsi_
    (
        IOobject
        (
            "ndsPsi",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        tpphi_
    ),
	
	wally_
    (
        IOobject
        (
            "wally",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        wallDist(mesh_).y()
    ),
	
	epsHat_
    (
        IOobject
        (
            "epsHat",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        (epsilon_/(k_ + k0_))
    ),
    
	kSqrt_
    (
        IOobject
        (
            "kSqrt",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (sqrt(k_))
    ),
	
	alpha_
    (
        IOobject
        (
            "alpha",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        (1.0/(1.0 + 1.5*tpphi_))
    ),
	
	gamma_
    (
        IOobject
        (
            "gamma",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        (1.0/(1.0 + cG_*nut_/nu()))
    ),
	
	upsilon_
    (
        IOobject
        (
            "upsilon",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        k_
    ),
	
	chi_
    (
        IOobject
        (
            "chi",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        k_
    ),
	
	phiSqrt_
    (
        IOobject
        (
            "phiSqrt",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (sqrt(tpphi_*k_))
    ),
	
	phiProd_
    (
        IOobject
        (
            "phiProd",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (R() && (skew(fvc::grad(U_))))
    ),
    
	gradkSqrt_
    (
        IOobject
        (
            "gradkSqrt",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (fvc::grad(kSqrt_))
    ),
    
	cEp2_
    (
        IOobject
        (
            "cEp2",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (cEp2con_ - 0.16*exp(-0.1*sqr(k_)/(nu()*epsilon_)))
    ),
    
	tpProd_
    (
        IOobject
        (
            "tpProd",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        ((2*nut_*magSqr(symm(fvc::grad(U_)))/k_))
    ),
	
	addedPhiProd_
    (
        IOobject
        (
            "addedPhiProd",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (tpProd_)
    ),
	
	iPPsi_
    (
        IOobject
        (
            "iPPsi",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (tpProd_)
    ),
    
	cP1eqn_
    (
        IOobject
        (
            "cP1eqn",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (2.0*(0.5+0.5*((tpProd_*k_)/epsilon_)))
    ),
    
	dimRat_
    (
        IOobject
        (
            "dimRat",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (psiReal() & psiReal())/(k_*phiReal())
    ),
    
	gradTpphi_
    (
        IOobject
        (
            "gradTpphi",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (fvc::grad(tpphi_))
    ),
    
	gradTppsi_
    (
        IOobject
        (
            "gradTppsi",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (fvc::grad(tppsi_))
    ),
	
	unitGradTpphi_
    (
        IOobject
        (
            "unitGradTpphi",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (gradTpphi_*pow(k_,1.5)/epsilon_)
    ),
    
	tpProdSqr_
    (
        IOobject
        (
            "tpProdSqr",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (sqr(tppsi_ & vorticity_))
    ),
    
	tpProd3d_
    (
        IOobject
        (
            "tpProd3d",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (mag(psiReal() ^ vorticity_))
    ),
	
	f_
    (
        IOobject
        (
            "f",
            runTime_.timeName(),
            U_.db(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

	wdamp_
    (
        IOobject
        (
            "wdamp",
            runTime_.timeName(),
            U_.db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        tpphi_
    ),
	
	Rev_
    (
        IOobject
        (
            "Rev",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (0.667*I*k_ - nut_*twoSymm(fvc::grad(U_)))
    ),

	D_
    (
        IOobject
        (
            "D",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        tpphi_
    ),

	J_
    (
        IOobject
        (
            "J",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        ((0.667*I*k_ - nut_*(fvc::grad(U_)))*epsilon_/pow(k_,5.0/2.0))
    ),
	
	A_
    (
        IOobject
        (
            "A",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::ddt(U_) + fvc::div(phi_, U_) - fvc::div(phi_)*U_
    ),
	
	B_
    (
        IOobject
        (
            "B",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        vorticity_
    ),
	
	C_
    (
        IOobject
        (
            "C",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        vorticity_
    ),
	
	funit_
    (
        IOobject
        (
            "funit",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U_/(mag(U_) + dimensionedScalar("UMinTemp", U_.dimensions(), SMALL))
    ),
	
	nunit_
    (
        IOobject
        (
            "nunit",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U_/(mag(U_) + dimensionedScalar("UMinTemp", U_.dimensions(), SMALL))
    ),
	
	punit_
    (
        IOobject
        (
            "punit",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U_/(mag(U_) + dimensionedScalar("UMinTemp", U_.dimensions(), SMALL))
    ),
	
	cK_
    (
        IOobject
        (
            "cK",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mag(vorticity_)/(mag(U_) + dimensionedScalar("UMinTemp", U_.dimensions(), SMALL))
    ),
	
	newRS_
    (
        IOobject
        (
            "newRS",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        tppsi_*k_
    ),

	ev_
    (
        IOobject
        (
            "ev",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (0.667*I*k_ - nut_*(fvc::grad(U_)))
    ),
	
	Rnfp_
    (
        IOobject
        (
            "Rnfp",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Rev_
    ),
	
	Rcart_
    (
        IOobject
        (
            "Rcart",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Rev_
    ),
	
	Pcart_
    (
        IOobject
        (
            "Pcart",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        symm(((Rev_ & fvc::grad(U_))))
    ),
	
	nfp_
    (
        IOobject
        (
            "nfp",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        uGrad_*k_/epsilon_
    ),

	uGradnfp_
    (
        IOobject
        (
            "uGradnfp",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        uGrad_
    ),
	
	Snfp_
    (
        IOobject
        (
            "Snfp",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        symm(uGradnfp_)
    ),
	
	Pnfp_
    (
        IOobject
        (
            "Pnfp",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        symm(((Rev_ & fvc::grad(U_))))
    ),

	Dnfp_
    (
        IOobject
        (
            "Dnfp",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        symm(((Rev_ & fvc::grad(U_))))
    ),
	
	fastGlm_
    (
        IOobject
        (
            "fastGlm",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        symm(fvc::grad(U_))
    ),
	
	fastGlmVec_
    (
        IOobject
        (
            "fastGlmVec",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        vorticity_
    ),

	tenPCart_
    (
        IOobject
        (
            "tenPCart",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        ((Rev_ & fvc::grad(U_)))
    ),
	
	tenP_
    (
        IOobject
        (
            "tenP",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        ((Rev_ & fvc::grad(U_)))
    )
{

    Info<< "Made it past constructors " << endl;

    // Calculate eddy viscosity
    if(solveNut_ == "true")
    {
		nut_ = cMu_*k_*tpphi_*Ts();	
        nut_ = min(nut_,nutRatMax_*nu());        
        nut_.correctBoundaryConditions();
        bound(nut_,dimensionedScalar("minNut", nut_.dimensions(), SMALL));       
    }
	
	
    //*************************************//	
    // Epsilon-tilda-hat
    //*************************************//
    epsHat_ = epsilon_/(k_ + (cEhmM_*nu()*mag(gradkSqrt_)));
    bound(epsHat_,dimensionedScalar("minEpsHat", epsHat_.dimensions(), SMALL));
	
    Info<< "solveK is: " <<solveK_ <<endl;
    Info<< "solveEps is: " <<solveEps_ <<endl;
    Info<< "solvePhi is: " <<solvePhi_ <<endl;
    Info<< "solvePsi is: " <<solvePsi_ <<endl;

    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Not used but necessary for RAS Model
tmp<volSymmTensorField> v2ftpFinal::R() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "R",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            ((2.0/3.0)*I)*k_ - nut_*twoSymm(fvc::grad(U_))
        )
    );
}

// Not used but necessary for RAS Model
tmp<volSymmTensorField> v2ftpFinal::devReff() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "devRhoReff",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            -nuEff()*dev(twoSymm(fvc::grad(U_)))
        )
    );
}

// Term that is directly added to the momentum equation
tmp<fvVectorMatrix> v2ftpFinal::divDevReff() const
{
    return
    (
       fvc::grad(phiReal())
     + fvc::curl(psiReal())
     + fvc::laplacian(nut_, U_, "laplacian(nuEff,U_)")
     - fvm::laplacian(nuEff(), U_)
    );
}


bool v2ftpFinal::read()
{
    if (RASModel::read())
    {
        cEp1_.readIfPresent(coeffDict());
        cEp2con_.readIfPresent(coeffDict());
        cEp3_.readIfPresent(coeffDict());
		cEp4_.readIfPresent(coeffDict());
        cP1_.readIfPresent(coeffDict());
        cP2_.readIfPresent(coeffDict());
        cP3_.readIfPresent(coeffDict());
		cP4_.readIfPresent(coeffDict());
		cP5_.readIfPresent(coeffDict());
		cL1_.readIfPresent(coeffDict());
		cL2_.readIfPresent(coeffDict());
		cN1_.readIfPresent(coeffDict());
		cN2_.readIfPresent(coeffDict());
        cMu_.readIfPresent(coeffDict());
		cEhmM_.readIfPresent(coeffDict());
		cPrK_.readIfPresent(coeffDict());
		cPrP_.readIfPresent(coeffDict());
		cD1_.readIfPresent(coeffDict());
		cD2_.readIfPresent(coeffDict());
		cG_.readIfPresent(coeffDict());
		cGw_.readIfPresent(coeffDict());
		cT_.readIfPresent(coeffDict());
		cA_.readIfPresent(coeffDict());
		cNF_.readIfPresent(coeffDict());
		sigmaK_.readIfPresent(coeffDict());
        sigmaEps_.readIfPresent(coeffDict());
        sigmaPhi_.readIfPresent(coeffDict());
		sigmaPsi_.readIfPresent(coeffDict());
		prodType_.readIfPresent(coeffDict());
		nutType_.readIfPresent(coeffDict());
		cp1Type_.readIfPresent(coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


void v2ftpFinal::correct()
{
	
	//Info << "Made it to correct" << endl;

    //**********************************************//	
    // Bounding values not already defined by model
    //**********************************************//
	
    const dimensionedScalar eH0("minEpsHat", epsHat_.dimensions(), ROOTVSMALL);
	const dimensionedScalar nut0("minNut", nut_.dimensions(), SMALL);
	const dimensionedScalar tph0("minTpphi", tpphi_.dimensions(), ROOTVSMALL);
	const dimensionedScalar gtph0("minGTpphi", gradTpphi_.dimensions(), ROOTVSMALL);
	const dimensionedScalar f0("fMin", f_.dimensions(), 0.0);
	const dimensionedScalar L0("lMin", dimensionSet(0,1,0,0,0,0,0), ROOTVSMALL);
	const dimensionedScalar pg0("pgMin", dimensionSet(0,1,-2,0,0,0,0), ROOTVSMALL);
	const dimensionedScalar v0("vMin", vorticity_.dimensions(), SMALL);
    const dimensionedScalar U0("UMin", U_.dimensions(), SMALL);

    if (mesh_.changing())
    {
        y_.correct();
        bound(k_, k0_);
        bound(epsilon_, epsilonSmall_);
		bound(tpphi_,tph0);
		bound(nut_,nut0);
    }
	
	
    RASModel::correct();

	
    if (!turbulence_)
    {
        return;
    }
	
	word sMMdebug = runTime_.controlDict().lookup("showMaxMin");


    //*************************************//	
    // Vorticity, Gradient, Invariants
    //*************************************//
    
	vorticity_ = fvc::curl(U_);
	uGrad_ = fvc::grad(U_);
	
	volSymmTensorField symmUgrad("symmUgrad", twoSymm(uGrad_));
	volTensorField skewUgrad("skewUgrad",2.0*skew(uGrad_));
	volTensorField roTen("roTen",2.0*skew(uGrad_.T()));
	
	D_ = ((symmUgrad && symmUgrad) - (skewUgrad && skewUgrad))/((symmUgrad && symmUgrad) + (skewUgrad && skewUgrad));
	
	// volScalarField uII("uII", -0.5*(pow(tr(symmUgrad),2.0) - tr(symmUgrad & symmUgrad)));
	// volScalarField Q("Q", 0.25*((vorticity_ & vorticity_) - 2.0*(symmUgrad && symmUgrad)));
	
	// //iPPsi_ = sqrt(((tppsi_ & tppsi_)+SMALL)*uII);
	
    // Rev_ = ((2.0/3.0)*I)*k_ - nut_*twoSymm(uGrad_);
	// volSymmTensorField aij("aij",-nut_*twoSymm(uGrad_));	
	// volSymmTensorField symP2("symP2",-twoSymm((aij & uGrad_)));
		
	// volScalarField inv2Rev("inv2Rev", pow(tr(Rev_),2.0) - tr(Rev_ & Rev_));		
	// volScalarField invP2("invP2", pow(tr(symP2),2.0) - tr(symP2 & symP2));
	
	// //Info << "Made it past grads" << endl;

	// A_ = fvc::div(phi_, U_) - fvc::div(phi_)*U_;
	// cK_ = mag((U_ ^ A_))/(pow(mag(U_),3.0) + U0*U0*U0);
	
	// funit_ = ( U_/(mag(U_) + U0));
	// volVectorField psicu("psicu", (tppsi_*k_) ^ U_);
	// nunit_ = psicu/(mag(psicu) + k0_*U0);
	// punit_ = ((tppsi_*k_)/(mag(tppsi_*k_) + k0_));
	
	// B_ = funit_ ^ (funit_ & uGrad_);
	// C_ = nunit_ ^ (nunit_ & uGrad_);
	// volScalarField Bs("Bs", funit_ & (funit_ & uGrad_));
	// volScalarField H("H", funit_ & (funit_ & S_));

	
    //*************************************//	
    // Length and Time Scales
    //*************************************//	
	
	const volScalarField L("Length",Ls());
	const volScalarField L2("Lsqr",sqr(L));
	const volScalarField T("Time",Ts());	


	
	//*************************************//	
    // Gradient and Misc Terms
    //*************************************//

	kSqrt_ = sqrt(mag(k_)+k0_);
    bound(kSqrt_,dimensionedScalar("minKsqrt", kSqrt_.dimensions(), sqrt(ROOTVSMALL)));

    gradk_ = fvc::grad(k_);
    gradkSqrt_ = fvc::grad(kSqrt_);
	
	phiSqrt_ = sqrt(tpphi_*k_ + k0_);
	
	const volVectorField gradPhiSqrt_("gradPhiSqrt",fvc::grad(phiSqrt_));
	const volVectorField gradPhi_("gradPhi", fvc::grad(phiReal()));		
	const volScalarField gradgradPhi_("gradgradPhi", fvc::laplacian(DphiEff(),phiReal()));

	gradTpphi_ = fvc::grad(tpphi_);
	tpphiSqrt_ = sqrt(tpphi_ + ROOTVSMALL);
	const volVectorField gradTpphiSqrt("gradTpphiSqrt",fvc::grad(tpphiSqrt_));
	
    gradTppsi_ = fvc::grad(tppsi_);


	
    //*************************************//	
    // K Production
    //*************************************//
 
	S_ = symm(uGrad_);
	const volScalarField S2 = 2*magSqr(dev(S_));
	const volScalarField magS = sqrt(S2);
	
	volScalarField G("RASModel::G", G_);
	volScalarField GdK("GdK", G_/(k_ + k0_));
	const volScalarField Gnut("Gnut", nut_*S2);


	if(prodType_.value() == 1.0){
		Info<< "Using psi-vorticity production term" <<endl;
		tpProd_ = mag(tppsi_ & vorticity_);		
		G_ = tpProd_*k_;
		GdK = tpProd_;		
	}
	
	if(prodType_.value() == 1.1){
		Info<< "Using psi-vorticity prod - no mag" <<endl;
		tpProd_ = tppsi_ & vorticity_;		   
		G_ = tpProd_*k_;
		GdK = tpProd_;		
	}
	
	if(prodType_.value() == 2.0){
		Info<< "Using strain production term" <<endl;
		G_ = Gnut;
		GdK = G_/(k_ + k0_);
		tpProd_ = GdK;
	} 
	
	if(prodType_.value() == 3.0){
		Info<< "Using mixed 3 production term" <<endl;
		tpProd_ = alpha_*mag(tppsi_ & vorticity_) + pMix_*(1.0-alpha_)*cPrK_*magS + (1.0-pMix_)*(1.0-alpha_)*cPrP_*tpphi_*magS;
		G_ = tpProd_*k_;
		GdK = tpProd_;	
    }
	
	if(prodType_.value() == 4.0){
		Info<< "Using mixed 4 production term" <<endl;
		tpProd_ = alpha_*mag(tppsi_ & vorticity_) + 0.88*(1.0-alpha_)*Gnut/(k_ + k0_);
		G_ = tpProd_*k_;
		GdK = tpProd_;
    }
	
	if(prodType_.value() == 5.0){
		Info<< "Using mixed upsilon production term" <<endl;
		tpProd_ = alpha_*mag(tppsi_ & vorticity_) + 0.5*(1.0-alpha_)*(0.33*(upsilon_/(k_ + k0_)) + 1.1*tpphi_)*magS;
		G_ = tpProd_*k_;
		GdK = tpProd_;	
    }
    
	tpProdSqr_ = sqr(tpProd_);
	tpProd3d_ = mag(psiReal() ^ vorticity_);
	
	const volScalarField pOD = G_/epsilon_; 	
	const volScalarField maxpOD = 2.0*epsilon_/epsilon_;

	    
	
	
	//*************************************//	
    // Update Alpha
    //*************************************//
    
	alpha_ = 1.0/(1.0 + 1.5*tpphi_);
	
	//volScalarField IIb("IIb", alpha_*(2.0*alpha_-1.0));	
	volScalarField IIb("IIb", 2.0*(0.5*sqr(2.0*alpha_-1.0) + 0.6*(tppsi_ & tppsi_)));
	bound(IIb, SMALL);
	
	volScalarField phiActual("phiActual",tpphi_*k_);
	volVectorField psiActual("psiActual",tppsi_*k_);
	
	volScalarField nutExact("nutExact", mag(psiActual)/(mag(vorticity_) + (cNL_/T)));
    volScalarField gammaNut("gammaNut", nut_);
	
	//volScalarField gammaNut("gammaNut", (alpha_*(psiActual & psiActual) + 0.57*(1.0-alpha_)*phiActual*phiActual)/(epsHat_*k_));
	
	if(nutType_.value() == 1.0){
		gammaNut = pow(IIb, 0.5)*nutExact + (1.0-pow(IIb, 0.5))*cMu_*phiActual/epsHat_;
	}
	
	if(nutType_.value() == 2.0){
		gammaNut = alpha_*nutExact + (1.0-alpha_)*cMu_*phiActual/epsHat_;
	}
	
	if(nutType_.value() == 3.0){
		gammaNut = alpha_*alpha_*nutExact + (1.0-alpha_*alpha_)*cMu_*phiActual/epsHat_;
	}
	
	if(nutType_.value() == 4.0){
		gammaNut = (alpha_*(psiActual & psiActual) + 0.57*(1.0-alpha_)*phiActual*phiActual)/(epsHat_*k_);
	}
	
	if(nutType_.value() == 5.0){
		gammaNut = alpha_*nutExact + (0.5/cMu_)*(1.0-alpha_)*tpphi_*nut_;
	}
	
	if(nutType_.value() == 6.0){
		gammaNut = IIb*nutExact + (0.5/cMu_)*(1.0-IIb)*tpphi_*nut_;
	}
	
	if(nutType_.value() == 7.0){
		gammaNut = nut_;
	}
	
	volScalarField gammaWall("gammaWall", 3.0*nu()*(gradTpphiSqrt & gradTpphiSqrt)*k_/epsilon_); 

	gamma_ = 1.0/(1.0 + cG_*gammaNut/nu() + cGw_*gammaWall);
	
	
    //*************************************//
    // Epsilon Constant Functions 
    //*************************************//	
    if(eqncEp2_ == "true")
    {
        cEp2_ = cEp2con_ - 0.16*exp(-0.25*sqr(k_)/(nu()*(epsilon_ + epsilonSmall_)));
    }
    else
    {
        cEp2_ =  cEp2con_;
    }
	
	volScalarField cEp1eqn("cEp1eqn",cEp1_*(epsilon_/epsilon_));
	
	if(eqncEp1_ == "true")
	{
		Info << "Using Eqn Cep1: " << cEp1_ << " Cep4: " << cEp4_ << endl;
		cEp1eqn = cEp1_ + cEp4_*(2.0*alpha_-1.0);		
	}
	
	

	//*************************************//	
    // Update Epsilon-hat
    //*************************************//
    
    epsHat_ = epsilon_/(k_ + (cEhmM_*nu()*mag(gradkSqrt_)));
    bound(epsHat_,eH0);	
		
	epsilon_.boundaryField().updateCoeffs();
	
	
	
    //*************************************//
    //Dissipation equation
    //*************************************//


    tmp<fvScalarMatrix> epsEqn  
    (
        fvm::ddt(epsilon_)
      + fvm::div(phi_, epsilon_)
      + fvm::SuSp(-fvc::div(phi_), epsilon_)
      - fvm::laplacian(DepsilonEff(), epsilon_)
     ==
       cEp1eqn*G_/T
     - fvm::Sp(cEp2_/T,epsilon_)
     + cEp3_*tpProd3d_/T
    );

    if(solveEps_ == "true")
    {
    epsEqn().relax();
    solve(epsEqn);
    bound(epsilon_,epsilonSmall_);
    }
	
	
	
	
    //*************************************//
    // Turbulent kinetic energy equation
    //*************************************//
    
    tmp<fvScalarMatrix> kEqn 
    (

        fvm::ddt(k_)
      + fvm::div(phi_, k_)
      + fvm::SuSp(-fvc::div(phi_), k_)
      - fvm::laplacian(DkEff(), k_)
     ==
        G_
      - fvm::Sp(epsilon_/(k_+k0_),k_)
    );


    if(solveK_ == "true")
    {
    kEqn().relax();
    solve(kEqn);
    bound(k_,k0_);
    }
	

	
	//*************************************//   
    // Axial Reynolds stress estimation
    //*************************************//
	
	const volScalarField rII("rII", sqrt(phiActual*phiActual + (psiActual & psiActual) + k0_*k0_));
	
	chi_ = 2.0*alpha_*rII;
	upsilon_ = (4.0/3.0)*(1.75*alpha_-0.375)*k_ ;

	const volScalarField uudk("uudk",upsilon_/(k_ + k0_));
	const volScalarField wwdk("wwdk",chi_/k_);
	
	const volScalarField bup("bup", upsilon_/k_ - (2.0/3.0));
    const volScalarField bph("bph", tpphi_ - (2.0/3.0));
	const volScalarField bch("bch", chi_/k_ - (2.0/3.0));
	
	volScalarField Det("Det", (27.0/8.0)*(uudk*tpphi_*wwdk - wwdk*(tppsi_ & tppsi_)));
	bound(Det, SMALL);
	

	
    //*************************************// 
    // f equation - with elliptic
    //*************************************//
	
	cP1eqn_ = cP1_ - gamma_;
	
	
	if(cp1Type_.value() == 1.0){
		cP1eqn_ = cP1_*tpphi_/tpphi_;
	}

	if(cp1Type_.value() == 2.0){
		cP1eqn_ = cP1_ - gamma_;
	}

	if(cp1Type_.value() == 3.0){
		cP1eqn_ = cP1_ - 0.2 + 0.2*((psiActual & psiActual)/(cMu_*phiActual*k_));
	}	
	
	if(cp1Type_.value() == 4.0){
		cP1eqn_ = cP1_ - sqr(gamma_);
	}
	
	if(cp1Type_.value() == 5.0){
		cP1eqn_ = cP1_*(1.0-0.5*gamma_*sqrt(IIb));    
	}
	
	if(cp1Type_.value() == 6.0){
		cP1eqn_ = cP1_*(1.0 - 0.4*gamma_ + 0.15*((psiActual & psiActual)/(cMu_*phiActual*k_)));
	}		
	
	if(cp1Type_.value() == 7.0){
		cP1eqn_ = cP1_*(1.4 - sqrt(IIb));    
	}
	
    //*************************************//
    // Fix for separated flows 
    //*************************************// 	

	volScalarField cTexp("cTexp", cT_*(1.0-gamma_)*sqrt((((nu()/100.0)+nut_)/nu())));
	
    volScalarField transPhi("transPhi", cTexp*cA_*((2.0/3.0) - tpphi_)*tpProd_);	
	volVectorField transPsi("transPsi", cTexp*((1.0 - alpha_)*vorticity_ - cA_*tppsi_*tpProd_));
	
	
	
	
	
    //*************************************//
    // Pressure strain terms phi 
    //*************************************// 		
	volScalarField slowPS 
    (
        "v2ftpFinal::slowPS",
        -(cP1eqn_-(1.0-gamma_))*bph/T 	 
	);
	
	volScalarField slowPSnonlin
    (
        "v2ftpFinal::slowPSnonlin",
        cD2_*(sqr(bph) + (tppsi_ & tppsi_) - (2.0/3.0)*IIb)/T	 
	);
	
	volScalarField fastPS
    (
        "v2ftpFinal::fastPS",
		0.667*cP2_*GdK
	);
	
	   
	volScalarField fwall 
    (
        "v2ftpFinal::fwall",
		IIb*epsHat_*tpphi_
	); 
		

    tmp<fvScalarMatrix> fEqn
    (
      - fvm::laplacian(f_)
     ==
      - fvm::Sp(1.0/L2, f_)
	  + fwall/L2
      + slowPS/L2
	  + slowPSnonlin/L2	  
	  + fastPS/L2
	  + transPhi/L2
    );
 
    fEqn().relax();
    solve(fEqn);
	
    
	


    //*************************************//
    // Phi/K equation - with elliptic
    //*************************************// 

    tmp<fvScalarMatrix> tpphiEqn
    (
        fvm::ddt(tpphi_)
      + fvm::div(phi_, tpphi_)
	  + fvm::SuSp(-fvc::div(phi_), tpphi_)
      - fvm::laplacian(DphiEff(), tpphi_)
      ==
	    //Production
		//phiProd_/(k_+k0_)
		
	    //Source - pressure strain
        min(f_,(slowPS + fastPS + fwall + transPhi))
		   
        //BC wall correct	
	  - fvm::Sp(fwall/(tpphi_+tph0),tpphi_)
	  
	    //From k eqn phi/k derivation
      - fvm::Sp(GdK, tpphi_)

	    // Dissipation included in pressure strain term
	    //+ (1.0-gamma_)*bph/T
    ); 
	

    tpphiEqn().relax();
    solve(tpphiEqn);
	bound(tpphi_,tph0);

	
	
	
    //*************************************//   
    // Psi Equation
    //*************************************//
	volVectorField vecProd("vecProd", phiActual*vorticity_);
	volVectorField curvProd("curvProd", upsilon_*B_);
	volVectorField addedPsiProd("addedPsiProd", tppsi_ & roTen);
	
    
	volVectorField psExtra("psExtra", cP3_*(1.0-sqrt(IIb))*vorticity_);
	
	if(psExtraType_.value() == 1.0){
		psExtra = cP3_*(1.0-sqrt(IIb))*vorticity_;
	}
	
	if(psExtraType_.value() == 2.0){
		psExtra = cP3_*(1.0/6.0)*(pow(Det,(2.0/3.0)))*vorticity_;
	}
	
	if(psExtraType_.value() == 3.0){
		psExtra = cP3_*(1.0-alpha_)*vorticity_;
	}	

	
	//volVectorField psiDisWall("psiDisWall", cD1_*sqr(gamma_)*psExtra);
    volVectorField psiDisWall("psiDisWall", cD1_*gamma_*(1.0-sqrt(IIb))*vorticity_);
	
    tmp<fvVectorMatrix> tppsiEqn
    (
        fvm::ddt(tppsi_)
      + fvm::div(phi_, tppsi_)
      + fvm::Sp(-fvc::div(phi_), tppsi_)
      - fvm::laplacian(DpsiEff(), tppsi_) 

      == 

	  // Production
	    vecProd/(k_+k0_)
	  //+ addedPsiProd //3d Psi production

	  // Slow Pressure Strain
      - fvm::Sp(cP1eqn_/T,tppsi_)
	  + cD2_*(bph + (uudk-(2.0/3.0)))*tppsi_/T
	      
	  // Fast Pressure Strain      
	  - cP2_*vecProd/(k_+k0_)
	  - (cP3_-cD1_*gamma_)*(1.0-sqrt(IIb))*vorticity_  
	  
	  // Extra term for fixing hump recirc zone
	  + cP4_*((1.0-gamma_)/sqrt(1.12-alpha_))*sqrt((epsilon_ + epsilonSmall_)/nu())*tppsi_
	  
	  // Dissipation 
	  + (1.0-gamma_)*tppsi_/T
	  // Near wall included in Fast Pressure Strain 
	  
	  // From K equation
	  - fvm::Sp(tpProd_,tppsi_)

	  // Transition Term
      + transPsi 
    );

    if(solvePsi_ == "true")
    {
		tppsiEqn().relax();
		solve(tppsiEqn);    
    }  

	
	psiActual = tppsi_*k_;
	phiActual = tpphi_*k_;

    
	//*************************************//
    // Calculate eddy viscosity
    //*************************************//

	//nut_ = cMu_*(0.6*phiActual + 2.2*(psiActual & tppsi_))*k_/epsilon_;	

	nut_ = (cN1_ + (1.0 - cN1_)*(psiActual & psiActual)/(cMu_*phiActual*k_))*cMu_*phiActual/epsHat_;  
	
	//nut_ = cMu_*(cN1_ + cN2_*sqrt(IIb))*tpphi_*k_*T;
    
	nut_ = min(nut_,nutRatMax_*nu()); 
	nut_.correctBoundaryConditions();
    bound(nut_,nut0); 
	
	ndsPhi_ = cND1_ + (1.0 - cND1_)*(psiActual & psiActual)/(cMu_*phiActual*k_);
	ndsPsi_ = cND2_ + (1.0 - cND2_)*(psiActual & psiActual)/(cMu_*phiActual*k_);
	
	
    //*************************************//   
    // Output some min/max debug values
    //*************************************//
	
	if(sMMdebug == "true")
	{
    
	volScalarField uTauSquared((nu() + nut_)*vorticity_.component(2));
	volVectorField addedPsiG(addedPsiProd*k_);
	volVectorField phiOmega(tpphi_*k_*vorticity_);
	

	Info << " * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *" << endl;
	Info << "Max psiDisWall: " << gMax(psiDisWall) << " Min psiDisWall: " << gMin(psiDisWall) << endl;
	Info << "Max phiOmega: " << gMax(phiOmega) << " Min phiOmega: " << gMin(phiOmega) << endl;
	Info << "Max cEps1: " << gMax(cEp1eqn) << " Min cEps1: " << gMin(cEp1eqn) << endl;
	Info<< "Max f: " << gMax(f_) << " Min f: " << gMin(f_) << " Max G: " << gMax(G) << " Min G: " << gMin(G) << " Max Gnut: " << gMax(Gnut) << endl;
    Info<< "Max nut: " << gMax(nut_) << " Max K: " << gMax(k_) << " Max Epsilon: " << gMax(epsilon_) << " Max Phi: " << gMax(phiActual) <<endl;
    Info<< "Max Psi: " << gMax(psiActual) << " Min Psi: " << gMin(psiActual) << endl;
	Info<< "Max vorticity: " << gMax(vorticity_) << " Min vorticity: " << gMin(vorticity_) << endl;
	Info << " * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *" << endl;  
	}
	  
	

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
