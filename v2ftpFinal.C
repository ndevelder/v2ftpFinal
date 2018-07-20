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

// Hello

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
            1.45
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
    cD1_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cD1",
       	    coeffDict_,
            0.5
        )
    ),
    cD2_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cD2",
       	    coeffDict_,
            0.12
        )
    ),
    cD3_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cD3",
       	    coeffDict_,
            0.5
        )
    ),
    cD4_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cD4",
       	    coeffDict_,
            1.12
        )
    ),
    cVv1_
    (
     	dimensionedScalar::lookupOrAddToDict
        (
            "cVv1",
            coeffDict_,
            0.0
        )
    ),
    cTv1_
    (
     	dimensionedScalar::lookupOrAddToDict
        (
            "cTv1",
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
            0.6
        )
    ),
    cP3_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cP3",
            coeffDict_,
            0.12
        )
    ),
    cP4_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cP4",
            coeffDict_,
            0.85714
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
    eC1_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "eC1",
            coeffDict_,
            1.4
        )
    ),
    eC2_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "eC2",
            coeffDict_,
            0.43
        )
    ),
    eC3_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "eC3",
            coeffDict_,
            0.0
        )
    ),
    eC4_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "eC4",
            coeffDict_,
            0.0
        )
    ),
    eC5_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "eC5",
            coeffDict_,
            0.0
        )
    ),
    eC6_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "eC6",
            coeffDict_,
            0.3
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
	cPa_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cPa",
            coeffDict_,
            1.6
        )
    ),
	cPp_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cPp",
            coeffDict_,
            0.0
        )
    ),
    cPr_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cPr",
            coeffDict_,
            0.33
        )
    ),
	cPw_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cPw",
            coeffDict_,
            25.0
        )
    ),
	cPD_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cPD",
            coeffDict_,
            0.0
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
    cEhmP_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cEhmP",
            coeffDict_,
            0.67
        )
    ),
    cEhmPK_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cEhmPK",
            coeffDict_,
            0.09
        )
    ),
    cEhmPK2_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cEhmPK",
            coeffDict_,
            0.67
        )
    ),
    cEhR_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cEhR",
            coeffDict_,
            1.0
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
	rS_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "rS",
            coeffDict_,
            1e-10
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
	rPr_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "rPr",
            coeffDict_,
            0.5
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
    sigmaKInit_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "sigmaKInit",
            coeffDict_,
            1.0
        )
    ),
    sigmaEpsInit_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "sigmaEpsInit",
            coeffDict_,
            0.833
        )
    ),
    sigmaEpsVisc_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "sigmaEpsVisc",
            coeffDict_,
            1.0
        )
    ),
    sigmaPhiInit_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "sigmaPhiInit",
            coeffDict_,
            0.33
        )
    ),
    sigmaPsiInit_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "sigmaPsiInit",
            coeffDict_,
            1.0
        )
    ),
    nutScale_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "nutScale",
            coeffDict_,
            1.0
        )
    ),
    nutBlend_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "nutBlend",
            coeffDict_,
            1.0
        )
    ),
    psiNuFrac_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "psiNuFrac",
            coeffDict_,
            1.0
        )
    ),
    ellipticSwitch_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "ellipticSwitch",
            coeffDict_,
            0.0
        )
    ),
	coordType_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "coordType",
            coeffDict_,
            2.0
        )
    ),
	pDir_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "pDir",
            coeffDict_,
            1.0
        )
    ),
	fastPsType_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "fastPsType",
            coeffDict_,
            1.0
        )
    ),
	slowPsType_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "slowPsType",
            coeffDict_,
            1.0
        )
    ),
	fWallType_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "fWallType",
            coeffDict_,
            1.0
        )
    ),	
	apType_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "apType",
            coeffDict_,
            1.0
        )
    ),	
	pr_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "pr",
            coeffDict_,
            0.7
        )
    ),
	prt_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "prt",
            coeffDict_,
            0.9
        )
    ),
	
   nutType_
   (
       coeffDict_.lookup("nutType")
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

   solveT_
   (
       coeffDict_.lookup("solveT")
   ),
   
   eqnSigmaK_
   (
       coeffDict_.lookup("eqnSigmaK")
   ),

   eqnSigmaEps_
   (
       coeffDict_.lookup("eqnSigmaEps")
   ),

   eqnSigmaPhi_
   (
       coeffDict_.lookup("eqnSigmaPhi")
   ),

   eqnSigmaPsi_
   (
       coeffDict_.lookup("eqnSigmaPsi")
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
   
   timeScaleEps_
   (
       coeffDict_.lookup("timeScaleEps")
   ),
   prodType_
   (
       coeffDict_.lookup("prodType")
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
   phiType_
   (
       coeffDict_.lookup("phiType")
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
	
	T_
    (
        IOobject
        (
            "T",
            runTime_.timeName(),
            U_.db(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
	
	S_
    (
        IOobject
        (
            "uGrad",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        symm(fvc::grad(U_))
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
            IOobject::NO_WRITE
        ),
        (1.0/(1.0 + 1.5*tpphi_))
    ),
	
	upsilon_
    (
        IOobject
        (
            "upsilon",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
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
            IOobject::AUTO_WRITE
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
            IOobject::AUTO_WRITE
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
    
	sigmaK_
    (
        IOobject
        (
            "sigmaK",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        (sigmaKInit_)
    ),
    
	sigmaEps_
    (
        IOobject
        (
            "sigmaEps",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        (sigmaEpsInit_)
    ),
    
	sigmaPhi_
    (
        IOobject
        (
            "sigmaPhi",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        (sigmaPhiInit_)
    ),
    
	sigmaPsi_
    (
        IOobject
        (
            "sigmaPsi",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        (sigmaPsiInit_)
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
            IOobject::AUTO_WRITE
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
            IOobject::AUTO_WRITE
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
            IOobject::AUTO_WRITE
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
            IOobject::AUTO_WRITE
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
            IOobject::AUTO_WRITE
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
            IOobject::AUTO_WRITE
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
            IOobject::AUTO_WRITE
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
            IOobject::AUTO_WRITE
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
            IOobject::AUTO_WRITE
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
		if(timeScaleEps_ == "epsilon" || timeScaleEps_ != "epsHat")
		{
			if(nutType_ == "strain"){
				nut_ = 0.09*k_*Ts();
			}else{
				nut_ = cMu_*k_*tpphi_*Ts();	
			}
        }
        
        if(timeScaleEps_ == "epsHat")
		{
            nut_ = cMu_*k_*tpphi_/epsHat_;
        }

        nut_ = min(nut_,nutRatMax_*nu());        
        nut_.correctBoundaryConditions();
        bound(nut_,dimensionedScalar("minNut", nut_.dimensions(), SMALL));       
    }
	
	
    //*************************************//	
    // Epsilon-tilda-hat
    //*************************************//
    
	if(eqnEpsHat_ == "mod")
	{
        epsHat_ = epsilon_/(k_ + (cEhmM_*nu()*mag(gradkSqrt_)));
        bound(epsHat_,dimensionedScalar("minEpsHat", epsHat_.dimensions(), SMALL));
	}
	else if(eqnEpsHat_ == "pk")
	{
        epsHat_ = epsilon_/(k_ + (cEhmPK_*nu()*mag(gradk_)/(tpphi_*kSqrt_ + sqrt(k0_))));
        bound(epsHat_,dimensionedScalar("minEpsHat", epsHat_.dimensions(), SMALL));
	}
	else if(eqnEpsHat_ == "pk2")
	{
        epsHat_ = epsilon_/(k_ + (cEhmPK2_*nu()*mag(gradk_)/(phiSqrt_ + sqrt(k0_)))); 
        bound(epsHat_,dimensionedScalar("minEpsHat", epsHat_.dimensions(), SMALL));
	}
	else if(eqnEpsHat_ == "phi")
	{		
		volVectorField gradphiSqrt("gradphiSqrt", fvc::grad(sqrt(tpphi_*k_))) ;
        epsHat_ = ((epsilon_)/(1.0 + cEhmP_*nu()*mag(gradphiSqrt)/(tpphi_*k_)))/(k_+k0_);
        bound(epsHat_,dimensionedScalar("minEpsHat", epsHat_.dimensions(), SMALL));
	}
	else if(eqnEpsHat_ == "dif")
	{
        epsHat_ = (epsilon_ - 2.0*nu()*sqr(mag(gradkSqrt_)))/k_;
        bound(epsHat_,dimensionedScalar("minEpsHat", epsHat_.dimensions(), SMALL));
	}
	else if(eqnEpsHat_ == "rough")
	{
        epsHat_ = epsilon_/(k_+k0_);
        bound(epsHat_,dimensionedScalar("minEpsHat", epsHat_.dimensions(), SMALL));
	}
	else
	{
        Info<< "No EpsHat Model Chosen" <<endl;
	    epsHat_ = (epsilon_)/(k_ + cEhmM_*nu()*mag(fvc::grad(kSqrt_)));
	    bound(epsHat_,dimensionedScalar("minEpsHat", epsHat_.dimensions(), SMALL));
	}

	
	
	
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
        cP1_.readIfPresent(coeffDict());
        cP2_.readIfPresent(coeffDict());
        cP3_.readIfPresent(coeffDict());
		cP4_.readIfPresent(coeffDict());
		cL1_.readIfPresent(coeffDict());
		cL2_.readIfPresent(coeffDict());
        cMu_.readIfPresent(coeffDict());
		eC1_.readIfPresent(coeffDict());
		eC2_.readIfPresent(coeffDict());
		eC3_.readIfPresent(coeffDict());
		eC4_.readIfPresent(coeffDict());
		eC5_.readIfPresent(coeffDict());
		cPa_.readIfPresent(coeffDict());
		cPp_.readIfPresent(coeffDict());
		apType_.readIfPresent(coeffDict());
		cEhmP_.readIfPresent(coeffDict());
		cEhmM_.readIfPresent(coeffDict());
		cEhmPK_.readIfPresent(coeffDict());
		cEhmPK2_.readIfPresent(coeffDict());
		cEhR_.readIfPresent(coeffDict());
		cPr_.readIfPresent(coeffDict());
		cPrK_.readIfPresent(coeffDict());
		cPrP_.readIfPresent(coeffDict());
		cD1_.readIfPresent(coeffDict());
		cD2_.readIfPresent(coeffDict());
		cD3_.readIfPresent(coeffDict());
		cD4_.readIfPresent(coeffDict());
        cVv1_.readIfPresent(coeffDict());
        cTv1_.readIfPresent(coeffDict());
		cT_.readIfPresent(coeffDict());
		cNF_.readIfPresent(coeffDict());
		sigmaKInit_.readIfPresent(coeffDict());
        sigmaEpsInit_.readIfPresent(coeffDict());
        sigmaEpsVisc_.readIfPresent(coeffDict());
        sigmaPhiInit_.readIfPresent(coeffDict());
		sigmaPsiInit_.readIfPresent(coeffDict());
		nutBlend_.readIfPresent(coeffDict());
		nutScale_.readIfPresent(coeffDict());
		apType_.readIfPresent(coeffDict());

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
	const dimensionedScalar nut0("minNut", nut_.dimensions(), ROOTVSMALL);
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
	
	
    //*************************************//	
    // Timestep - for use in elliptic switch
    //*************************************//
	
    dimensionedScalar cTime = U_.mesh().time().value();
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
	volScalarField uII("uII", -0.5*(pow(tr(symmUgrad),2.0) - tr(symmUgrad & symmUgrad)));
	volScalarField Q("Q", 0.25*((vorticity_ & vorticity_) - 2.0*(symmUgrad && symmUgrad)));
	
	//iPPsi_ = sqrt(((tppsi_ & tppsi_)+SMALL)*uII);
	
    Rev_ = ((2.0/3.0)*I)*k_ - nut_*twoSymm(uGrad_);
	volSymmTensorField aij("aij",-nut_*twoSymm(uGrad_));	
	volSymmTensorField symP2("symP2",-twoSymm((aij & uGrad_)));
		
	volScalarField inv2Rev("inv2Rev", pow(tr(Rev_),2.0) - tr(Rev_ & Rev_));		
	volScalarField invP2("invP2", pow(tr(symP2),2.0) - tr(symP2 & symP2));
	
	//Info << "Made it past grads" << endl;

	A_ = fvc::div(phi_, U_) - fvc::div(phi_)*U_;
	cK_ = mag((U_ ^ A_))/(pow(mag(U_),3.0) + U0*U0*U0);
	
	funit_ = ( U_/(mag(U_) + U0));
	volVectorField psicu("psicu", (tppsi_*k_) ^ U_);
	nunit_ = psicu/(mag(psicu) + k0_*U0);
	punit_ = ((tppsi_*k_)/(mag(tppsi_*k_) + k0_));
	
	B_ = funit_ ^ (funit_ & uGrad_);
	C_ = nunit_ ^ (nunit_ & uGrad_);
	volScalarField Bs("Bs", funit_ & (funit_ & uGrad_));
	volScalarField H("H", funit_ & (funit_ & S_));

	
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
	volScalarField G("RASModel::G", nut_*S2); 
	volScalarField GdK("GdK", G/(k_ + k0_));
	const volScalarField Gnut("Gnut", nut_*S2);
	

	if(prodType_ == "strain"){
		Info<< "Using strain production term" <<endl;
		tpProd_ = GdK;
	} else if(prodType_ == "mixed3"){
		Info<< "Using mixed 3 production term" <<endl;
		tpProd_ = alpha_*mag(tppsi_ & vorticity_) + pMix_*(1.0-alpha_)*cPrK_*alpha_*magS + (1.0 - pMix_)*(1.0 - alpha_)*cPrP_*tpphi_*magS;
		G = tpProd_*k_;
		GdK = tpProd_;	
    } else if(prodType_ == "curvp"){
		Info<< "Using curvp production term" <<endl;
		tpProd_ = mag((tppsi_ & vorticity_) + (tppsi_ & B_));
		G = tpProd_*k_;
		GdK = tpProd_;	
    } else if(prodType_ == "psic"){
		Info<< "Using psic production term" <<endl;
		tpProd_ = mag(tppsi_ & C_);
		G = tpProd_*k_;
		GdK = tpProd_;	
    } else if(prodType_ == "twocomp"){
		Info<< "Using twocomp production term" <<endl;
		tpProd_ = mag((tppsi_ & vorticity_) + upsilon_*mag((funit_ & uGrad_) & funit_)/(k_ + k0_));
		G = tpProd_*k_;
		GdK = tpProd_;	
    } else if(prodType_ == "mixed4"){
		Info<< "Using mixed 4 production term" <<endl;
		tpProd_ = alpha_*mag(tppsi_ & vorticity_) + 0.88*(1.0-alpha_)*GdK;
		G = tpProd_*k_;
		GdK = tpProd_;
    } else if(prodType_ == "tenp"){
		Info<< "Using tensor P production term" <<endl;
		G = 0.5*mag(tr(Pnfp_));
		tpProd_ = G/(k_+k0_);		
		GdK = tpProd_;		
	} else{
		Info<< "Using psi-vorticity production term" <<endl;
		tpProd_ = tppsi_ & vorticity_;
		G = tpProd_*k_;
		GdK = tpProd_;		
	}
    
	tpProdSqr_ = sqr(tpProd_);
	tpProd3d_ = mag(psiReal() ^ vorticity_);
	
	const volScalarField pOD = G/epsilon_; 
	
	const volScalarField maxpOD = 2.0*epsilon_/epsilon_;
	
	

	
	
    //*************************************//
    // Sigma Functions 
    //*************************************//
    
    if(eqnSigmaK_ == "true")
    {
	    sigmaK_ = 0.67+ 0.4*max(pOD, maxpOD);
	}

    if(eqnSigmaEps_ == "true")
    {
	    sigmaEps_ = 0.33 + 0.4*max(pOD, maxpOD); 
    }

    if(eqnSigmaPhi_ == "true")
    {
	    sigmaPhi_ = 0.5 + 0.4*max(pOD, maxpOD);
	}

    if(eqnSigmaPsi_ == "true")
    {
	    sigmaEps_ = 0.67 + 0.4*max(pOD, maxpOD);
    }

	
	
	
	//*************************************//	
    // Update Alpha
    //*************************************//
    
	alpha_ = 1.0/(1.0 + 1.5*tpphi_);
	volScalarField a1("a1", 1.0 - alpha_);
	volScalarField a2("a2", 2.0*alpha_ - 1.0);
	volScalarField aa("aa", alpha_*alpha_);
	
	
	
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
		cEp1eqn = min(2.0*(tpphi_/tpphi_),(cEp1_)*(1.0+0.09*(alpha_-0.42)));
	}
	
	
	

	//*************************************//	
    // Update Epsilon-hat
    //*************************************//
    
    epsHat_ = epsilon_/(k_ + (cEhmM_*nu()*mag(gradkSqrt_)));
    bound(epsHat_,eH0);	
	
  
	if(eqnEpsHat_ == "pk")
	{
        epsHat_ = epsilon_/(k_ + (2.0*cEhmPK_*nu()*mag(gradkSqrt_)/(tpphi_ + SMALL)));
        bound(epsHat_,eH0);
	}
	else
	{
        epsHat_ = epsilon_/(k_ + (cEhmM_*nu()*mag(gradkSqrt_)));
        bound(epsHat_,eH0);
	}

		
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
       cEp1eqn*G*epsHat_
     - fvm::Sp(cEp2_*epsHat_,epsilon_)
     + cEp3_*tpProd3d_*epsHat_
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
        G
      - fvm::Sp(epsilon_/(k_+k0_),k_)
    );


    if(solveK_ == "true")
    {
    kEqn().relax();
    solve(kEqn);
    bound(k_,k0_);
    }
	

	//*************************************//   
    // Psi Specific Constants
    //*************************************//
	
	const volScalarField phiNow("phiNow", tpphi_*k_);
	const volVectorField psiNow("psiNow", tppsi_*k_);
	chi_ = 2.0*alpha_*sqrt(phiNow*phiNow + (psiNow & psiNow));
	upsilon_ = 2.0*k_ - tpphi_*k_ - chi_ ;

	const volScalarField bup("bup", upsilon_/k_ - (2.0/3.0));
    const volScalarField bph("bph", tpphi_ - (2.0/3.0));
	const volScalarField bch("bch", chi_/k_ - (2.0/3.0));

	volVectorField ucv("ucv", U_ ^ vorticity_);
	volVectorField ucph("ucph", U_ ^ gradPhi_);
	
	volVectorField psiunit("psiunit", psiNow/(mag(psiNow) + k0_));
	volVectorField vunit("vunit", (vorticity_/(mag(vorticity_)+v0))); 

	volScalarField fini("fini", funit_.component(0)*nunit_.component(0));
    volScalarField finj("finj", funit_.component(0)*nunit_.component(1));
	volScalarField fjnj("fjnj", funit_.component(1)*nunit_.component(1));
	volScalarField fjni("fjni", funit_.component(1)*nunit_.component(0));	
	volScalarField nini("nini", nunit_.component(0)*nunit_.component(0));
    volScalarField fifi("fifi", funit_.component(0)*funit_.component(0));
	volScalarField njnj("njnj", nunit_.component(1)*nunit_.component(1));
	volScalarField fjfj("fjfj", funit_.component(1)*funit_.component(1));	
	volScalarField fifj("fifj", funit_.component(0)*funit_.component(1));
	volScalarField ninj("ninj", nunit_.component(0)*nunit_.component(1));
	
	volVectorField vecAdd("vecAdd", punit_/T);
	vecAdd.replace(0,0.0*uGrad_.component(2));
	vecAdd.replace(1,0.0*uGrad_.component(2));
	vecAdd.replace(2,uGrad_.component(1));
	
	volVectorField vecPrim("vecPrim", punit_/T);
	vecPrim.replace(0,0.0*uGrad_.component(2));
	vecPrim.replace(1,0.0*uGrad_.component(2));
	vecPrim.replace(2,uGrad_.component(3));
		
	//volVectorField vecProd("vecProd", (phiNow*vorticity_ + cPp_*(upsilon_ - phiNow)*C_));
	//volVectorField vecProd("vecProd", (phiNow*C_));
	//volVectorField vecProd("vecProd", (phiNow*vorticity_ + (upsilon_ - phiNow)*vecAdd));
	//volVectorField vecProd("vecProd", -1.0*(njnj*phiNow + fjfj*upsilon_)*vecPrim + (nini*phiNow + fifi*upsilon_)*vecAdd);
	//volVectorField vecProd("vecProd", phiNow*C_ - upsilon_*B_);


	
    //*************************************// 
    // f equation - with elliptic
    //*************************************//
	
	phiProd_ = -2.0*((psiNow & B_) + phiNow*Bs);
	const volScalarField c1Comb("c1Comb", cP1_ + (3.0*eC4_)*cP1_*max(pOD, maxpOD)); 
		
	volScalarField slowPS 
    (
        "v2ftpFinal::slowPS",
        -c1Comb*epsHat_*bph + (eC5_/12.0)*epsHat_*((tppsi_ & tppsi_) - bup*bup + 2.0*bph*bph)	
	);
	
	volScalarField fastPS
    (
        "v2ftpFinal::fastPS",
		cP2_*(2.0/3.0)*GdK
	);
	
	   
	volScalarField fwall
    (
        "v2ftpFinal::fwall",
		(2.0*alpha_-1.0)*epsHat_*tpphi_
	);
		

    tmp<fvScalarMatrix> fEqn
    (
      - fvm::laplacian(f_)
     ==
      - fvm::Sp(1.0/L2, f_)
	  + fwall/L2
      + slowPS/L2 
	  + fastPS/L2        
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
        min(f_,(slowPS + fastPS + fwall))
		
        //BC wall correct	  
	  - fvm::Sp(fwall/(tpphi_+tph0),tpphi_)
	  
	    //From k eqn phi/k derivation
      - fvm::Sp(GdK, tpphi_)
	  + epsHat_*tpphi_

	    // Dissipation 
	  - fvm::Sp((cD1_ + 2.0)*alpha_*epsHat_, tpphi_)
	  - cD2_*(1.0-alpha_)*epsHat_
	  //+ (cVv1_*nu())*(gradk_ & gradTpphi_)/(k_+k0_) Cross diffusion left out
    );
	

    tpphiEqn().relax();
    solve(tpphiEqn);
	bound(tpphi_,tph0);

	
	
    //*************************************//   
    // Psi Equation
    //*************************************//
	volVectorField vecProd("vecProd", phiNow*vorticity_);
	volVectorField addedPsiProd("addedPsiProd", tppsi_ & roTen);
	wdamp_ = nutFrac();
	
    
    tmp<fvVectorMatrix> tppsiEqn
    (
        fvm::ddt(tppsi_)
      + fvm::div(phi_, tppsi_)
      + fvm::Sp(-fvc::div(phi_), tppsi_)
      - fvm::laplacian(DpsiEff(), tppsi_)

      ==

	  // Production
	    vecProd/(k_+k0_) 
	  + addedPsiProd //3d Psi production

	  // Slow Pressure Strain
      - fvm::Sp(c1Comb*epsHat_,tppsi_)
	  
	  // Fast Pressure Strain      
	  - cP2_*vecProd/(k_+k0_)
	  
	  // Dissipation
	  - fvm::Sp(cD3_*alpha_*epsHat_,tppsi_)
	  
	  // From K equation
	  + epsHat_*tppsi_
	  - fvm::Sp(tpProd_,tppsi_) 
	  
	  // Gradients
      //+ (cTv1_*nut_)*(gradk_ & gradTppsi_)/(k_+k0_)

	  // Transition Term
      + cT_*sqrt((((nu()/100.0)+nut_)/nu()))*vorticity_
    );

    if(solvePsi_ == "true")
    {
		tppsiEqn().relax();
		solve(tppsiEqn);   
    }  


	
	//*************************************//
    // Calculate eddy viscosity
    //*************************************//
    
    if(solveNut_ == "true")
    {       
		nut_ = cMu_*k_*tpphi_*T;	  
        nut_ = min(nut_,nutRatMax_*nu()); 
		nut_.correctBoundaryConditions();
        bound(nut_,nut0);
    }	
	
	//*************************************//
    // Solve temperature as passive scalar
    //*************************************//
	
	volScalarField kappaEff
	(
		"kappaEff",
		nu()/pr_ + nut_/prt_
	);

	tmp<fvScalarMatrix> TEqn
	(
		fvm::ddt(T_)
	  + fvm::div(phi_, T_)
	  - fvm::Sp(fvc::div(phi_), T_)
	  - fvm::laplacian(kappaEff, T_)         
	);
    
	if(solveT_ == "true")	
	{
		TEqn().relax();
		solve(TEqn);
	}
	
	
    //*************************************//   
    // Output some min/max debug values
    //*************************************//
	
	if(sMMdebug == "true")
	{
    
    volScalarField phiActual("phiActual",tpphi_*k_);
	volVectorField psiActual("psiActual",tppsi_*k_);
	volScalarField uTauSquared((nu() + nut_)*vorticity_.component(2));
	volVectorField addedPsiG(addedPsiProd*k_);
	volVectorField phiOmega(tpphi_*k_*vorticity_);
	volVectorField prodDiff(vecProd - phiOmega); 
	volScalarField moreProd(upsilon_*mag((funit_ & uGrad_) & funit_));
	volVectorField transTerm(0.0033*sqrt((((nu()/100.0)+nut_)/nu()))*vorticity_*k_);
	

	Info << "***************************" << endl;
	Info << "Max addedPsiG: " << gMax(addedPsiG) << " Min addedPsiG: " << gMin(addedPsiG) << endl;
	Info << "Max transTerm: " << gMax(transTerm) << " Min transTerm: " << gMin(transTerm) << endl;
	Info << "Max moreProd: " << gMax(moreProd) << " Min moreProd: " << gMin(moreProd) << endl;
	Info << "Max phiProd: " << gMax(phiProd_) << " Min phiProd: " << gMin(phiProd_) << endl;
	Info << "Max phiOmega: " << gMax(phiOmega) << " Min phiOmega: " << gMin(phiOmega) << endl;
	Info << "Max vecProd: " << gMax(vecProd) << " Min vecProd: " << gMin(vecProd) << endl;
	Info << "Max prodDiff: " << gMax(prodDiff) << " Min prodDiff: " << gMin(prodDiff) << endl;
	Info << "Max cEps1: " << gMax(cEp1eqn) << " Min cEps1: " << gMin(cEp1eqn) << endl;
	Info<< "Max f: " << gMax(f_) << " Min f: " << gMin(f_) << " Max G: " << gMax(G) << " Max Gnut: " << gMax(Gnut) << endl;
    Info<< "Max nut: " << gMax(nut_) << " Max K: " << gMax(k_) << " Max Epsilon: " << gMax(epsilon_) << " Max Phi: " << gMax(phiActual) <<endl;
    Info<< "Max Psi: " << gMax(psiActual) << " Min Psi: " << gMin(psiActual) << endl;
	//Info<< "XX: " << tensor::XX << " YY: " << tensor::YY << " XY: " << tensor::XY << endl;
    Info<< "Max 3D Production: " << gMax(tpProd3d_) << " Max uTauSquared: " << gMax(uTauSquared) << endl;
	Info<< "Max vorticity: " << gMax(vorticity_) << " Min vorticity: " << gMin(vorticity_) << endl;
	Info << "***************************" << endl;
	}
	
	

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
