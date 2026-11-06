/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2026 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is derived from OpenFOAM (Fickian, FickianFourier).

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::laminarThermophysicalTransportModels::mixtureAveragedFourier

Description
    Mixture-averaged (Hirschfelder-Curtiss) specie diffusion with a
    correction velocity and Fourier based temperature gradient heat flux for
    laminar flow, as the mixture-averaged model of Cantera/Chemkin:

        j_i = -rho D_i (W_i/W) grad(X_i) + Y_i sum_k rho D_k (W_k/W) grad(X_k)

        D_i = (1 - Y_i)/sum_{j != i} (X_j/D_ij)

    with (W_i/W) grad(X_i) = grad(Y_i) + Y_i grad(ln W).

    Time discretisation (deferred correction, as the explicit/implicit split
    of MaxwellStefan): with Gamma_i an implicit diffusivity,

        j_i = -rho Gamma_i grad(Y_i) + [j_i^MA + rho Gamma_i grad(Y_i)]^old

    where the bracket is evaluated at the old-time mass fractions Y_i^0 and
    the D_i of the last predict() (called once at the start of each time
    step), so that j_i = j_i^MA at the start of the step. To keep the memory
    at that of FickianFourier only the D_i and four fields shared by all
    species (W, sum_k Y_k D_k, snGrad(ln W) and the sum of the uncorrected
    fluxes) are cached; Gamma_i and the bracket are evaluated per specie
    when needed. Gamma_i =
    max(D_i, Dii_i) with Dii_i the derivative -dj_i/d(rho grad(Y_i)) of the
    mixture-averaged flux at fixed gradients of the other non-default
    species (the default specie closes the mass fractions):

        Dii_i = D_i (1 - Y_i) + Y_i D_d + Y_i g_i (D_i - sum_k Y_k D_k),
        g_i = -W (1/W_i - 1/W_d)

    which reduces to the binary D for two species. The time-lagged part
    then contains no first-order dependence on grad(Y_i), which makes the
    scheme stable for Dii_i > 0 at any D dt/dx^2. Treating the grad(ln W)
    term and the correction velocity simply as explicit (or as convection
    with a lagged correction flux) is unstable at D dt/dx^2 >> 1: for a
    dominant specie they cancel most of the implicit -rho D_i grad(Y_i).

    Differences to FickianFourier with mixtureDiffusionCoefficients no:
      - D_i uses (1 - Y_i) instead of (1 - X_i) in the numerator
      - the Y_i grad(ln W) term is included
      - the sum of the specie fluxes is removed by the correction velocity
        instead of being absorbed by the default specie

    Soret diffusion is not included.

    Optional unityLewisSpecies: for the listed species the mixture
    diffusion coefficient is replaced by the thermal diffusivity,
    D_i = kappa/(rho Cp) (Le_i = 1 in the mixture-averaged sense, i.e. with
    the (W_i/W) grad(X_i) driving force and the correction velocity as for
    the other species; Cantera's Le_k = alpha/D_k,mix); all other species
    keep the mixture-averaged D_i. The binary coefficients of these species
    are still needed for the D_j of the others. Default: none.

Usage
    \verbatim
    laminar
    {
        model           mixtureAveragedFourier;

        unityLewisSpecies (nc7h16); // optional

        D // [m^2/s] binary mass diffusion coefficient functions
        {
            O2-O2 1e-2;
            O2-N2 1e-2;
            N2-N2 1e-2;
        }
    }
    \endverbatim

SourceFiles
    mixtureAveragedFourier.C

\*---------------------------------------------------------------------------*/

#ifndef mixtureAveragedFourier_H
#define mixtureAveragedFourier_H

#include "unityLewisFourier.H"
#include "Function2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace laminarThermophysicalTransportModels
{

/*---------------------------------------------------------------------------*\
                Class mixtureAveragedFourier Declaration
\*---------------------------------------------------------------------------*/

template<class laminarThermophysicalTransportModel>
class mixtureAveragedFourier
:
    public unityLewisFourier<laminarThermophysicalTransportModel>,
    public TopoChangeableMeshObject<fvMesh>
{
    // Private data

        //- Array of specie binary mass diffusion coefficient functions
        //  [m^2/s]
        List<PtrList<Function2<scalar>>> DFuncs_;

        //- Species with D_i = kappa/(rho Cp) (unityLewisSpecies)
        boolList unityLewis_;

        //- Mass diffusion coefficients in the mixture D_i [m^2/s]
        mutable PtrList<volScalarField> Dm_;

        //- Mixture molecular weight at the old-time mass fractions [kg/kmol]
        mutable autoPtr<volScalarField> W0_;

        //- sum_k Y_k D_k at the old-time mass fractions [m^2/s]
        mutable autoPtr<volScalarField> sumYD_;

        //- snGrad(ln W) at the old-time mass fractions [1/m]
        mutable autoPtr<surfaceScalarField> snGradLnW_;

        //- Sum of the uncorrected mixture-averaged fluxes [kg/m^2/s]
        mutable autoPtr<surfaceScalarField> sumJ_;


    // Private Member Functions

        //- Update Dm_ and W0_ from the old-time mass fractions
        void updateDm() const;

        //- Update Dm_ if not yet constructed and return
        const PtrList<volScalarField>& Dm() const;

        //- Update sumYD_, snGradLnW_ and sumJ_ from Dm_ and the old-time
        //  mass fractions
        void updateSums() const;

        //- Update the cached fields if not yet constructed
        void update() const;

        //- Uncorrected mixture-averaged flux of specie i at the old-time
        //  mass fractions [kg/m^2/s]
        tmp<surfaceScalarField> jMA(const label i) const;

        //- Implicit diffusivity Gamma_i [m^2/s]
        tmp<volScalarField> Gamma(const label i) const;

        //- Implicit diffusivity Gamma_i for patch [m^2/s]
        tmp<scalarField> Gamma(const label i, const label patchi) const;

        //- Time-lagged part of the mass flux of specie i [kg/m^2/s]
        tmp<surfaceScalarField> jexp(const label i) const;

        //- Delete the cached fields, re-created in predict
        void clear() const;


public:

    typedef typename laminarThermophysicalTransportModel::alphaField
        alphaField;

    typedef typename
        laminarThermophysicalTransportModel::momentumTransportModel
        momentumTransportModel;

    typedef typename laminarThermophysicalTransportModel::thermoModel
        thermoModel;


    //- Runtime type information
    TypeName("mixtureAveragedFourier");


    // Constructors

        //- Construct from a momentum transport model and a thermo model
        mixtureAveragedFourier
        (
            const momentumTransportModel& momentumTransport,
            const thermoModel& thermo
        );


    //- Destructor
    virtual ~mixtureAveragedFourier()
    {}


    // Member Functions

        //- Read thermophysicalTransport dictionary
        virtual bool read();

        //- Mass diffusivity (implicit part, rho Gamma_i)
        //  for a given specie mass-fraction [kg/m/s]
        virtual tmp<volScalarField> D(const volScalarField& Yi) const;

        //- Mass diffusivity
        //  for a given specie mass-fraction for patch [kg/m/s]
        virtual tmp<scalarField> D
        (
            const volScalarField& Yi,
            const label patchi
        ) const;

        //- Effective mass diffusion coefficient
        //  for a given specie mass-fraction [kg/m/s]
        virtual tmp<volScalarField> DEff(const volScalarField& Yi) const;

        //- Effective mass diffusion coefficient
        //  for a given specie mass-fraction for patch [kg/m/s]
        virtual tmp<scalarField> DEff
        (
            const volScalarField& Yi,
            const label patchi
        ) const;

        //- Return the heat flux [W/m^2]
        virtual tmp<surfaceScalarField> q() const;

        //- Return the patch heat flux [W/m^2]
        virtual tmp<scalarField> q(const label patchi) const;

        //- Return the source term for the energy equation
        virtual tmp<fvScalarMatrix> divq(volScalarField& he) const;

        //- Return the specie flux for the given specie mass-fraction [kg/m^2/s]
        virtual tmp<surfaceScalarField> j(const volScalarField& Yi) const;

        //- Return the specie flux
        //  for the given specie mass-fraction for patch [kg/m^2/s]
        virtual tmp<scalarField> j
        (
            const volScalarField& Yi,
            const label patchi
        ) const;

        //- Return the source term for the given specie mass-fraction equation
        virtual tmp<fvScalarMatrix> divj(volScalarField& Yi) const;

        //- Update the diffusion coefficients and the time-lagged fluxes
        virtual void predict();


        // Mesh changes

            //- Update for mesh motion
            virtual bool movePoints();

            //- Update topology using the given map
            virtual void topoChange(const polyTopoChangeMap& map);

            //- Update from another mesh using the given map
            virtual void mapMesh(const polyMeshMap& map);

            //- Redistribute or update using the given distribution map
            virtual void distribute(const polyDistributionMap& map);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace laminarThermophysicalTransportModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "mixtureAveragedFourier.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
