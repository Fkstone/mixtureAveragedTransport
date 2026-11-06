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

\*---------------------------------------------------------------------------*/

#include "mixtureAveragedFourier.H"
#include "fvcDiv.H"
#include "fvcLaplacian.H"
#include "fvcSnGrad.H"
#include "fvmSup.H"
#include "surfaceInterpolate.H"
#include "Function2Evaluate.H"
#include "linear.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace laminarThermophysicalTransportModels
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class laminarThermophysicalTransportModel>
void mixtureAveragedFourier<laminarThermophysicalTransportModel>::
updateDm() const
{
    const PtrList<volScalarField>& Y = this->thermo().Y();
    const volScalarField& p = this->thermo().p();
    const volScalarField& T = this->thermo().T();

    Dm_.setSize(Y.size());

    // Mixture molecular weight at the old-time mass fractions
    volScalarField sumYbyW
    (
        volScalarField::New
        (
            "sumYbyW",
            T.mesh(),
            dimensionedScalar(dimless/this->thermo().Wi(0).dimensions(), 0)
        )
    );

    forAll(Y, j)
    {
        sumYbyW += Y[j].oldTime()/this->thermo().Wi(j);
    }

    W0_.reset(volScalarField::New("W0", 1/sumYbyW).ptr());

    volScalarField sumXbyD
    (
        volScalarField::New
        (
            "sumXbyD",
            T.mesh(),
            dimless/dimKinematicViscosity/W0_().dimensions()
        )
    );

    // Thermal diffusivity for the unityLewisSpecies
    tmp<volScalarField> tDth;
    if (findIndex(unityLewis_, true) != -1)
    {
        tDth = volScalarField::New
        (
            "Dth",
            this->thermo().kappa()
           /(this->thermo().Cp()*this->momentumTransport().rho())
        );
    }

    forAll(Dm_, i)
    {
        if (unityLewis_[i])
        {
            Dm_.set(i, volScalarField::New("Dth", tDth()));
            continue;
        }

        sumXbyD = Zero;

        forAll(Y, j)
        {
            if (j != i)
            {
                sumXbyD +=
                    Y[j].oldTime()
                   /(
                        this->thermo().Wi(j)
                       *(
                            i < j
                          ? evaluate
                            (
                                DFuncs_[i][j],
                                dimKinematicViscosity,
                                p,
                                T
                            )
                          : evaluate
                            (
                                DFuncs_[j][i],
                                dimKinematicViscosity,
                                p,
                                T
                            )
                       )
                   );
            }
        }

        // sumXbyD = sum_j Y_j/(W_j D_ij) = sum_j (X_j/D_ij)/W, so
        // D_i = (1 - Y_i)/sum_j (X_j/D_ij) = (1 - Y_i)/(W sumXbyD)
        Dm_.set
        (
            i,
            (1 - Y[i].oldTime())
           /(
               W0_()
              *max(sumXbyD, dimensionedScalar(sumXbyD.dimensions(), small))
            )
        );
    }
}


template<class laminarThermophysicalTransportModel>
const PtrList<volScalarField>&
mixtureAveragedFourier<laminarThermophysicalTransportModel>::Dm() const
{
    if (!Dm_.size())
    {
        updateDm();
    }

    return Dm_;
}


template<class laminarThermophysicalTransportModel>
void mixtureAveragedFourier<laminarThermophysicalTransportModel>::
updateSums() const
{
    const PtrList<volScalarField>& Y = this->thermo().Y();
    const PtrList<volScalarField>& D = Dm();
    const volScalarField& W0 = W0_();
    const fvMesh& mesh = W0.mesh();
    const linear<scalar> interp(mesh);

    // grad(ln W) on the faces
    snGradLnW_.reset
    (
        surfaceScalarField::New
        (
            "snGradLnW",
            fvc::snGrad(W0)/interp.interpolate(W0)
        ).ptr()
    );

    // Sum of the uncorrected mixture-averaged fluxes of all species
    // (including the default specie)
    sumJ_.reset
    (
        surfaceScalarField::New
        (
            "sumJ",
            mesh,
            dimensionedScalar(dimMass/dimArea/dimTime, 0)
        ).ptr()
    );

    forAll(Y, i)
    {
        sumJ_() += jMA(i);
    }

    // sum_k Y_k D_k
    sumYD_.reset
    (
        volScalarField::New
        (
            "sumYD",
            mesh,
            dimensionedScalar(D[0].dimensions(), 0)
        ).ptr()
    );

    forAll(Y, k)
    {
        sumYD_() += Y[k].oldTime()*D[k];
    }
}


template<class laminarThermophysicalTransportModel>
void mixtureAveragedFourier<laminarThermophysicalTransportModel>::
update() const
{
    if (!Dm_.size())
    {
        updateDm();
    }

    if (!sumJ_.valid())
    {
        updateSums();
    }
}


template<class laminarThermophysicalTransportModel>
tmp<surfaceScalarField>
mixtureAveragedFourier<laminarThermophysicalTransportModel>::jMA
(
    const label i
) const
{
    const volScalarField& Y0i = this->thermo().Y()[i].oldTime();
    const volScalarField& rho = this->momentumTransport().rho();
    const linear<scalar> interp(Y0i.mesh());

    return surfaceScalarField::New
    (
        "jMA" + this->thermo().Y()[i].name(),
       -interp.interpolate(this->alpha()*rho*Dm()[i])
       *(fvc::snGrad(Y0i) + interp.interpolate(Y0i)*snGradLnW_())
    );
}


template<class laminarThermophysicalTransportModel>
tmp<volScalarField>
mixtureAveragedFourier<laminarThermophysicalTransportModel>::Gamma
(
    const label i
) const
{
    update();

    const PtrList<volScalarField>& D = Dm_;
    const label d = this->thermo().defaultSpecie();
    const word name("Gamma" + this->thermo().Y()[i].name());

    if (i == d)
    {
        return volScalarField::New(name, D[i]);
    }
    else
    {
        const volScalarField& Y0i = this->thermo().Y()[i].oldTime();

        return volScalarField::New
        (
            name,
            max
            (
                D[i],
                D[i]*(1 - Y0i) + Y0i*D[d]
              - Y0i*W0_()*(1/this->thermo().Wi(i) - 1/this->thermo().Wi(d))
               *(D[i] - sumYD_())
            )
        );
    }
}


template<class laminarThermophysicalTransportModel>
tmp<scalarField>
mixtureAveragedFourier<laminarThermophysicalTransportModel>::Gamma
(
    const label i,
    const label patchi
) const
{
    update();

    const label d = this->thermo().defaultSpecie();
    const scalarField& Di = Dm_[i].boundaryField()[patchi];

    if (i == d)
    {
        return tmp<scalarField>(new scalarField(Di));
    }
    else
    {
        const scalarField& Y0i =
            this->thermo().Y()[i].oldTime().boundaryField()[patchi];

        return max
        (
            Di,
            Di*(1 - Y0i) + Y0i*Dm_[d].boundaryField()[patchi]
          - Y0i*W0_().boundaryField()[patchi]
           *(1/this->thermo().Wi(i).value() - 1/this->thermo().Wi(d).value())
           *(Di - sumYD_().boundaryField()[patchi])
        );
    }
}


template<class laminarThermophysicalTransportModel>
tmp<surfaceScalarField>
mixtureAveragedFourier<laminarThermophysicalTransportModel>::jexp
(
    const label i
) const
{
    update();

    // Corrected mixture-averaged flux plus the implicit flux, both at the
    // old-time state
    const volScalarField& Y0i = this->thermo().Y()[i].oldTime();
    const volScalarField& rho = this->momentumTransport().rho();
    const linear<scalar> interp(Y0i.mesh());

    tmp<surfaceScalarField> tjexp(jMA(i));

    tjexp.ref() +=
      - interp.interpolate(Y0i)*sumJ_()
      + interp.interpolate(this->alpha()*rho*Gamma(i))*fvc::snGrad(Y0i);

    tjexp.ref().rename("jexp" + this->thermo().Y()[i].name());

    return tjexp;
}


template<class laminarThermophysicalTransportModel>
void mixtureAveragedFourier<laminarThermophysicalTransportModel>::
clear() const
{
    Dm_.clear();
    W0_.clear();
    sumYD_.clear();
    snGradLnW_.clear();
    sumJ_.clear();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class laminarThermophysicalTransportModel>
mixtureAveragedFourier<laminarThermophysicalTransportModel>::
mixtureAveragedFourier
(
    const momentumTransportModel& momentumTransport,
    const thermoModel& thermo
)
:
    unityLewisFourier<laminarThermophysicalTransportModel>
    (
        typeName,
        momentumTransport,
        thermo
    ),

    TopoChangeableMeshObject(*this),

    DFuncs_(this->thermo().species().size()),

    unityLewis_(this->thermo().species().size(), false)
{
    read();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class laminarThermophysicalTransportModel>
bool mixtureAveragedFourier<laminarThermophysicalTransportModel>::read()
{
    if (unityLewisFourier<laminarThermophysicalTransportModel>::read())
    {
        const speciesTable& species = this->thermo().species();

        const dictionary& coeffDict = this->typeDict();

        const dictionary& Ddict = coeffDict.subDict("D");

        // Species with unity Lewis number
        unityLewis_ = false;
        const wordList uLe
        (
            coeffDict.lookupOrDefault<wordList>("unityLewisSpecies", wordList())
        );
        forAll(uLe, k)
        {
            if (!species.found(uLe[k]))
            {
                FatalIOErrorInFunction(coeffDict)
                    << "unityLewisSpecies: unknown specie " << uLe[k]
                    << exit(FatalIOError);
            }
            unityLewis_[species[uLe[k]]] = true;
            Info<< "    " << uLe[k]
                << ": unity Lewis number, D = kappa/(rho Cp)" << endl;
        }

        // Read the array of specie binary mass diffusion coefficient
        // functions
        forAll(species, i)
        {
            DFuncs_[i].setSize(species.size());

            forAll(species, j)
            {
                if (j >= i)
                {
                    const word nameij(species[i] + '-' + species[j]);
                    const word nameji(species[j] + '-' + species[i]);

                    word Dname;

                    if (Ddict.found(nameij) && Ddict.found(nameji))
                    {
                        if (i != j)
                        {
                            WarningInFunction
                                << "Binary mass diffusion coefficients "
                                   "for Both " << nameij
                                << " and " << nameji << " provided, using "
                                << nameij << endl;
                        }

                        Dname = nameij;
                    }
                    else if (Ddict.found(nameij))
                    {
                        Dname = nameij;
                    }
                    else if (Ddict.found(nameji))
                    {
                        Dname = nameji;
                    }
                    else
                    {
                        FatalIOErrorInFunction(Ddict)
                            << "Binary mass diffusion coefficient for pair "
                            << nameij << " or " << nameji << " not provided"
                            << exit(FatalIOError);
                    }

                    DFuncs_[i].set
                    (
                        j,
                        Function2<scalar>::New
                        (
                            Dname,
                            dimPressure,
                            dimTemperature,
                            dimKinematicViscosity,
                            Ddict
                        ).ptr()
                    );
                }
            }
        }

        return true;
    }
    else
    {
        return false;
    }
}


template<class laminarThermophysicalTransportModel>
tmp<volScalarField>
mixtureAveragedFourier<laminarThermophysicalTransportModel>::D
(
    const volScalarField& Yi
) const
{
    return volScalarField::New
    (
        "D",
        this->momentumTransport().rho()
       *Gamma(this->thermo().specieIndex(Yi))
    );
}


template<class laminarThermophysicalTransportModel>
tmp<scalarField>
mixtureAveragedFourier<laminarThermophysicalTransportModel>::D
(
    const volScalarField& Yi,
    const label patchi
) const
{
    return
        this->momentumTransport().rho().boundaryField()[patchi]
       *Gamma(this->thermo().specieIndex(Yi), patchi);
}


template<class laminarThermophysicalTransportModel>
tmp<volScalarField>
mixtureAveragedFourier<laminarThermophysicalTransportModel>::DEff
(
    const volScalarField& Yi
) const
{
    return volScalarField::New
    (
        "DEff",
        this->momentumTransport().rho()
       *Gamma(this->thermo().specieIndex(Yi))
    );
}


template<class laminarThermophysicalTransportModel>
tmp<scalarField>
mixtureAveragedFourier<laminarThermophysicalTransportModel>::DEff
(
    const volScalarField& Yi,
    const label patchi
) const
{
    return
        this->momentumTransport().rho().boundaryField()[patchi]
       *Gamma(this->thermo().specieIndex(Yi), patchi);
}


template<class laminarThermophysicalTransportModel>
tmp<surfaceScalarField>
mixtureAveragedFourier<laminarThermophysicalTransportModel>::q() const
{
    tmp<surfaceScalarField> tmpq
    (
        surfaceScalarField::New
        (
            IOobject::groupName
            (
                "q",
                this->thermo().phaseName()
            ),
           -fvc::interpolate(this->alpha()*this->kappaEff())
           *fvc::snGrad(this->thermo().T())
        )
    );

    const PtrList<volScalarField>& Y = this->thermo().Y();
    const volScalarField& p = this->thermo().p();
    const volScalarField& T = this->thermo().T();

    if (Y.size())
    {
        surfaceScalarField sumJ
        (
            surfaceScalarField::New
            (
                "sumJ",
                Y[0].mesh(),
                dimensionedScalar(dimMass/dimArea/dimTime, 0)
            )
        );

        surfaceScalarField sumJh
        (
            surfaceScalarField::New
            (
                "sumJh",
                Y[0].mesh(),
                dimensionedScalar(sumJ.dimensions()*dimEnergy/dimMass, 0)
            )
        );

        forAll(Y, i)
        {
            if (i != this->thermo().defaultSpecie())
            {
                const volScalarField hi(this->thermo().hsi(i, p, T));

                const surfaceScalarField ji(this->j(Y[i]));

                sumJ += ji;

                sumJh += ji*fvc::interpolate(hi);
            }
        }

        {
            const label i = this->thermo().defaultSpecie();

            const volScalarField hi(this->thermo().hsi(i, p, T));

            sumJh -= sumJ*fvc::interpolate(hi);
        }

        tmpq.ref() += sumJh;
    }

    return tmpq;
}


template<class laminarThermophysicalTransportModel>
tmp<scalarField>
mixtureAveragedFourier<laminarThermophysicalTransportModel>::q
(
    const label patchi
) const
{
    tmp<scalarField> tmpq
    (
      - (
            this->alpha().boundaryField()[patchi]
           *this->kappaEff(patchi)
           *this->thermo().T().boundaryField()[patchi].snGrad()
        )
    );

    const PtrList<volScalarField>& Y = this->thermo().Y();
    const scalarField& p = this->thermo().p().boundaryField()[patchi];
    const scalarField& T = this->thermo().T().boundaryField()[patchi];

    if (Y.size())
    {
        scalarField sumJ(tmpq->size(), scalar(0));
        scalarField sumJh(tmpq->size(), scalar(0));

        forAll(Y, i)
        {
            if (i != this->thermo().defaultSpecie())
            {
                const scalarField hi(this->thermo().hsi(i, p, T));

                const scalarField ji(this->j(Y[i], patchi));

                sumJ += ji;

                sumJh += ji*hi;
            }
        }

        {
            const label i = this->thermo().defaultSpecie();

            const scalarField hi(this->thermo().hsi(i, p, T));

            sumJh -= sumJ*hi;
        }

        tmpq.ref() += sumJh;
    }

    return tmpq;
}


template<class laminarThermophysicalTransportModel>
tmp<fvScalarMatrix>
mixtureAveragedFourier<laminarThermophysicalTransportModel>::divq
(
    volScalarField& he
) const
{
    tmp<fvScalarMatrix> tmpDivq
    (
        fvm::Su
        (
            -fvc::laplacian(this->alpha()*this->kappaEff(), this->thermo().T()),
            he
        )
    );

    const PtrList<volScalarField>& Y = this->thermo().Y();
    const volScalarField& p = this->thermo().p();
    const volScalarField& T = this->thermo().T();

    tmpDivq.ref() -=
        fvm::laplacianCorrection(this->alpha()*this->alphaEff(), he);

    surfaceScalarField sumJ
    (
        surfaceScalarField::New
        (
            "sumJ",
            he.mesh(),
            dimensionedScalar(dimMass/dimArea/dimTime, 0)
        )
    );

    surfaceScalarField sumJh
    (
        surfaceScalarField::New
        (
            "sumJh",
            he.mesh(),
            dimensionedScalar(sumJ.dimensions()*he.dimensions(), 0)
        )
    );

    forAll(Y, i)
    {
        if (i != this->thermo().defaultSpecie())
        {
            const volScalarField hi(this->thermo().hsi(i, p, T));

            const surfaceScalarField ji(this->j(Y[i]));

            sumJ += ji;

            sumJh += ji*fvc::interpolate(hi);
        }
    }

    {
        const label i = this->thermo().defaultSpecie();

        const volScalarField hi(this->thermo().hsi(i, p, T));

        sumJh -= sumJ*fvc::interpolate(hi);
    }

    tmpDivq.ref() += fvc::div(sumJh*he.mesh().magSf());

    return tmpDivq;
}


template<class laminarThermophysicalTransportModel>
tmp<surfaceScalarField>
mixtureAveragedFourier<laminarThermophysicalTransportModel>::j
(
    const volScalarField& Yi
) const
{
    return
        unityLewisFourier<laminarThermophysicalTransportModel>::j(Yi)
      + jexp(this->thermo().specieIndex(Yi));
}


template<class laminarThermophysicalTransportModel>
tmp<scalarField>
mixtureAveragedFourier<laminarThermophysicalTransportModel>::j
(
    const volScalarField& Yi,
    const label patchi
) const
{
    return
        unityLewisFourier<laminarThermophysicalTransportModel>::j(Yi, patchi)
      + jexp(this->thermo().specieIndex(Yi))().boundaryField()[patchi];
}


template<class laminarThermophysicalTransportModel>
tmp<fvScalarMatrix>
mixtureAveragedFourier<laminarThermophysicalTransportModel>::divj
(
    volScalarField& Yi
) const
{
    return
        unityLewisFourier<laminarThermophysicalTransportModel>::divj(Yi)
      + fvc::div
        (
            jexp(this->thermo().specieIndex(Yi))*Yi.mesh().magSf()
        );
}


template<class laminarThermophysicalTransportModel>
void mixtureAveragedFourier<laminarThermophysicalTransportModel>::predict()
{
    unityLewisFourier<laminarThermophysicalTransportModel>::predict();
    updateDm();
    updateSums();
}


template<class laminarThermophysicalTransportModel>
bool mixtureAveragedFourier<laminarThermophysicalTransportModel>::movePoints()
{
    return true;
}


template<class laminarThermophysicalTransportModel>
void mixtureAveragedFourier<laminarThermophysicalTransportModel>::topoChange
(
    const polyTopoChangeMap& map
)
{
    clear();
}


template<class laminarThermophysicalTransportModel>
void mixtureAveragedFourier<laminarThermophysicalTransportModel>::mapMesh
(
    const polyMeshMap& map
)
{
    clear();
}


template<class laminarThermophysicalTransportModel>
void mixtureAveragedFourier<laminarThermophysicalTransportModel>::distribute
(
    const polyDistributionMap& map
)
{
    clear();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace laminarThermophysicalTransportModels
} // End namespace Foam

// ************************************************************************* //
