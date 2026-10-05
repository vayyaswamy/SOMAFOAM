/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

Class
    Foam::multiSpeciesPlasmaModel

Description
    Abstract base class for multispecies plasma fluid models.

SourceFiles
    multiSpeciesPlasmaModelI.H
    multiSpeciesPlasmaModel.C

\*---------------------------------------------------------------------------*/

#ifndef multiSpeciesPlasmaModel_H
#define multiSpeciesPlasmaModel_H

#include "fvc.H"
#include "IOdictionary.H"
#include "psiChemistryModel.H"
#include "reactingMixture.H"
#include "hsCombustionThermo.H"
#include "IFstream.H"
#include "graph.H"
#include "plasmaInterpolateXY.H"
#include "dimensionedConstants.H"
#include "plasmaConstants.H"
#include "relaxFinal.H"
#include "scharfetterGummel.H"
#include <math.h>
#include "foamTime.H"
#include "speciesTable.H"
#include "HashTable.H"
#include <algorithm>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class multiSpeciesPlasmaModel Declaration
\*---------------------------------------------------------------------------*/

class multiSpeciesPlasmaModel
:
    public IOdictionary
{

protected:

        hsCombustionThermo& thermo_;

        const fvMesh& mesh_;

        Time& runTime_;

        const int activeSpecies_;

        const bool multiTimeStep;

        const scalar STScounter;

        const scalar MTScounter;

        const scalar LTScounter;

        const word eSpecie;

        const label eIndex_;

        const word bGas;

        const label bIndex_;

        word plasmaChemistry;

        Switch restartcapabale;

        Switch eonCalculation;

        PtrList<graph> graphDiffData_;

        PtrList<graph> graphMuData_;

        PtrList<scalar> z_;

        PtrList<word> diffusionModel_;

        PtrList<word> mobilityModel_;

        PtrList<word> transportModel_;

        PtrList<word> collisionFrequency_;

        PtrList<Switch> speciesSolution_;

        PtrList<word> timeS_;

        PtrList<volScalarField> Sy_;

        PtrList<volScalarField> D_;

        PtrList<volScalarField> mu_;

        PtrList<volScalarField> T_;

        PtrList<volVectorField> U_;

        PtrList<volVectorField> F_;

        PtrList<volVectorField> J_;

        PtrList<volScalarField> N_;

        PtrList<volScalarField> collFrequency_;

        volScalarField NG_;

        // Inner (nonlinear) iterations per time step, optional dictionary
        // innerIterations in plasmaProperties

            //- Maximum number of passes for charged species, neutral
            //  species and the electron temperature
            label nCorrCharged_;
            label nCorrNeutral_;
            label nCorrTe_;

            //- Initial residual below which the passes stop
            scalar innerTolerance_;

            //- Report loops that end without reaching the tolerance
            Switch reportUnconverged_;

        //- Per species: Scharfetter-Gummel fluxes instead of the convection
        //  and Laplacian schemes of fvSchemes. Optional entry fluxScheme in
        //  plasmaProperties, and in the dictionary of a species for that
        //  species alone
        boolList scharfetterGummel_;

        //- Fluxes of the drift-diffusion species through the faces [1/s],
        //  as discretised in their equations; set after each species has
        //  been solved. The electron energy equation and the fluxes in the
        //  cells use them
        PtrList<surfaceScalarField> faceFlux_;

        // Limits, optional dictionary limits in plasmaProperties

            //- Lower limit of the number densities of the transported
            //  species [1/m3]
            scalar densityFloor_;

            //- Lower and upper limit of the electron temperature [K]
            scalar TeMin_;
            scalar TeMax_;

            //- Lower and upper limits of individual fields, by field name
            HashTable<scalar> fieldMin_;
            HashTable<scalar> fieldMax_;

            //- Modification time of plasmaProperties when the limits were
            //  read
            time_t limitsTime_;

        // Acceleration of slow neutral species (e.g. metastables), optional
        // dictionary slowSpeciesAcceleration in plasmaProperties

            Switch accelerate_;

            //- Indices of the accelerated species
            labelList accSpecies_;

            //- Period of the applied voltage [s]
            scalar accPeriod_;

            //- Time after which the acceleration is active [s]
            scalar accStartTime_;

            //- Periods of full simulation between two advances, and periods
            //  at the end of them over which the sources are averaged
            label accFullCycles_;
            label accAverageCycles_;

            //- Solve for the steady state of the averaged equation instead
            //  of advancing in time
            Switch accSteady_;

            //- Time step and number of steps of the advance
            scalar accDeltaT_;
            label accNSteps_;

            //- Fraction of the computed change that is applied
            scalar accRelaxation_;

            //- Largest factor by which the density may change in one advance
            scalar accMaxChangeFactor_;

            //- Relative change per advance below which the advances stop
            scalar accTolerance_;

            //- Start time of the current block, averaged time so far
            scalar accBlockStart_;
            scalar accAveragedTime_;
            scalar accLastTime_;
            bool accConverged_;

            //- Time integrals of the explicit and implicit source parts
            PtrList<volScalarField> accSu_;
            PtrList<volScalarField> accSp_;

        //- Read the innerIterations and slowSpeciesAcceleration dictionaries
        void readNumericalControls();

        //- Read the limits dictionary of plasmaProperties
        void readLimits(const dictionary& properties);

        //- Read the limits again if plasmaProperties has been modified
        //  since they were last read
        void readLimitsIfModified();

        //- Average the sources of the slow species over the applied period
        //  and advance these species alone with a large time step
        void accelerateSlowSpecies(psiChemistryModel& chemistry);


private:

        multiSpeciesPlasmaModel(const multiSpeciesPlasmaModel&);

        void operator=(const multiSpeciesPlasmaModel&);

        tmp<Foam::volScalarField> RR
        (
            const word Yname,
            const psiChemistryModel& chemistry,
            const label i
        ) const;

public:

    TypeName("multiSpeciesPlasmaModel");


#ifndef SWIG
    declareRunTimeSelectionTable
    (
        autoPtr,
        multiSpeciesPlasmaModel,
        fvMesh,
        (
            hsCombustionThermo& thermo
        ),
        (thermo)
    );
#endif

    // Constructors

        //- Construct from components
        multiSpeciesPlasmaModel
        (
            hsCombustionThermo& thermo
        );

    // Selectors

        //- Return a reference to the selected multispeciesplasma model

        static autoPtr<multiSpeciesPlasmaModel> New
        (
            hsCombustionThermo& thermo
        );

    //- Destructor

        virtual ~multiSpeciesPlasmaModel()
        {}

    // Member Functions

        // Access

        inline void LTSset();

        inline void MTSset();

        inline void STSset();

        inline const speciesTable& species() const;

        inline const volScalarField& N(const label i) const;

        inline const volScalarField& N(const word& specieName) const;

        inline const volVectorField& U(const label i) const;

        inline const volVectorField& J(const label i) const;

        inline const volVectorField& J(const word& specieName) const;

        inline const volScalarField& mu(const label i) const;

        virtual scalar W(label i) const = 0;

        inline void updateTemperature();

        void updateChemistry(psiChemistryModel& chemistry);

        inline void updateChemistryCollFreq(psiChemistryModel& chemistry);

        tmp<volScalarField> electronTempSource
        (
            const psiChemistryModel& chemistry
        );

        tmp<volScalarField> dElectronTempSourceDTe
        (
            const psiChemistryModel& chemistry
        );

        tmp<volScalarField> ionTempSource
        (
            const psiChemistryModel& chemistry,
            const volVectorField& E
        );

        scalar divFe();

        surfaceScalarField meshParameter();

        tmp<volVectorField> gradTe();

        tmp<volScalarField> netCharge();

        tmp<volScalarField> potentialImpSource();

        tmp<volScalarField> potentialExpSource();

        tmp<volVectorField> netChargeFlux();

        tmp<volVectorField> totalIonFlux();

        tmp<volVectorField> electronConvectiveFlux();

        tmp<volScalarField> electronConductivity(const psiChemistryModel& chemistry);

        virtual scalar correct
        (
            psiChemistryModel& chemistry,
            const volVectorField& E,
            multivariateSurfaceInterpolationScheme<scalar>::fieldTable& fields
        ) = 0;

        scalar correct
        (
            PtrList<volScalarField>& Y,
            const volVectorField& E,
            psiChemistryModel& chemistry,
            multivariateSurfaceInterpolationScheme<scalar>::fieldTable& fields
        );

        virtual void input();

        // Inner iteration controls

            label nCorrCharged() const
            {
                return nCorrCharged_;
            }

            label nCorrNeutral() const
            {
                return nCorrNeutral_;
            }

            label nCorrTe() const
            {
                return nCorrTe_;
            }

            scalar innerTolerance() const
            {
                return innerTolerance_;
            }

        // Flux scheme

            //- Does species i use Scharfetter-Gummel fluxes?
            bool scharfetterGummel(const label i) const
            {
                return scharfetterGummel_[i];
            }

            //- Implicit drift and diffusion terms of the equation of
            //  species i, div(F Ni) - laplacian(D, Ni), for a drift
            //  velocity F
            tmp<fvScalarMatrix> driftDiffusionTerms
            (
                const label i,
                const volVectorField& F,
                const volScalarField& D,
                const volScalarField& Ni
            ) const;

            //- Store the flux of species i through the faces, as
            //  discretised in its equation, and reconstruct its flux J in
            //  the cells from it, so that the electron energy equation and
            //  the currents use the flux of the species equation. Called
            //  after species i has been solved and J has been set.
            void correctDriftDiffusionFlux
            (
                const label i,
                const volVectorField& F,
                const volScalarField& D,
                const volScalarField& Ni
            );

            //- Is the flux of species i through the faces available? Only
            //  for drift-diffusion species
            bool hasFaceFlux(const label i) const
            {
                return faceFlux_.set(i);
            }

            //- Flux of species i through the faces
            const surfaceScalarField& faceFlux(const label i) const
            {
                return faceFlux_[i];
            }

        // Limits

            scalar densityFloor() const
            {
                return densityFloor_;
            }

            scalar TeMin() const
            {
                return TeMin_;
            }

            scalar TeMax() const
            {
                return TeMax_;
            }

            //- Apply the limits of a field: those given for it by name in
            //  the limits dictionary, and the defaults passed here otherwise
            void limitField
            (
                volScalarField& field,
                const scalar defaultMin = -GREAT,
                const scalar defaultMax = GREAT
            ) const;

            //- Report an inner loop that ended above the tolerance
            void reportInnerIterations
            (
                const word& name,
                const label maxPasses,
                const scalar initialResidual
            ) const;

        virtual bool read() = 0;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "multiSpeciesPlasmaModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
