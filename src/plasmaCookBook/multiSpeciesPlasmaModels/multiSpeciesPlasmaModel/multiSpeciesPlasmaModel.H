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
#include <math.h>
#include "foamTime.H"
#include "speciesTable.H"
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
        
        scalar meshParameter();

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

        virtual bool read() = 0;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "multiSpeciesPlasmaModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
