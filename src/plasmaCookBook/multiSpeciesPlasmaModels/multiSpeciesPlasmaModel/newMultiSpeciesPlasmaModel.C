/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.
\*---------------------------------------------------------------------------*/

#include "multiSpeciesPlasmaModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::multiSpeciesPlasmaModel>
Foam::multiSpeciesPlasmaModel::New
(
    hsCombustionThermo& thermo
)
{
    word modelName;
    {
        IOdictionary dict
        (
            IOobject
            (
                "plasmaProperties", 
                thermo.T().mesh().time().constant(),
                thermo.T().mesh(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );
        dict.lookup("plasmaModel") >> modelName;
    }
  
    fvMeshConstructorTable::iterator cstrIter =
        fvMeshConstructorTablePtr_->find(modelName);

    if (cstrIter == fvMeshConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "plasmaModel"
        )   << "Unknown type "
            << modelName << endl << endl
            << fvMeshConstructorTablePtr_->toc()
            << abort(FatalError);
  	}

  return autoPtr<multiSpeciesPlasmaModel>
      (cstrIter()(thermo));
}

// ************************************************************************* //
