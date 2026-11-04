/*---------------------------------------------------------------------------*\
Class
    Foam::convergenceStop

Description
    Optional stop of a run when it has reached a steady or periodic state,
    set by the subdictionary convergenceStop of system/controlDict:

    \verbatim
    convergenceStop
    {
        fields      (N_electron N_Arm Te);  // volScalarFields to test
        window      1.2e-7;     // averaging time [s]; whole RF periods for a
                                // periodic discharge
        tolerance   1e-4;       // largest relative change between windows
        nWindows    3;          // consecutive windows below the tolerance
                                // (optional, default 3)
        startTime   0;          // first window starts here (optional)
    }
    \endverbatim

    Each field is averaged in time over successive windows, with each step
    weighted by its length. At the end of a window its change from the
    previous window is measured, in each field, as the volume-weighted L2
    norm of the difference relative to that of the previous average. When the
    largest change over the fields stays below the tolerance for nWindows
    consecutive windows the run writes its fields and ends. Without the
    subdictionary nothing is computed and the run ends at endTime.

    The test measures the change per window, not the distance to the final
    state: a field relaxing with time constant tau changes by about
    (window/tau) times its remaining error in each window. Short windows and
    a loose tolerance therefore stop a slow transient early. Choose the
    tolerance as the accepted error times window/tau, and test the slowest
    field.

    A slowly relaxing species must be among the fields, and the window must
    be long enough for its change to show: with slowSpeciesAcceleration the
    state jumps at each acceleration step, so use a window of at least
    blockTime. Averages restart when the mesh changes.

\*---------------------------------------------------------------------------*/

#ifndef convergenceStop_H
#define convergenceStop_H

#include "fvCFD.H"

namespace Foam
{

class convergenceStop
{
    const fvMesh& mesh_;
    bool active_;
    wordList fields_;
    scalar window_;
    scalar tolerance_;
    label nWindows_;
    scalar startTime_;

    // Time integrals over the current window and averages of the last one
    List<scalarField> integral_;
    List<scalarField> previous_;
    scalar integratedTime_;
    bool hasPrevious_;
    label nBelow_;

    void reset()
    {
        forAll(fields_, i)
        {
            integral_[i].setSize(mesh_.nCells());
            integral_[i] = 0;
        }
        integratedTime_ = 0;
        hasPrevious_ = false;
        nBelow_ = 0;
    }

public:

    explicit convergenceStop(const fvMesh& mesh)
    :
        mesh_(mesh),
        active_(mesh.time().controlDict().found("convergenceStop")),
        window_(GREAT),
        tolerance_(0),
        nWindows_(3),
        startTime_(-GREAT),
        integratedTime_(0),
        hasPrevious_(false),
        nBelow_(0)
    {
        if (!active_)
        {
            return;
        }

        const dictionary& dict =
            mesh.time().controlDict().subDict("convergenceStop");

        fields_ = wordList(dict.lookup("fields"));
        window_ = readScalar(dict.lookup("window"));
        tolerance_ = readScalar(dict.lookup("tolerance"));
        nWindows_ = dict.lookupOrDefault<label>("nWindows", 3);
        startTime_ = dict.lookupOrDefault<scalar>("startTime", -GREAT);

        if (fields_.empty() || window_ <= 0 || tolerance_ <= 0 || nWindows_ < 1)
        {
            FatalIOErrorIn("convergenceStop::convergenceStop", dict)
                << "convergenceStop needs at least one field, a positive"
                << " window and tolerance, and nWindows of at least 1"
                << exit(FatalIOError);
        }

        forAll(fields_, i)
        {
            if (!mesh.foundObject<volScalarField>(fields_[i]))
            {
                FatalIOErrorIn("convergenceStop::convergenceStop", dict)
                    << "convergenceStop: no field " << fields_[i]
                    << " on mesh " << mesh.name()
                    << exit(FatalIOError);
            }
        }

        integral_.setSize(fields_.size());
        previous_.setSize(fields_.size());
        reset();

        Info<< "convergenceStop: " << fields_ << ", window " << window_
            << " s, tolerance " << tolerance_ << ", " << nWindows_
            << " windows" << endl;
    }

    //- Call once per accepted step; true when the run has converged
    bool converged()
    {
        if (!active_)
        {
            return false;
        }

        const Time& runTime = mesh_.time();
        const scalar dt = runTime.deltaT0Value();

        if (runTime.value() - dt < startTime_)
        {
            return false;
        }

        bool sizeChanged = integral_[0].size() != mesh_.nCells();
        reduce(sizeChanged, orOp<bool>());

        if (sizeChanged || mesh_.changing())
        {
            reset();
            Info<< "convergenceStop: mesh changed, averages restarted" << endl;
        }

        forAll(fields_, i)
        {
            integral_[i] +=
                dt*mesh_.lookupObject<volScalarField>(fields_[i]).internalField();
        }
        integratedTime_ += dt;

        if (integratedTime_ < window_*(1 - SMALL) - 0.5*dt)
        {
            return false;
        }

        const scalarField& V = mesh_.V();
        scalar largest = 0;

        forAll(fields_, i)
        {
            const scalarField average(integral_[i]/integratedTime_);

            if (hasPrevious_)
            {
                const scalar denominator =
                    sqrt(gSum(V*sqr(previous_[i])));
                const scalar change =
                    sqrt(gSum(V*sqr(average - previous_[i])))
                   /max(denominator, VSMALL);

                largest = max(largest, change);

                Info<< "convergenceStop: " << fields_[i] << " change "
                    << change << endl;
            }

            previous_[i] = average;
            integral_[i] = 0;
        }

        const bool compared = hasPrevious_;
        hasPrevious_ = true;
        integratedTime_ = 0;

        if (!compared)
        {
            return false;
        }

        nBelow_ = largest < tolerance_ ? nBelow_ + 1 : 0;

        Info<< "convergenceStop: largest change " << largest << " at "
            << runTime.timeName() << " s, " << nBelow_ << " of " << nWindows_
            << " windows below " << tolerance_ << endl;

        if (nBelow_ >= nWindows_)
        {
            Info<< "convergenceStop: converged at " << runTime.timeName()
                << " s" << endl;

            return true;
        }

        return false;
    }
};

} // End namespace Foam

#endif

// ************************************************************************* //
