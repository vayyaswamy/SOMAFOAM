/*---------------------------------------------------------------------------*\
Copyright (C) 2018 by the LUEUR authors

License
This project is licensed under The 3-Clause BSD License. For further information
look for license file include with distribution.

Class
    Foam::ionParticle

Description
    Test ion for the ionTracker utility: a particle of the basic Lagrangian
    library with a velocity. Reaching a wall or a general patch ends the
    track and records the patch and the face; symmetry planes reflect.

\*---------------------------------------------------------------------------*/

#ifndef ionParticle_H
#define ionParticle_H

#include "Particle.H"
#include "CloudTemplate.H"
#include "wallPolyPatch.H"
#include "symmetryPolyPatch.H"
#include "wedgePolyPatch.H"
#include "cyclicPolyPatch.H"
#include "processorPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class ionParticle Declaration
\*---------------------------------------------------------------------------*/

class ionParticle
:
    public Particle<ionParticle>
{
    // Private data

        //- Velocity [m/s]
        vector U_;


public:

    friend class Cloud<ionParticle>;

    //- Tracking data: where the ion left the domain
    class trackData
    :
        public Particle<ionParticle>::trackData
    {
    public:

        //- Patch and patch face that the ion reached, -1 if none
        label hitPatch;
        label hitFace;

        trackData(Cloud<ionParticle>& cloud)
        :
            Particle<ionParticle>::trackData(cloud),
            hitPatch(-1),
            hitFace(-1)
        {}
    };


    // Constructors

        //- Construct from components
        ionParticle
        (
            const Cloud<ionParticle>& c,
            const vector& position,
            const label celli,
            const vector& U
        )
        :
            Particle<ionParticle>(c, position, celli),
            U_(U)
        {}

        //- Construct from Istream
        ionParticle
        (
            const Cloud<ionParticle>& c,
            Istream& is,
            bool readFields = true
        )
        :
            Particle<ionParticle>(c, is, readFields),
            U_(vector::zero)
        {}

        //- Construct and return a clone
        autoPtr<ionParticle> clone() const
        {
            return autoPtr<ionParticle>(new ionParticle(*this));
        }


    // Member functions

        const vector& U() const
        {
            return U_;
        }

        vector& U()
        {
            return U_;
        }


    // Patch interactions

        bool hitPatch(const polyPatch&, trackData&, const label)
        {
            return false;
        }

        void hitWedgePatch(const wedgePolyPatch&, trackData&)
        {}

        void hitSymmetryPatch(const symmetryPolyPatch& pp, trackData&)
        {
            vector n = pp.faceAreas()[pp.whichFace(face())];
            n /= mag(n);

            U_ -= 2*(U_ & n)*n;
        }

        void hitCyclicPatch(const cyclicPolyPatch&, trackData&)
        {}

        void hitProcessorPatch(const processorPolyPatch&, trackData& td)
        {
            td.switchProcessor = true;
        }

        void hitWallPatch(const wallPolyPatch& pp, trackData& td)
        {
            td.hitPatch = pp.index();
            td.hitFace = pp.whichFace(face());
            td.keepParticle = false;
        }

        void hitPatch(const polyPatch& pp, trackData& td)
        {
            td.hitPatch = pp.index();
            td.hitFace = pp.whichFace(face());
            td.keepParticle = false;
        }

        // The same for tracking without tracking data, which the library
        // instantiates as well; not used

        bool hitPatch(const polyPatch&, int&, const label)
        {
            return false;
        }

        void hitWedgePatch(const wedgePolyPatch&, int&)
        {}

        void hitSymmetryPatch(const symmetryPolyPatch&, int&)
        {}

        void hitCyclicPatch(const cyclicPolyPatch&, int&)
        {}

        void hitProcessorPatch(const processorPolyPatch&, int&)
        {}

        void hitWallPatch(const wallPolyPatch&, int&)
        {}

        void hitPatch(const polyPatch&, int&)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
