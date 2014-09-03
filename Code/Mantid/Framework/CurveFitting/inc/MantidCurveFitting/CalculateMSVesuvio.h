#ifndef MANTID_CURVEFITTING_CALCULATEMSVESUVIO_H_
#define MANTID_CURVEFITTING_CALCULATEMSVESUVIO_H_
//-----------------------------------------------------------------------------
// Includes
//-----------------------------------------------------------------------------
#include "MantidAPI/Algorithm.h"
#include "MantidGeometry/Instrument/ReferenceFrame.h"
#include "MantidKernel/V3D.h"

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real.hpp>

namespace Mantid
{
  namespace CurveFitting
  {
    //-----------------------------------------------------------------------------
    // CurveFitting forward declarations
    //-----------------------------------------------------------------------------
    struct DetectorParams;
    struct ResolutionParams;

    /**
      Calculates the multiple scattering & total scattering contributions
      for a flat-plate or cylindrical sample.

      Copyright &copy; 2014 ISIS Rutherford Appleton Laboratory &
      NScD Oak Ridge National Laboratory

      This file is part of Mantid.

      Mantid is free software; you can redistribute it and/or modify
      it under the terms of the GNU General Public License as published by
      the Free Software Foundation; either version 3 of the License, or
      (at your option) any later version.

      Mantid is distributed in the hope that it will be useful,
      but WITHOUT ANY WARRANTY; without even the implied warranty of
      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
      GNU General Public License for more details.

      You should have received a copy of the GNU General Public License
      along with this program.  If not, see <http://www.gnu.org/licenses/>.

      File change history is stored at: <https://github.com/mantidproject/mantid>
      Code Documentation is available at: <http://doxygen.mantidproject.org>
    */
    class DLLExport CalculateMSVesuvio : public API::Algorithm
    {
    private:
      // Produces random numbers with various probability distributions
      class RandomNumberGenerator
      {
        typedef boost::uniform_real<double> uniform_double;
        typedef boost::normal_distribution<double> gaussian_double;
      public:
        RandomNumberGenerator(const int seed);
        /// Returns a flat random number between 0.0 & 1.0
        double flat();
        /// Returns a random number distributed  by a normal distribution
        double gaussian(const double mean, const double sigma);

      private:
        boost::mt19937 m_generator;
      };
      // Stores counts for each scatter order
      // for a "run" of a given number of events
      struct Simulation
      {
        Simulation(const size_t order, const size_t ntimes);

        std::vector<std::vector<double>> counts;
        size_t maxorder;
        double weight;
        size_t nmscat;
      };
      // Accumulates and averages the results of each simulation
      struct SimulationAggregator
      {
        SimulationAggregator(const size_t nruns,
                             const size_t order,
                             const size_t ntimes);
        // Adds a result as part of the average
        void add(const Simulation & result);

        Simulation averaged;
        double prefactor;
      };

    public:
      CalculateMSVesuvio();
      ~CalculateMSVesuvio();

      virtual const std::string name() const;
      virtual int version() const;
      virtual const std::string category() const;
      virtual const std::string summary() const;

    private:
      void init();
      void exec();

      void cacheInputGeometry();
      void calculateMS(const size_t wsIndex, API::ISpectrum & totalsc,
                       API::ISpectrum & multsc) const;
      Simulation simulate(const size_t nevents, const size_t nscatters,
                            const DetectorParams & detpar,
                            const ResolutionParams &respar) const;
      double calculateCounts(const size_t nscatters,
                             const DetectorParams & detpar,
                             const ResolutionParams &respar,
                             Simulation & counts) const;

      // single-event helpers
      double generateSrcPos(const double l, Kernel::V3D & srcPos) const;
      double generateE0(const double l1, const double t2, double &weight) const;
      double generateTOF(const double gaussTOF, const double en0, const double dl1) const;
      Kernel::V3D generateScatter(const Kernel::V3D &startPos, const Kernel::V3D &direc,
                                  double &weight) const;

      // Member Variables
      RandomNumberGenerator *m_randgen; // random number generator

      size_t m_acrossDir, m_upDir, m_beamDir; // indices of each direction
      double m_srcR1; // beam umbra radius (m)
      double m_srcR2; // beam penumbra radius (m)
      double m_halfSampleHeight; // half-height of sample (m)
      double m_halfSampleWidth; // half-width of sample (m)
      double m_halfSampleThick; // half-thickness of sample(m)
      double m_maxWidthSampleFrame; // Maximum width in sample frame (m)
      Kernel::DblMatrix const *m_goniometer; // sample rotation
      Geometry::Object const *m_sampleShape; // sample shape
      Kernel::V3D m_samplePos;

      double m_mu; //attenuation xsec

      double m_tmin; // minimum tof value
      double m_tmax; // maximum tof value
      double m_dt; // tof value step

      API::Progress *m_progress;
      API::MatrixWorkspace_sptr m_inputWS;
    };

  } // namespace CurveFitting
} // namespace Mantid

#endif  /* MANTID_CURVEFITTING_CALCULATEMSVESUVIO_H_ */
