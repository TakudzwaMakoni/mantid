//----------------------------------------------------------------------
// Includes
//----------------------------------------------------------------------
#include "MantidCurveFitting/IkedaCarpenterPV.h"
#include "MantidCurveFitting/BoundaryConstraint.h"
#include "MantidCurveFitting/SpecialFunctionSupport.h"
#include "MantidAPI/MatrixWorkspace.h"
#include "MantidKernel/UnitFactory.h"
#include <cmath>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_multifit_nlin.h>
#include <limits>
#include "MantidGeometry/Instrument/DetectorGroup.h"
#include "MantidGeometry/Instrument/ParameterMap.h"
#include "MantidGeometry/Instrument/Component.h"
#include "MantidGeometry/Instrument/DetectorGroup.h"
#include "MantidGeometry/Instrument/FitParameter.h"
#include <limits>

namespace Mantid
{
namespace CurveFitting
{

using namespace Kernel;
using namespace SpecialFunctionSupport;
using namespace Geometry;

DECLARE_FUNCTION(IkedaCarpenterPV)

// Get a reference to the logger
Kernel::Logger& IkedaCarpenterPV::g_log = Kernel::Logger::get("IkedaCarpenterPV");

double IkedaCarpenterPV::centre()const 
{
  return getParameter("X0");
}


void IkedaCarpenterPV::setHeight(const double h) 
{
  // calculate height of peakshape function corresponding to intensity = 1
  setParameter("I",1);
  double h0 = height();

  // to avoid devide by zero and to insane value for I to be set
  double minCutOff = 100.0*std::numeric_limits<double>::min();
  if ( h0 > 0 && h0 < minCutOff )
    h0 = minCutOff;
  if ( h0 < 0 && h0 > -minCutOff )
    h0 = -minCutOff;

  // The intensity is then estimated to be h/h0
  setParameter("I", h/h0);
};


double IkedaCarpenterPV::height()const 
{
  // return the function value at centre()
  double h0;
  double toCentre = centre();
  constFunction(&h0, &toCentre, 1);
  return h0;
};

double IkedaCarpenterPV::width()const 
{
  double sigmaSquared = getParameter("SigmaSquared");
  double gamma = getParameter("Gamma");

  if ( sigmaSquared < 0 )
  {
    g_log.warning() << "SigmaSquared NEGATIVE!.\n"
                    << "Likely due to a fit not converging properly\n"
                    << "If this is frequent problem please report to Mantid team.\n"
                    << "For now to calculate width force SigmaSquared positive.\n";
    sigmaSquared = - sigmaSquared;
  }
  if ( gamma < 0 )
  {
    g_log.warning() << "Gamma NEGATIVE!.\n"
                    << "Likely due to a fit not converging properly\n"
                    << "If this is frequent problem please report to Mantid team.\n"
                    << "For now to calculate width force Gamma positive.\n";
    gamma = - gamma;;
  }
  return sqrt(8.0*M_LN2*sigmaSquared)+gamma;
};

void IkedaCarpenterPV::setWidth(const double w) 
{
  setParameter("SigmaSquared", w*w/(32.0*M_LN2));  // used 4.0 * 8.0 = 32.0
  setParameter("Gamma", w/2.0);
};

void IkedaCarpenterPV::setCentre(const double c) 
{
  setParameter("X0",c);
};


void IkedaCarpenterPV::init()
{
  declareParameter("I", 0.0);
  declareParameter("Alpha0",1.6);
  declareParameter("Alpha1",1.5);
  declareParameter("Beta0",31.9);
  declareParameter("Kappa",46.0);
  declareParameter("SigmaSquared",1.0);
  declareParameter("Gamma",1.0);
  declareParameter("X0",0.0);
}


/** Method for updating m_waveLength.
 *  If size of m_waveLength is equal to number of data (for a new instance of this 
 *  class this vector is empty initially) then don't recalculate it.
 *
 *  @param xValues :: x values
 *  @param nData :: length of xValues
 */
void IkedaCarpenterPV::calWavelengthAtEachDataPoint(const double* xValues, const size_t& nData) const
{
    // if wavelength vector already have the right size no need for resizing it
    // further we make the assumption that no need to recalculate this vector if
    // it already has the right size

    if (m_waveLength.size() != nData)
    {
      m_waveLength.resize(nData);

      Mantid::Kernel::Unit_sptr wavelength = Mantid::Kernel::UnitFactory::Instance().create("Wavelength");
      for (size_t i = 0; i < nData; i++)
      {
        m_waveLength[i] = xValues[i]; 
      }

      // note if a version of convertValue was added which allows a double* as first argument
      // then could avoid copying above plus only have to resize m_wavelength when 
      // its size smaller than nData
      if ( m_workspace != 0 )
      {
        IInstrument_const_sptr instrument = m_workspace->getInstrument();
        Geometry::IObjComponent_const_sptr sample = instrument->getSample();
        if (sample != NULL)
        {
	        convertValue(m_waveLength, wavelength, m_workspace, m_workspaceIndex);
        }
        else
        {
          g_log.warning() << "No sample set for instrument in workspace.\n"
                          << "Can't calculate wavelength in IkedaCarpenter.\n"
                          << "Default all wavelengths to one.\n"
                          << "Solution is to load appropriate instrument into workspace.\n";
          for (size_t i = 0; i < nData; i++)
            m_waveLength[i] = 1.0; 
        }
      }
      else
      {
        g_log.warning() << "Workspace not set.\n"
                        << "Can't calculate wavelength in IkedaCarpenter.\n"
                        << "Default all wavelengths to one.\n"
                        << "Solution call setMatrixWorkspace() for function.\n";
        for (size_t i = 0; i < nData; i++)
          m_waveLength[i] = 1.0; 
      
      }
    }
}


/** convert voigt params to pseudo voigt params
 *
 *  @param voigtSigmaSq :: voigt param
 *  @param voigtGamma :: voigt param
 *  @param H :: pseudo voigt param
 *  @param eta :: pseudo voigt param
 */
void IkedaCarpenterPV::convertVoigtToPseudo(const double& voigtSigmaSq, const double& voigtGamma, 
  double& H, double& eta) const
{
  double fwhmGsq = 8.0 * M_LN2 * voigtSigmaSq;
  double fwhmG = sqrt(fwhmGsq);
  double fwhmG4 = fwhmGsq*fwhmGsq;
  double fwhmL = voigtGamma;
  double fwhmLsq = voigtGamma*voigtGamma;
  double fwhmL4 = fwhmLsq*fwhmLsq;

  H = pow(fwhmG4*fwhmG+2.69269*fwhmG4*fwhmL+2.42843*fwhmGsq*fwhmG*fwhmLsq
    +4.47163*fwhmGsq*fwhmLsq*fwhmL+0.07842*fwhmG*fwhmL4+fwhmL4*fwhmL, 0.2);

  if (H == 0.0)
    H = std::numeric_limits<double>::epsilon()*1000.0;

  double tmp = fwhmL/H;

  eta = 1.36603*tmp - 0.47719*tmp*tmp + 0.11116*tmp*tmp*tmp;
}

void IkedaCarpenterPV::constFunction(double* out, const double* xValues, const int& nData) const
{
    const double& I = getParameter("I");
    const double& alpha0 =getParameter("Alpha0");
    const double& alpha1 = getParameter("Alpha1");
    const double& beta0 = getParameter("Beta0");
    const double& kappa = getParameter("Kappa");
    const double& voigtsigmaSquared = getParameter("SigmaSquared");
    const double& voigtgamma = getParameter("Gamma");
    const double& X0 = getParameter("X0");

    // cal pseudo voigt sigmaSq and gamma and eta
    double gamma = 1.0; // dummy initialization
    double eta = 0.5;   // dummy initialization
    convertVoigtToPseudo(voigtsigmaSquared, voigtgamma, gamma, eta);
    double sigmaSquared = gamma*gamma/(8.0*M_LN2); 

    const double beta = 1/beta0;

    // equations taken from Fullprof manual

    const double k = 0.05;   

    double u,v,s,r;
    double yu, yv, ys, yr;

    // Not entirely sure what to do if sigmaSquared ever negative
    // for now just post a warning
    double someConst = std::numeric_limits<double>::max() / 100.0;
    if ( sigmaSquared > 0 )
      someConst = 1/sqrt(2.0*sigmaSquared);
    else if (sigmaSquared < 0 )
    {
      g_log.warning() << "sigmaSquared negative in functionLocal.\n";
    }

    double R, Nu, Nv, Ns, Nr, N;

    std::complex<double> zs, zu, zv, zr;

    double alpha, a_minus, a_plus, x, y, z;

    // update wavelength vector
    calWavelengthAtEachDataPoint(xValues, nData);

    for (int i = 0; i < nData; i++) {
        double diff=xValues[i]-X0;

        R = exp(-81.799/(m_waveLength[i]*m_waveLength[i]*kappa));
        alpha = 1.0 / (alpha0+m_waveLength[i]*alpha1);

        a_minus = alpha*(1-k);
        a_plus = alpha*(1+k);
        x=a_minus-beta;
        y=alpha-beta;
        z=a_plus-beta; 

        Nu=1-R*a_minus/x;
        Nv=1-R*a_plus/z;
        Ns=-2*(1-R*alpha/y);
        Nr=2*R*alpha*alpha*beta*k*k/(x*y*z);

        u=a_minus*(a_minus*sigmaSquared-2*diff)/2.0;
        v=a_plus*(a_plus*sigmaSquared-2*diff)/2.0;
        s=alpha*(alpha*sigmaSquared-2*diff)/2.0;
        r=beta*(beta*sigmaSquared-2*diff)/2.0;

        yu = (a_minus*sigmaSquared-diff)*someConst;
        yv = (a_plus*sigmaSquared-diff)*someConst;
        ys = (alpha*sigmaSquared-diff)*someConst;
        yr = (beta*sigmaSquared-diff)*someConst;

        zs = std::complex<double>(-alpha*diff, 0.5*alpha*gamma);
        zu = (1-k)*zs;
        zv = (1-k)*zs;
        zr = std::complex<double>(-beta*diff, 0.5*beta*gamma);

        N = 0.25*alpha*(1-k*k)/(k*k);

        out[i] = I*N*( (1-eta)*(Nu*exp(u+gsl_sf_log_erfc(yu))+Nv*exp(v+gsl_sf_log_erfc(yv)) + 
                        Ns*exp(s+gsl_sf_log_erfc(ys))+Nr*exp(r+gsl_sf_log_erfc(yr))) -
                 eta*2.0/M_PI*(Nu*exponentialIntegral(zu).imag()+Nv*exponentialIntegral(zv).imag()
                              +Ns*exponentialIntegral(zs).imag()+Nr*exponentialIntegral(zr).imag()) );
    }

}

void IkedaCarpenterPV::functionLocal(double* out, const double* xValues, const size_t nData)const
{
    const double& I = getParameter("I");
    const double& alpha0 =getParameter("Alpha0");
    const double& alpha1 = getParameter("Alpha1");
    const double& beta0 = getParameter("Beta0");
    const double& kappa = getParameter("Kappa");
    const double& voigtsigmaSquared = getParameter("SigmaSquared");
    const double& voigtgamma = getParameter("Gamma");
    const double& X0 = getParameter("X0");

    // cal pseudo voigt sigmaSq and gamma and eta
    double gamma = 1.0; // dummy initialization
    double eta = 0.5;   // dummy initialization
    convertVoigtToPseudo(voigtsigmaSquared, voigtgamma, gamma, eta);
    double sigmaSquared = gamma*gamma/(8.0*M_LN2); // pseudo voigt sigma^2

    const double beta = 1/beta0;

    // equations taken from Fullprof manual

    const double k = 0.05;   

    double u,v,s,r;
    double yu, yv, ys, yr;

    // Not entirely sure what to do if sigmaSquared ever negative
    // for now just post a warning
    double someConst = std::numeric_limits<double>::max() / 100.0;
    if ( sigmaSquared > 0 )
      someConst = 1/sqrt(2.0*sigmaSquared);
    else if (sigmaSquared < 0 )
    {
      g_log.warning() << "sigmaSquared negative in functionLocal.\n";
    }

    double R, Nu, Nv, Ns, Nr, N;

    std::complex<double> zs, zu, zv, zr;

    double alpha, a_minus, a_plus, x, y, z;

    // update wavelength vector
    calWavelengthAtEachDataPoint(xValues, nData);

    for (size_t i = 0; i < nData; i++) {
        double diff=xValues[i]-X0;


        R = exp(-81.799/(m_waveLength[i]*m_waveLength[i]*kappa));
        alpha = 1.0 / (alpha0+m_waveLength[i]*alpha1);


        a_minus = alpha*(1-k);
        a_plus = alpha*(1+k);
        x=a_minus-beta;
        y=alpha-beta;
        z=a_plus-beta; 

        Nu=1-R*a_minus/x;
        Nv=1-R*a_plus/z;
        Ns=-2*(1-R*alpha/y);
        Nr=2*R*alpha*alpha*beta*k*k/(x*y*z);

        u=a_minus*(a_minus*sigmaSquared-2*diff)/2.0;
        v=a_plus*(a_plus*sigmaSquared-2*diff)/2.0;
        s=alpha*(alpha*sigmaSquared-2*diff)/2.0;
        r=beta*(beta*sigmaSquared-2*diff)/2.0;

        yu = (a_minus*sigmaSquared-diff)*someConst;
        yv = (a_plus*sigmaSquared-diff)*someConst;
        ys = (alpha*sigmaSquared-diff)*someConst;
        yr = (beta*sigmaSquared-diff)*someConst;

        zs = std::complex<double>(-alpha*diff, 0.5*alpha*gamma);
        zu = (1-k)*zs;
        zv = (1-k)*zs;
        zr = std::complex<double>(-beta*diff, 0.5*beta*gamma);

        N = 0.25*alpha*(1-k*k)/(k*k);

        out[i] = I*N*( (1-eta)*(Nu*exp(u+gsl_sf_log_erfc(yu))+Nv*exp(v+gsl_sf_log_erfc(yv)) + 
                        Ns*exp(s+gsl_sf_log_erfc(ys))+Nr*exp(r+gsl_sf_log_erfc(yr))) -
                 eta*2.0/M_PI*(Nu*exponentialIntegral(zu).imag()+Nv*exponentialIntegral(zv).imag()
                              +Ns*exponentialIntegral(zs).imag()+Nr*exponentialIntegral(zr).imag()) );
    }
}

void IkedaCarpenterPV::functionDerivLocal(API::Jacobian* out, const double* xValues, const size_t nData)
{
  calNumericalDeriv(out, xValues, nData);
}


} // namespace CurveFitting
} // namespace Mantid
