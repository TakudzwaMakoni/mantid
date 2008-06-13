//----------------------------------------------------------------------
// Includes
//----------------------------------------------------------------------
#include "MantidDataHandling/LoadInstrumentFromRaw.h"
#include "MantidAPI/Instrument.h"

#include "MantidKernel/ConfigService.h"
#include "MantidGeometry/Detector.h"
#include "MantidGeometry/CompAssembly.h"
#include "MantidGeometry/Component.h"
#include "LoadRaw/isisraw.h"

#include <fstream>


namespace Mantid
{
namespace DataHandling
{

DECLARE_ALGORITHM(LoadInstrumentFromRaw)

using namespace Kernel;
using namespace API;

Logger& LoadInstrumentFromRaw::g_log = Logger::get("LoadInstrumentFromRaw");

/// Empty default constructor
LoadInstrumentFromRaw::LoadInstrumentFromRaw()
{}

/// Initialisation method.
void LoadInstrumentFromRaw::init()
{
  // When used as a sub-algorithm the workspace name is not used - hence the "Anonymous" to satisfy the validator
  declareProperty(new WorkspaceProperty<Workspace>("Workspace","Anonymous",Direction::InOut));
  declareProperty("Filename","",new MandatoryValidator);
}

/** Executes the algorithm. Reading in the file and creating and populating
 *  the output workspace
 * 
 *  @throw FileError Thrown if unable to parse XML file
 */
void LoadInstrumentFromRaw::exec()
{
  // Retrieve the filename from the properties
  m_filename = getPropertyValue("Filename");

  // Get the input workspace
  const Workspace_sptr localWorkspace = getProperty("Workspace");

  // open raw file
  ISISRAW iraw(NULL);
  if (iraw.readFromFile(m_filename.c_str()) != 0)
  {
    g_log.error("Unable to open file " + m_filename);
    throw Exception::FileError("Unable to open File:" , m_filename);	  
  }

  // Get reference to Instrument and set its name
  boost::shared_ptr<API::Instrument> instrument = (localWorkspace->getInstrument());
  instrument->setName(iraw.i_inst);

  // Add dummy source and samplepos to instrument
  // The L2 and 2-theta values from Raw file assumed to be relative to sample position

  Geometry::ObjComponent *samplepos = new Geometry::ObjComponent;
  samplepos->setParent(instrument.get());
  instrument->add(samplepos);
  samplepos->setName("Unknown");
  instrument->markAsSamplePos(samplepos);
  samplepos->setPos(0.0,0.0,0.0);

  Geometry::ObjComponent *source = new Geometry::ObjComponent;
  source->setParent(instrument.get());
  instrument->add(source);
  source->setName("Unknown");
  instrument->markAsSource(source);
  double l1;
  if ( ! Kernel::ConfigService::Instance().getValue("instrument.L1", l1) )
  {
    l1 = 10.0;
  }
  source->setPos(0.0,-1.0*l1,0.0); 

  // add detectors

  int numDetector = iraw.i_det; // number of detectors
  int* detID = iraw.udet;      // detector IDs
  float* r = iraw.len2;         // distance from sample
  float* angle = iraw.tthe;     // angle between indicent beam and direction from sample to detector (two-theta)

  Geometry::Detector detector("det",samplepos);

  for (int i = 0; i < numDetector; i++)
  {

    Geometry::V3D pos;
    pos.spherical(r[i], angle[i], 0.0);
    detector.setPos(pos);

    // set detector ID, add copy to instrument and mark it

    detector.setID(detID[i]);
    int toGetHoldOfDetectorCopy = instrument->addCopy(&detector);
    Geometry::Detector* temp = dynamic_cast<Geometry::Detector*>((*instrument)[toGetHoldOfDetectorCopy-1]);
    instrument->markAsDetector(temp);
  }

  // Information to the user about what info is extracted from raw file

  g_log.information() << "SamplePos component added with position set to (0,0,0).\n"
    << "Detector components added with position coordinates assumed to be relative to the position of the sample; \n"
    << "L2 and two-theta values were read from raw file and used to set the r and theta spherical coordinates; \n"
    << "the remaining spherical coordinate phi was set to zero.\n"
    << "Source component added with position set to (0,-" << l1 << ",0). In standard configuration, with \n"
    << "the beam along y-axis pointing from source to sample, this implies the source is " << l1 << "m in front \n"
    << "of the sample. This value can be changed via the 'instrument.l1' configuration property.\n";

  return;
}


} // namespace DataHandling
} // namespace Mantid
