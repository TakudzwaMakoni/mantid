//----------------------------------------------------------------------
// Includes
//----------------------------------------------------------------------
#include "MantidGeometry/Crystal/IPeak.h"

#include "MantidAlgorithms/CheckWorkspacesMatch.h"
#include "MantidAPI/IMDWorkspace.h"
#include "MantidAPI/IMDEventWorkspace.h"
#include "MantidAPI/IMDHistoWorkspace.h"
#include "MantidAPI/IPeaksWorkspace.h"
#include "MantidAPI/WorkspaceGroup.h"
#include "MantidAPI/NumericAxis.h"
#include "MantidAPI/TableRow.h"
#include "MantidDataObjects/EventWorkspace.h"
#include "MantidDataObjects/Events.h"
#include "MantidGeometry/MDGeometry/IMDDimension.h"
#include <sstream>

//
namespace Mantid {
namespace Algorithms {

// Register the algorithm into the AlgorithmFactory
DECLARE_ALGORITHM(CheckWorkspacesMatch)

/// Constructor
CheckWorkspacesMatch::CheckWorkspacesMatch() : API::Algorithm() {}

/// Virtual destructor
CheckWorkspacesMatch::~CheckWorkspacesMatch() {}

using namespace Kernel;
using namespace API;
using namespace DataObjects;
using namespace Geometry;

//----------------------------------------------------------------------------------------------
void CheckWorkspacesMatch::init() {
  declareProperty(
      new WorkspaceProperty<Workspace>("Workspace1", "", Direction::Input),
      "The name of the first input workspace.");
  declareProperty(
      new WorkspaceProperty<Workspace>("Workspace2", "", Direction::Input),
      "The name of the second input workspace.");

  declareProperty(
      "Tolerance", 0.0,
      "The maximum amount by which values may differ between the workspaces.");

  declareProperty("CheckType", true, "Whether to check that the data types "
                                     "(Workspace2D vs EventWorkspace) match.");
  declareProperty("CheckAxes", true, "Whether to check that the axes match.");
  declareProperty("CheckSpectraMap", true,
                  "Whether to check that the spectra-detector maps match. ");
  declareProperty("CheckInstrument", true,
                  "Whether to check that the instruments match. ");
  declareProperty("CheckMasking", true,
                  "Whether to check that the bin masking matches. ");
  declareProperty(
      "CheckSample", false,
      "Whether to check that the sample (e.g. logs)."); // Have this one false
                                                        // by default - the logs
                                                        // are brittle

  declareProperty("Result", "", Direction::Output);

  declareProperty(
      "ToleranceRelErr", false,
      "Treat tolerance as relative error rather then the absolute error.\n"
      "This is only applicable to Matrix workspaces.");
  declareProperty("CheckAllData", false,
                  "Usually checking data ends when first mismatch occurs. This "
                  "forces algorithm to check all data and print mismatch to "
                  "the debug log.\n"
                  "Very often such logs are huge so making it true should be "
                  "the last option.");
  // Have this one false by default - it can be a lot of printing.

  declareProperty("NumberMismatchedSpectraToPrint", 1,
                  "Number of mismatched spectra from lowest to be listed. ");

  declareProperty("DetailedPrintIndex", EMPTY_INT(),
                  "Mismatched spectra that will be printed out in details. ");
}

//----------------------------------------------------------------------------------------------
void CheckWorkspacesMatch::exec() {
  // Run new algorithm
  auto result = runCompareWorkspaces();

  // Output as per previous behaviour
  if (result != successString()) {
    g_log.notice() << "The workspaces did not match: " << result << std::endl;
  }

  setProperty("Result", result);
}

//----------------------------------------------------------------------------------------------
/**
 * Process two groups and ensure the Result string is set properly on the final
 * algorithm
 *
 * returns True if everything executed correctly
 */
bool CheckWorkspacesMatch::processGroups() {
  // Run new algorithm
  auto result = runCompareWorkspaces();

  // Output as per previous behaviour
  if (result != successString()) {
    g_log.notice() << result << "\n";
  }

  setProperty("Result", result);

  setExecuted(true);
  notificationCenter().postNotification(
      new FinishedNotification(this, this->isExecuted()));

  return true;
}

//----------------------------------------------------------------------------------------------
/**
 * Run new CompareWorkspaces algorithm as a child algorithm.
 *
 * Result string formatted the same way as before; "Success!" when workspaces
 * match or a newline separated list of mismatch messages.
 *
 * @return A string containing either successString() or mismatch messages
 */
std::string CheckWorkspacesMatch::runCompareWorkspaces() {
  // This algorithm produces a single result string
  std::string result;

  // Use new CompareWorkspaces algorithm to perform comparison
  Algorithm_sptr compare = this->createChildAlgorithm("CompareWorkspaces");
  compare->setRethrows(true);
  compare->setLogging(false);

  // Forward workspace properties
  Workspace_sptr ws1 = getProperty("Workspace1");
  Workspace_sptr ws2 = getProperty("Workspace2");
  compare->setProperty("Workspace1", ws1);
  compare->setProperty("Workspace2", ws2);

  // Copy any other non-default properties
  const std::vector<Property *> &allProps = this->getProperties();
  auto propCount = allProps.size();
  for (size_t i = 0; i < propCount; ++i) {
    Property *prop = allProps[i];
    const std::string &pname = prop->name();

    if (!prop->isDefault() && pname != "Workspace1" && pname != "Workspace2" &&
        pname != "Result")
      compare->setPropertyValue(pname, prop->value());
  }

  // Execute comparison
  compare->execute();

  // Generate result string
  if (!compare->getProperty("Result")) {
    ITableWorkspace_sptr table = compare->getProperty("Messages");
    auto rowcount = table->rowCount();
    for (size_t i = 0; i < rowcount; ++i) {
      result += table->cell<std::string>(i, 0);
      if (i < (rowcount - 1))
        result += "\n";
    }
  } else {
    result = successString();
  }

  return result;
}

} // namespace Algorithms
} // namespace Mantid
