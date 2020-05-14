// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#pragma once

#include "MantidAPI/AnalysisDataService.h"
#include "MantidAPI/ISpectrum.h"
#include "MantidAPI/MatrixWorkspace.h"
#include "MantidAPI/NumericAxis.h"
#include "MantidAPI/Run.h"
#include "MantidAPI/Sample.h"
#include "MantidAPI/SpectraAxis.h"
#include "MantidAPI/SpectrumDetectorMapping.h"
#include "MantidAPI/SpectrumInfo.h"
#include "MantidAPI/WorkspaceFactory.h"
#include "MantidGeometry/Crystal/OrientedLattice.h"
#include "MantidGeometry/Instrument.h"
#include "MantidGeometry/Instrument/ComponentInfo.h"
#include "MantidGeometry/Instrument/Detector.h"
#include "MantidGeometry/Instrument/DetectorInfo.h"
#include "MantidGeometry/Instrument/ReferenceFrame.h"
#include "MantidGeometry/MDGeometry/MDImplicitFunction.h"
#include "MantidHistogramData/Histogram.h"
#include "MantidHistogramData/LinearGenerator.h"
#include "MantidHistogramData/LogarithmicGenerator.h"
#include "MantidIndexing/IndexInfo.h"
#include "MantidKernel/TimeSeriesProperty.h"
#include "MantidKernel/VMD.h"
#include "MantidKernel/make_cow.h"
#include "MantidTestHelpers/ComponentCreationHelper.h"
#include "MantidTestHelpers/FakeObjects.h"
#include "MantidTestHelpers/InstrumentCreationHelper.h"
#include "MantidTestHelpers/NexusTestHelper.h"
#include "MantidTestHelpers/ParallelRunner.h"
#include "MantidTypes/SpectrumDefinition.h"
#include "PropertyManagerHelper.h"

#include <cxxtest/TestSuite.h>

#include <memory>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <functional>
#include <limits>
#include <numeric>

using std::size_t;
using namespace Mantid;
using namespace Mantid::Kernel;
using namespace Mantid::API;
using namespace Mantid::Geometry;
using Mantid::Indexing::IndexInfo;
using Mantid::Types::Core::DateAndTime;

// Declare into the factory.
DECLARE_WORKSPACE(WorkspaceTester)

namespace {
/** Create a workspace with numSpectra, with
 * each spectrum having one detector, at id = workspace index.
 * @param numSpectra
 * @return
 */
std::shared_ptr<MatrixWorkspace> makeWorkspaceWithDetectors(size_t numSpectra,
                                                            size_t numBins) {
  std::shared_ptr<MatrixWorkspace> ws2 = std::make_shared<WorkspaceTester>();
  ws2->initialize(numSpectra, numBins, numBins);

  auto inst = std::make_shared<Instrument>("TestInstrument");
  // We get a 1:1 map by default so the detector ID should match the spectrum
  // number
  for (size_t i = 0; i < ws2->getNumberHistograms(); ++i) {
    // Create a detector for each spectra
    Detector *det = new Detector("pixel", static_cast<detid_t>(i), inst.get());
    det->setShape(
        ComponentCreationHelper::createSphere(0.01, V3D(0, 0, 0), "1"));
    inst->add(det);
    inst->markAsDetector(det);
    ws2->getSpectrum(i).addDetectorID(static_cast<detid_t>(i));
  }
  ws2->setInstrument(inst);
  return ws2;
}

void run_legacy_setting_spectrum_numbers_with_MPI(
    const Parallel::Communicator &comm) {
  using namespace Parallel;
  for (const auto storageMode : {StorageMode::MasterOnly, StorageMode::Cloned,
                                 StorageMode::Distributed}) {
    WorkspaceTester ws;
    if (comm.rank() == 0 || storageMode != StorageMode::MasterOnly) {
      Indexing::IndexInfo indexInfo(1000, storageMode, comm);
      ws.initialize(indexInfo,
                    HistogramData::Histogram(HistogramData::Points(1)));
    }
    if (storageMode == StorageMode::Distributed && comm.size() > 1) {
      TS_ASSERT_THROWS_EQUALS(ws.getSpectrum(0).setSpectrumNo(42),
                              const std::logic_error &e, std::string(e.what()),
                              "Setting spectrum numbers in MatrixWorkspace via "
                              "ISpectrum::setSpectrumNo is not possible in MPI "
                              "runs for distributed workspaces. Use "
                              "IndexInfo.");
    } else {
      if (comm.rank() == 0 || storageMode != StorageMode::MasterOnly) {
        TS_ASSERT_THROWS_NOTHING(ws.getSpectrum(0).setSpectrumNo(42));
      }
    }
  }
}
} // namespace

class MatrixWorkspaceTest : public CxxTest::TestSuite {
public:
  // This pair of boilerplate methods prevent the suite being created statically
  // This means the constructor isn't called when running other tests
  static MatrixWorkspaceTest *createSuite() {
    return new MatrixWorkspaceTest();
  }
  static void destroySuite(MatrixWorkspaceTest *suite) { delete suite; }

  MatrixWorkspaceTest() : ws(std::make_shared<WorkspaceTester>()) {
    ws->initialize(1, 1, 1);
  }

  void test_indexInfo() {
    WorkspaceTester ws;
    ws.initialize(3, 1, 1);
    const auto &indexInfo = ws.indexInfo();
    // IndexInfo contains spectrum numbers generated by WorkspaceTester.
    TS_ASSERT_EQUALS(indexInfo.spectrumNumber(0), 1);
    TS_ASSERT_EQUALS(indexInfo.spectrumNumber(1), 2);
    TS_ASSERT_EQUALS(indexInfo.spectrumNumber(2), 3);
    // Workspace tester contains detector IDs...
    TS_ASSERT_EQUALS(ws.getSpectrum(0).getDetectorIDs().size(), 1);
    TS_ASSERT_EQUALS(ws.getSpectrum(1).getDetectorIDs().size(), 1);
    TS_ASSERT_EQUALS(ws.getSpectrum(2).getDetectorIDs().size(), 1);
    // ... but spectrum definitions in IndexInfo are empty...
    TS_ASSERT_EQUALS((*indexInfo.spectrumDefinitions())[0],
                     SpectrumDefinition{});
    TS_ASSERT_EQUALS((*indexInfo.spectrumDefinitions())[1],
                     SpectrumDefinition{});
    TS_ASSERT_EQUALS((*indexInfo.spectrumDefinitions())[2],
                     SpectrumDefinition{});
    // ... since there is no instrument, i.e., all detector IDs are invalid.
    TS_ASSERT_EQUALS(ws.detectorInfo().size(), 0);
  }

  void test_setIndexInfo_size_failure() {
    WorkspaceTester ws;
    ws.initialize(3, 1, 1);
    IndexInfo bad(2);
    TS_ASSERT_THROWS_EQUALS(ws.setIndexInfo(std::move(bad)),
                            const std::invalid_argument &e,
                            std::string(e.what()),
                            "MatrixWorkspace::setIndexInfo: IndexInfo size "
                            "does not match number of histograms in workspace");
  }

  void test_setIndexInfo_bad_detector_index() {
    WorkspaceTester ws;
    ws.initialize(1, 1, 1);
    IndexInfo indices(1);
    indices.setSpectrumNumbers({2});
    std::vector<SpectrumDefinition> specDefs(1);
    specDefs[0].add(0);
    indices.setSpectrumDefinitions(specDefs);
    TS_ASSERT_THROWS_EQUALS(ws.setIndexInfo(std::move(indices)),
                            const std::invalid_argument &e,
                            std::string(e.what()),
                            "MatrixWorkspace: SpectrumDefinition contains an "
                            "out-of-range detector index, i.e., the spectrum "
                            "definition does not match the instrument in the "
                            "workspace.");
  }

  void test_setIndexInfo_bad_detector_time_index() {
    WorkspaceTester ws;
    ws.initialize(1, 1, 1);
    ws.setInstrument(
        ComponentCreationHelper::createTestInstrumentRectangular(1, 1));
    IndexInfo indices(1);
    indices.setSpectrumNumbers({2});
    std::vector<SpectrumDefinition> specDefs(1);
    specDefs[0].add(0, 1);
    indices.setSpectrumDefinitions(specDefs);
    TS_ASSERT_THROWS_EQUALS(ws.setIndexInfo(std::move(indices)),
                            const std::invalid_argument &e,
                            std::string(e.what()),
                            "MatrixWorkspace: SpectrumDefinition contains an "
                            "out-of-range time index for a detector, i.e., the "
                            "spectrum definition does not match the instrument "
                            "in the workspace.");
  }

  void test_setIndexInfo_default_mapping_failure() {
    WorkspaceTester ws;
    ws.initialize(3, 1, 1);
    // 2x2 = 4 pixels
    ws.setInstrument(
        ComponentCreationHelper::createTestInstrumentRectangular(1, 2));
    IndexInfo indices(3);
    TS_ASSERT_THROWS_EQUALS(ws.setIndexInfo(std::move(indices)),
                            const std::invalid_argument &e,
                            std::string(e.what()),
                            "MatrixWorkspace: IndexInfo does not contain "
                            "spectrum definitions so building a 1:1 mapping "
                            "from spectra to detectors was attempted, but the "
                            "number of spectra in the workspace is not equal "
                            "to the number of detectors in the instrument.");
  }

  void test_setIndexInfo_default_mapping() {
    WorkspaceTester ws;
    ws.initialize(4, 1, 1);
    // 2x2 = 4 pixels
    ws.setInstrument(
        ComponentCreationHelper::createTestInstrumentRectangular(1, 2));
    IndexInfo indices(4);

    TS_ASSERT_THROWS_NOTHING(ws.setIndexInfo(std::move(indices)));

    const auto &info = ws.indexInfo();
    TS_ASSERT(info.spectrumDefinitions());
    std::vector<SpectrumDefinition> specDefs(4);
    specDefs[0].add(0);
    specDefs[1].add(1);
    specDefs[2].add(2);
    specDefs[3].add(3);
    TS_ASSERT_EQUALS(*info.spectrumDefinitions(), specDefs);
    TS_ASSERT_EQUALS(ws.getSpectrum(0).getSpectrumNo(), 1);
    TS_ASSERT_EQUALS(ws.getSpectrum(1).getSpectrumNo(), 2);
    TS_ASSERT_EQUALS(ws.getSpectrum(2).getSpectrumNo(), 3);
    TS_ASSERT_EQUALS(ws.getSpectrum(3).getSpectrumNo(), 4);
    TS_ASSERT_EQUALS(ws.getSpectrum(0).getDetectorIDs(),
                     (std::set<detid_t>{4}));
    TS_ASSERT_EQUALS(ws.getSpectrum(1).getDetectorIDs(),
                     (std::set<detid_t>{5}));
    TS_ASSERT_EQUALS(ws.getSpectrum(2).getDetectorIDs(),
                     (std::set<detid_t>{6}));
    TS_ASSERT_EQUALS(ws.getSpectrum(3).getDetectorIDs(),
                     (std::set<detid_t>{7}));
  }

  void test_setIndexInfo_updates_ISpectrum() {
    // NOTE: This test checks if the IndexInfo set via
    // MatrixWorkspace::setIndexInfo() affects data stored in ISpectrum and
    // obtained via the legacy interface. THIS TEST SHOULD BE REMOVED ONCE THAT
    // INTERFACE IS BEING REMOVED.
    WorkspaceTester ws;
    ws.initialize(3, 1, 1);
    ws.setInstrument(
        ComponentCreationHelper::createTestInstrumentRectangular(1, 2));
    TS_ASSERT_EQUALS(ws.detectorInfo().size(), 4);
    IndexInfo indices(3);
    indices.setSpectrumNumbers({2, 4, 6});
    std::vector<SpectrumDefinition> specDefs(3);
    specDefs[0].add(0);
    specDefs[1].add(1);
    specDefs[2].add(2);
    specDefs[2].add(3);
    indices.setSpectrumDefinitions(specDefs);
    TS_ASSERT_THROWS_NOTHING(ws.setIndexInfo(std::move(indices)));
    TS_ASSERT_EQUALS(ws.getSpectrum(0).getSpectrumNo(), 2);
    TS_ASSERT_EQUALS(ws.getSpectrum(1).getSpectrumNo(), 4);
    TS_ASSERT_EQUALS(ws.getSpectrum(2).getSpectrumNo(), 6);
    TS_ASSERT_EQUALS(ws.getSpectrum(0).getDetectorIDs(),
                     (std::set<detid_t>{4}));
    TS_ASSERT_EQUALS(ws.getSpectrum(1).getDetectorIDs(),
                     (std::set<detid_t>{5}));
    TS_ASSERT_EQUALS(ws.getSpectrum(2).getDetectorIDs(),
                     (std::set<detid_t>{6, 7}));
  }

  void test_indexInfo_legacy_compatibility() {
    // NOTE: This test checks if the IndexInfo reference returned by
    // MatrixWorkspace::indexInfo() reflects changes done via the legacy
    // interface of ISpectrum. THIS TEST SHOULD BE REMOVED ONCE THAT INTERFACE
    // IS BEING REMOVED.
    WorkspaceTester ws;
    ws.initialize(1, 1, 1);
    ws.setInstrument(
        ComponentCreationHelper::createTestInstrumentRectangular(1, 2));
    const auto &indexInfo = ws.indexInfo();
    TS_ASSERT_EQUALS(indexInfo.spectrumNumber(0), 1);
    TS_ASSERT_EQUALS((*indexInfo.spectrumDefinitions())[0],
                     SpectrumDefinition{});
    ws.getSpectrum(0).setSpectrumNo(7);
    ws.getSpectrum(0).setDetectorID(7);
    // No changes -- old and new interface should not be mixed!
    TS_ASSERT_EQUALS(indexInfo.spectrumNumber(0), 1);
    TS_ASSERT_EQUALS((*indexInfo.spectrumDefinitions())[0],
                     SpectrumDefinition{});
    // After getting a new reference we should see the changes.
    const auto &indexInfo2 = ws.indexInfo();
    TS_ASSERT_EQUALS(indexInfo2.spectrumNumber(0), 7);
    SpectrumDefinition specDef;
    specDef.add(ws.detectorInfo().indexOf(7));
    TS_ASSERT_EQUALS((*indexInfo.spectrumDefinitions())[0], specDef);
  }

  void test_IndexInfo_copy() {
    WorkspaceTester ws;
    ws.initialize(3, 1, 1);
    ws.setInstrument(
        ComponentCreationHelper::createTestInstrumentRectangular(1, 2));
    IndexInfo indices(3);
    indices.setSpectrumNumbers({2, 4, 6});
    std::vector<SpectrumDefinition> specDefs(3);
    specDefs[0].add(0);
    specDefs[1].add(1);
    specDefs[2].add(2);
    specDefs[2].add(3);
    indices.setSpectrumDefinitions(specDefs);
    ws.setIndexInfo(std::move(indices));

    // Internally this references data in ISpectrum
    const auto &indexInfo = ws.indexInfo();
    // This should create a copy, dropping any links to MatrixWorkspace or
    // ISpectrum
    const auto copy(indexInfo);

    TS_ASSERT_EQUALS(copy.spectrumNumber(0), 2);
    TS_ASSERT_EQUALS(copy.spectrumNumber(1), 4);
    TS_ASSERT_EQUALS(copy.spectrumNumber(2), 6);
    TS_ASSERT_EQUALS(*copy.spectrumDefinitions(), specDefs);
    // Changing data in workspace affects indexInfo, but not copy
    ws.getSpectrum(0).setSpectrumNo(7);
    ws.getSpectrum(0).addDetectorID(7);
    const auto &indexInfo2 = ws.indexInfo();
    TS_ASSERT_EQUALS(indexInfo2.spectrumNumber(0), 7);
    SpectrumDefinition specDef;
    specDef.add(ws.detectorInfo().indexOf(4));
    specDef.add(ws.detectorInfo().indexOf(7));
    TS_ASSERT_EQUALS((*indexInfo2.spectrumDefinitions())[0], specDef);
    TS_ASSERT_EQUALS(copy.spectrumNumber(0), 2);
    TS_ASSERT_EQUALS((*copy.spectrumDefinitions())[0], specDefs[0]);
  }

  void test_setIndexInfo_shares_spectrumDefinition() {
    WorkspaceTester ws;
    ws.initialize(3, 1, 1);
    ws.setInstrument(
        ComponentCreationHelper::createTestInstrumentRectangular(1, 2));
    IndexInfo indices(3);
    indices.setSpectrumNumbers({2, 4, 6});

    auto defs = Kernel::make_cow<std::vector<SpectrumDefinition>>(3);
    TS_ASSERT_THROWS_NOTHING(indices.setSpectrumDefinitions(defs));
    TS_ASSERT_THROWS_NOTHING(ws.setIndexInfo(std::move(indices)));
    TS_ASSERT_EQUALS(ws.indexInfo().spectrumDefinitions().get(), defs.get());
  }

  void test_clone_shares_data_in_IndexInfo() {
    WorkspaceTester ws;
    ws.initialize(3, 1, 1);
    ws.setInstrument(
        ComponentCreationHelper::createTestInstrumentRectangular(1, 2));
    const auto clone = ws.clone();
    const auto &info1 = ws.indexInfo();
    const auto &info2 = clone->indexInfo();
    // Spectrum numbers should also be shared, but there is no access by
    // reference, so we cannot check.
    TS_ASSERT_EQUALS(info1.spectrumDefinitions(), info2.spectrumDefinitions());
    TS_ASSERT(info1.spectrumDefinitions()); // should not be nullptr
    TS_ASSERT_EQUALS(info1.spectrumDefinitions(), info2.spectrumDefinitions());
  }

  void test_WorkspaceFactory_shares_data_in_IndexInfo() {
    const auto ws =
        WorkspaceFactory::Instance().create("WorkspaceTester", 3, 1, 1);
    ws->setInstrument(
        ComponentCreationHelper::createTestInstrumentRectangular(1, 2));
    const auto copy = WorkspaceFactory::Instance().create(ws);
    const auto &info1 = ws->indexInfo();
    const auto &info2 = copy->indexInfo();
    // Spectrum numbers should also be shared, but there is no access by
    // reference, so we cannot check.
    TS_ASSERT_EQUALS(info1.spectrumDefinitions(), info2.spectrumDefinitions());
    TS_ASSERT(info1.spectrumDefinitions()); // should not be nullptr
    TS_ASSERT_EQUALS(info1.spectrumDefinitions(), info2.spectrumDefinitions());
  }

  void test_toString_Produces_Expected_Contents() {
    auto testWS = std::make_shared<WorkspaceTester>();
    testWS->initialize(1, 2, 1);
    testWS->setTitle("A test run");
    testWS->getAxis(0)->setUnit("TOF");
    testWS->setYUnitLabel("Counts");

    std::string expected = "WorkspaceTester\n"
                           "Title: A test run\n"
                           "Histograms: 1\n"
                           "Bins: 1\n"
                           "Histogram\n"
                           "X axis: Time-of-flight / microsecond\n"
                           "Y axis: Counts\n"
                           "Distribution: False\n"
                           "Instrument: None\n"
                           "Run start: not available\n"
                           "Run end:  not available\n";

    TS_ASSERT_EQUALS(expected, testWS->toString());
  }

  void test_initialize_with_IndexInfo_does_not_set_default_detectorIDs() {
    WorkspaceTester ws;
    Indexing::IndexInfo indexInfo(1);
    ws.initialize(indexInfo,
                  HistogramData::Histogram(HistogramData::Points(1)));
    TS_ASSERT_EQUALS(ws.getSpectrum(0).getDetectorIDs().size(), 0);
  }

  void testCloneClearsWorkspaceName() {
    auto ws = std::make_shared<WorkspaceTester>();
    ws->initialize(1, 1, 1);
    const std::string name{"MatrixWorkspace_testCloneClearsWorkspaceName"};
    AnalysisDataService::Instance().add(name, ws);
    TS_ASSERT_EQUALS(ws->getName(), name)
    auto cloned = ws->clone();
    TS_ASSERT(cloned->getName().empty())
    AnalysisDataService::Instance().clear();
  }

  void testGetSetTitle() {
    TS_ASSERT_EQUALS(ws->getTitle(), "");
    ws->setTitle("something");
    TS_ASSERT_EQUALS(ws->getTitle(), "something");
    ws->setTitle("");
  }

  void testGetSetComment() {
    TS_ASSERT_EQUALS(ws->getComment(), "");
    ws->setComment("commenting");
    TS_ASSERT_EQUALS(ws->getComment(), "commenting");
    ws->setComment("");
  }

  void test_getIndicesFromDetectorIDs() {
    WorkspaceTester ws;
    ws.initialize(10, 1, 1);
    for (size_t i = 0; i < 10; i++)
      ws.getSpectrum(i).setDetectorID(detid_t(i * 10));
    std::vector<detid_t> dets;
    dets.emplace_back(60);
    dets.emplace_back(20);
    dets.emplace_back(90);
    std::vector<size_t> indices = ws.getIndicesFromDetectorIDs(dets);
    TS_ASSERT_EQUALS(indices.size(), 3);
    TS_ASSERT_EQUALS(indices[0], 6);
    TS_ASSERT_EQUALS(indices[1], 2);
    TS_ASSERT_EQUALS(indices[2], 9);
  }

  void
  test_That_A_Workspace_Gets_SpectraMap_When_Initialized_With_NVector_Elements() {
    WorkspaceTester testWS;
    const size_t nhist(10);
    testWS.initialize(nhist, 1, 1);
    for (size_t i = 0; i < testWS.getNumberHistograms(); i++) {
      TS_ASSERT_EQUALS(testWS.getSpectrum(i).getSpectrumNo(), specnum_t(i + 1));
      TS_ASSERT(testWS.getSpectrum(i).hasDetectorID(detid_t(i)));
    }
  }

  void testEmptyWorkspace() {
    WorkspaceTester ws;
    TS_ASSERT(ws.isCommonBins());
    TS_ASSERT_EQUALS(ws.isCommonLogBins(), false);
    TS_ASSERT_EQUALS(ws.blocksize(), 0);
    TS_ASSERT_EQUALS(ws.size(), 0);
  }

  void testIsCommonBins() {
    WorkspaceTester ws;
    ws.initialize(10, 10, 10);
    // After initialization, ws should contain a shared HistogramX.
    TS_ASSERT(ws.isCommonBins());
    TS_ASSERT_EQUALS(ws.isCommonLogBins(), false);
    // Modifying the value of one Spectrum will cause this Histogram to detach.
    // Since the value is identical, isCommonBins is still true.
    ws.mutableX(0)[0] = 1.;
    TS_ASSERT(ws.isCommonBins());
    // Once we change the the value, however, isCommonsBins should return false
    ws.mutableX(0)[0] = 2.;
    TS_ASSERT_EQUALS(ws.isCommonBins(), false);
  }

  void testIsCommonLogAxis() {
    WorkspaceTester ws;
    ws.initialize(10, 10, 10);
    TS_ASSERT_EQUALS(ws.isCommonLogBins(), false);

    auto logAxis = Kernel::make_cow<Mantid::HistogramData::HistogramX>(
        10, Mantid::HistogramData::LogarithmicGenerator(1., 0.1));
    for (size_t i = 0; i < ws.getNumberHistograms(); ++i) {
      ws.setSharedX(i, logAxis);
    }
    TS_ASSERT_EQUALS(ws.isCommonLogBins(), true);

    auto linearAxis = Kernel::make_cow<Mantid::HistogramData::HistogramX>(
        10, Mantid::HistogramData::LinearGenerator(1., 0.1));
    for (size_t i = 0; i < ws.getNumberHistograms(); ++i) {
      ws.setSharedX(i, linearAxis);
    }
    TS_ASSERT_EQUALS(ws.isCommonLogBins(), false);
  }

  void test_updateSpectraUsing() {
    WorkspaceTester testWS;
    testWS.initialize(3, 1, 1);

    specnum_t specs[] = {1, 2, 2, 3};
    detid_t detids[] = {10, 99, 20, 30};
    TS_ASSERT_THROWS_NOTHING(
        testWS.updateSpectraUsing(SpectrumDetectorMapping(specs, detids, 4)));

    TS_ASSERT(testWS.getSpectrum(0).hasDetectorID(10));
    TS_ASSERT(testWS.getSpectrum(1).hasDetectorID(20));
    TS_ASSERT(testWS.getSpectrum(1).hasDetectorID(99));
    TS_ASSERT(testWS.getSpectrum(2).hasDetectorID(30));
  }

  void testDetectorMappingCopiedWhenAWorkspaceIsCopied() {
    std::shared_ptr<MatrixWorkspace> parent =
        std::make_shared<WorkspaceTester>();
    parent->initialize(1, 1, 1);
    parent->getSpectrum(0).setSpectrumNo(99);
    parent->getSpectrum(0).setDetectorID(999);

    MatrixWorkspace_sptr copied = WorkspaceFactory::Instance().create(parent);
    // Has it been copied?
    TS_ASSERT_EQUALS(copied->getSpectrum(0).getSpectrumNo(), 99);
    TS_ASSERT(copied->getSpectrum(0).hasDetectorID(999));
  }

  void testGetMemorySize() { TS_ASSERT_THROWS_NOTHING(ws->getMemorySize()); }

  void testHistory() { TS_ASSERT_THROWS_NOTHING(ws->history()); }

  void testAxes() { TS_ASSERT_EQUALS(ws->axes(), 2); }

  void testGetAxis() {
    TS_ASSERT_THROWS(ws->getAxis(-1), const Exception::IndexError &);
    TS_ASSERT_THROWS_NOTHING(ws->getAxis(0));
    TS_ASSERT(ws->getAxis(0));
    TS_ASSERT(ws->getAxis(0)->isNumeric());
    TS_ASSERT_THROWS(ws->getAxis(2), const Exception::IndexError &);
  }

  void testReplaceAxis() {
    auto ax = std::make_unique<SpectraAxis>(ws.get());
    auto ax1 = std::make_unique<SpectraAxis>(ws.get());
    TS_ASSERT_THROWS(ws->replaceAxis(2, std::move(ax)),
                     const Exception::IndexError &);
    TS_ASSERT_THROWS_NOTHING(ws->replaceAxis(0, std::move(ax1)));
    TS_ASSERT(ws->getAxis(0)->isSpectra());
  }

  void testIsDistribution() {
    TS_ASSERT(!ws->isDistribution());
    ws->setDistribution(true);
    TS_ASSERT(ws->isDistribution());
  }

  void testGetSetYUnit() {
    TS_ASSERT_EQUALS(ws->YUnit(), "");
    TS_ASSERT_THROWS_NOTHING(ws->setYUnit("something"));
    TS_ASSERT_EQUALS(ws->YUnit(), "something");
  }

  void testGetSpectrum() {
    WorkspaceTester ws;
    ws.initialize(4, 1, 1);
    TS_ASSERT_THROWS_NOTHING(ws.getSpectrum(0));
    TS_ASSERT_THROWS_NOTHING(ws.getSpectrum(3));
  }

  /** Get a detector sptr for each spectrum */
  void testGetDetector() {
    // Workspace has 3 spectra, each 1 in length
    const int numHist(3);
    std::shared_ptr<MatrixWorkspace> workspace(
        makeWorkspaceWithDetectors(3, 1));

    // Initially un masked
    for (int i = 0; i < numHist; ++i) {
      IDetector_const_sptr det;
      TS_ASSERT_THROWS_NOTHING(det = workspace->getDetector(i));
      if (det) {
        TS_ASSERT_EQUALS(det->getID(), i);
      } else {
        TS_FAIL("No detector defined");
      }
    }

    // Now a detector group
    auto &spec = workspace->getSpectrum(0);
    spec.addDetectorID(1);
    spec.addDetectorID(2);
    IDetector_const_sptr det;
    TS_ASSERT_THROWS_NOTHING(det = workspace->getDetector(0));
    TS_ASSERT(det);

    // Now an empty (no detector) pixel
    auto &spec2 = workspace->getSpectrum(1);
    spec2.clearDetectorIDs();
    IDetector_const_sptr det2;
    TS_ASSERT_THROWS_ANYTHING(det2 = workspace->getDetector(1));
    TS_ASSERT(!det2);
  }

  void testWholeSpectraMasking() {
    // Workspace has 3 spectra, each 1 in length
    const int numHist(3);
    std::shared_ptr<MatrixWorkspace> workspace(
        makeWorkspaceWithDetectors(3, 1));

    // Initially un masked
    const auto &spectrumInfo = workspace->spectrumInfo();
    for (int i = 0; i < numHist; ++i) {
      TS_ASSERT_EQUALS(workspace->readY(i)[0], 1.0);
      TS_ASSERT_EQUALS(workspace->readE(i)[0], 1.0);
      TS_ASSERT(spectrumInfo.hasDetectors(i));
      TS_ASSERT_EQUALS(spectrumInfo.isMasked(i), false);
    }

    // Mask a spectra
    workspace->getSpectrum(1).clearData();
    workspace->getSpectrum(2).clearData();
    workspace->mutableSpectrumInfo().setMasked(1, true);
    workspace->mutableSpectrumInfo().setMasked(2, true);

    const auto &spectrumInfo2 = workspace->spectrumInfo();
    for (int i = 0; i < numHist; ++i) {
      double expectedValue(0.0);
      bool expectedMasked(false);
      if (i == 0) {
        expectedValue = 1.0;
        expectedMasked = false;
      } else {
        expectedMasked = true;
      }
      TS_ASSERT_EQUALS(workspace->readY(i)[0], expectedValue);
      TS_ASSERT_EQUALS(workspace->readE(i)[0], expectedValue);
      TS_ASSERT(spectrumInfo2.hasDetectors(i));
      TS_ASSERT_EQUALS(spectrumInfo2.isMasked(i), expectedMasked);
    }
  }

  void testWholeSpectraMasking_SpectrumInfo() {
    // Workspace has 3 spectra, each 1 in length
    const int numHist(3);
    auto workspace = makeWorkspaceWithDetectors(numHist, 1);
    workspace->getSpectrum(1).clearData();
    workspace->getSpectrum(2).clearData();
    workspace->mutableSpectrumInfo().setMasked(1, true);
    workspace->mutableSpectrumInfo().setMasked(2, true);

    const auto &spectrumInfo = workspace->spectrumInfo();
    for (int i = 0; i < numHist; ++i) {
      bool expectedMasked(false);
      if (i == 0) {
        expectedMasked = false;
      } else {
        expectedMasked = true;
      }
      TS_ASSERT_EQUALS(spectrumInfo.isMasked(i), expectedMasked);
    }
  }

  void test_spectrumInfo_works_unthreaded() {
    const int numHist(3);
    auto workspace = makeWorkspaceWithDetectors(numHist, 1);
    std::atomic<bool> parallelException{false};
    for (int i = 0; i < numHist; ++i) {
      try {
        static_cast<void>(workspace->spectrumInfo());
      } catch (...) {
        parallelException = true;
      }
    }
    TS_ASSERT(!parallelException);
  }

  void test_spectrumInfo_works_threaded() {
    const int numHist(3);
    auto workspace = makeWorkspaceWithDetectors(numHist, 1);
    std::vector<const SpectrumInfo *> spectrumInfos(numHist);
    std::atomic<bool> parallelException{false};
    PARALLEL_FOR_IF(Kernel::threadSafe(*workspace))
    for (int i = 0; i < numHist; ++i) {
      try {
        spectrumInfos[i] = &(workspace->spectrumInfo());
      } catch (...) {
        parallelException = true;
      }
    }
    TS_ASSERT(!parallelException);
    for (int i = 0; i < numHist; ++i)
      TS_ASSERT_EQUALS(spectrumInfos[0], spectrumInfos[i]);
  }

  void testFlagMasked() {
    auto ws = makeWorkspaceWithDetectors(2, 2);
    // Now do a valid masking
    TS_ASSERT_THROWS_NOTHING(ws->flagMasked(0, 1, 0.75));
    TS_ASSERT(ws->hasAnyMaskedBins());
    TS_ASSERT(ws->hasMaskedBins(0));
    TS_ASSERT_EQUALS(ws->maskedBins(0).size(), 1);
    TS_ASSERT_EQUALS(ws->maskedBins(0).begin()->first, 1);
    TS_ASSERT_EQUALS(ws->maskedBins(0).begin()->second, 0.75);
    // flagMasked() shouldn't change the y-value maskBins() tested below does
    // that
    TS_ASSERT_EQUALS(ws->dataY(0)[1], 1.0);

    // Now mask a bin earlier than above and check it's sorting properly
    TS_ASSERT_THROWS_NOTHING(ws->flagMasked(1, 1))
    TS_ASSERT_EQUALS(ws->maskedBins(1).size(), 1)
    TS_ASSERT_EQUALS(ws->maskedBins(1).begin()->first, 1)
    TS_ASSERT_EQUALS(ws->maskedBins(1).begin()->second, 1.0)
    // Check the previous masking is still OK
    TS_ASSERT_EQUALS(ws->maskedBins(0).rbegin()->first, 1)
    TS_ASSERT_EQUALS(ws->maskedBins(0).rbegin()->second, 0.75)
  }

  void testMasking() {
    auto ws2 = makeWorkspaceWithDetectors(1, 2);

    TS_ASSERT(!ws2->hasAnyMaskedBins());
    TS_ASSERT(!ws2->hasMaskedBins(0));
    // Doesn't throw on invalid spectrum number, just returns false
    TS_ASSERT(!ws2->hasMaskedBins(1));
    TS_ASSERT(!ws2->hasMaskedBins(-1));

    // Will throw if nothing masked for spectrum
    TS_ASSERT_THROWS(ws2->maskedBins(0),
                     const Mantid::Kernel::Exception::IndexError &);
    // Will throw if attempting to mask invalid spectrum
    TS_ASSERT_THROWS(ws2->maskBin(-1, 1),
                     const Mantid::Kernel::Exception::IndexError &);
    TS_ASSERT_THROWS(ws2->maskBin(1, 1),
                     const Mantid::Kernel::Exception::IndexError &);
    // ...or an invalid bin
    TS_ASSERT_THROWS(ws2->maskBin(0, -1),
                     const Mantid::Kernel::Exception::IndexError &);
    TS_ASSERT_THROWS(ws2->maskBin(0, 2),
                     const Mantid::Kernel::Exception::IndexError &);

    // Now do a valid masking
    TS_ASSERT_THROWS_NOTHING(ws2->maskBin(0, 1, 0.5));
    TS_ASSERT(ws2->hasMaskedBins(0));
    TS_ASSERT_EQUALS(ws2->maskedBins(0).size(), 1);
    TS_ASSERT_EQUALS(ws2->maskedBins(0).begin()->first, 1);
    TS_ASSERT_EQUALS(ws2->maskedBins(0).begin()->second, 0.5);
    TS_ASSERT_EQUALS(ws2->dataY(0)[1], 0.5);

    // Now mask a bin earlier than above and check it's sorting properly
    TS_ASSERT_THROWS_NOTHING(ws2->maskBin(0, 0));
    TS_ASSERT_EQUALS(ws2->maskedBins(0).begin()->first, 0);
    TS_ASSERT_EQUALS(ws2->maskedBins(0).begin()->second, 1.0);
    TS_ASSERT_EQUALS(ws2->dataY(0)[0], 0.0);
    // Check the previous masking is still OK
    TS_ASSERT_EQUALS(ws2->maskedBins(0).rbegin()->first, 1);
    TS_ASSERT_EQUALS(ws2->maskedBins(0).rbegin()->second, 0.5);
    TS_ASSERT_EQUALS(ws2->dataY(0)[1], 0.5);
  }

  void testMaskingNaNInf() {
    const size_t s = 4;
    const double y[s] = {NAN, INFINITY, -INFINITY, 2.};
    WorkspaceTester ws;
    ws.initialize(1, s + 1, s);

    // initialize and mask first with 0 weights
    // masking with 0 weight should be equiavalent to flagMasked
    // i.e. values should not change, even Inf and NaN
    for (size_t i = 0; i < s; ++i) {
      ws.mutableY(0)[i] = y[i];
      ws.maskBin(0, i, 0);
    }

    TS_ASSERT(std::isnan(ws.y(0)[0]));
    TS_ASSERT(std::isinf(ws.y(0)[1]));
    TS_ASSERT(std::isinf(ws.y(0)[2]));
    TS_ASSERT_EQUALS(ws.y(0)[3], 2.);

    // now mask w/o specifying weight (e.g. 1 by default)
    // in this case everything should be 0, even NaN and Inf
    for (size_t i = 0; i < s; ++i) {
      ws.maskBin(0, i);
      TS_ASSERT_EQUALS(ws.y(0)[i], 0.);
    }
  }

  void testSetMaskedBins() {
    auto ws = makeWorkspaceWithDetectors(2, 2);
    ws->flagMasked(0, 1);
    ws->flagMasked(1, 0);
    ws->setMaskedBins(1, ws->maskedBins(0));
    TS_ASSERT(ws->hasMaskedBins(1));
    TS_ASSERT_EQUALS(ws->maskedBins(1).size(), 1);
    TS_ASSERT_EQUALS(ws->maskedBins(0).begin()->first, 1);
  }

  void testSize() {
    WorkspaceTester wkspace;
    wkspace.initialize(1, 4, 3);
    TS_ASSERT_EQUALS(wkspace.blocksize(), 3);
    TS_ASSERT_EQUALS(wkspace.size(), 3);
  }

  void test_hasGroupedDetectors() {
    auto ws = makeWorkspaceWithDetectors(5, 1);
    TS_ASSERT_EQUALS(ws->hasGroupedDetectors(), false);

    ws->getSpectrum(0).addDetectorID(3);
    TS_ASSERT_EQUALS(ws->hasGroupedDetectors(), true);
  }

  void test_getSpectrumToWorkspaceIndexMap() {
    WorkspaceTester ws;
    ws.initialize(2, 1, 1);
    auto map = ws.getSpectrumToWorkspaceIndexMap();
    TS_ASSERT_EQUALS(0, map[1]);
    TS_ASSERT_EQUALS(1, map[2]);
    TS_ASSERT_EQUALS(map.size(), 2);

    // Check it throws for non-spectra axis
    ws.replaceAxis(1, std::make_unique<NumericAxis>(1));
    TS_ASSERT_THROWS(ws.getSpectrumToWorkspaceIndexMap(),
                     const std::runtime_error &);
  }

  void test_getDetectorIDToWorkspaceIndexMap() {
    auto ws = makeWorkspaceWithDetectors(5, 1);
    detid2index_map idmap = ws->getDetectorIDToWorkspaceIndexMap(true);

    TS_ASSERT_EQUALS(idmap.size(), 5);
    int i = 0;
    for (auto it = idmap.begin(); it != idmap.end(); ++it, ++i) {
      TS_ASSERT_EQUALS(idmap.count(i), 1);
      TS_ASSERT_EQUALS(idmap[i], i);
    }

    ws->getSpectrum(2).addDetectorID(99); // Set a second ID on one spectrum
    TS_ASSERT_THROWS(ws->getDetectorIDToWorkspaceIndexMap(true),
                     const std::runtime_error &);
    detid2index_map idmap2 = ws->getDetectorIDToWorkspaceIndexMap();
    TS_ASSERT_EQUALS(idmap2.size(), 6);
  }

  void test_getDetectorIDToWorkspaceIndexVector() {
    auto ws = makeWorkspaceWithDetectors(100, 10);
    std::vector<size_t> out;
    detid_t offset = -1234;
    TS_ASSERT_THROWS_NOTHING(
        out = ws->getDetectorIDToWorkspaceIndexVector(offset));
    TS_ASSERT_EQUALS(offset, 0);
    TS_ASSERT_EQUALS(out.size(), 100);
    TS_ASSERT_EQUALS(out[0], 0);
    TS_ASSERT_EQUALS(out[1], 1);
    TS_ASSERT_EQUALS(out[99], 99);

    // Create some discontinuities and check that the default value is there
    // Have to create a whole new instrument to keep things consistent, since
    // the detector ID is stored in at least 3 places
    auto inst = std::make_shared<Instrument>("TestInstrument");
    // We get a 1:1 map by default so the detector ID should match the spectrum
    // number
    for (size_t i = 0; i < ws->getNumberHistograms(); ++i) {
      detid_t detid = static_cast<detid_t>(i);
      // Create a detector for each spectra
      if (i == 0)
        detid = -1;
      if (i == 99)
        detid = 110;
      Detector *det = new Detector("pixel", detid, inst.get());
      inst->add(det);
      inst->markAsDetector(det);
      ws->getSpectrum(i).addDetectorID(detid);
    }
    ws->setInstrument(inst);
    ws->getSpectrum(66).clearDetectorIDs();

    TS_ASSERT_THROWS_NOTHING(
        out = ws->getDetectorIDToWorkspaceIndexVector(offset));
    TS_ASSERT_EQUALS(offset, 1);
    TS_ASSERT_EQUALS(out.size(), 112);
    TS_ASSERT_EQUALS(out[66 + offset], std::numeric_limits<size_t>::max());
    TS_ASSERT_EQUALS(out[99 + offset], 99);
    TS_ASSERT_EQUALS(out[105 + offset], std::numeric_limits<size_t>::max());
    TS_ASSERT_EQUALS(out[110 + offset], 99);
  }

  void test_getSpectrumToWorkspaceIndexVector() {
    auto ws = makeWorkspaceWithDetectors(100, 10);
    std::vector<size_t> out;
    detid_t offset = -1234;
    TS_ASSERT_THROWS_NOTHING(out =
                                 ws->getSpectrumToWorkspaceIndexVector(offset));
    TS_ASSERT_EQUALS(offset, -1);
    TS_ASSERT_EQUALS(out.size(), 100);
    TS_ASSERT_EQUALS(out[0], 0);
    TS_ASSERT_EQUALS(out[1], 1);
    TS_ASSERT_EQUALS(out[99], 99);
  }

  void test_getSignalAtCoord_histoData() {
    // Create a test workspace
    const auto ws = createTestWorkspace(4, 6, 5);

    // Get signal at coordinates
    std::vector<coord_t> coords = {0.5, 1.0};
    TS_ASSERT_DELTA(
        ws.getSignalAtCoord(coords.data(), Mantid::API::NoNormalization), 0.0,
        1e-5);
    coords[0] = 1.5;
    TS_ASSERT_DELTA(
        ws.getSignalAtCoord(coords.data(), Mantid::API::NoNormalization), 1.0,
        1e-5);
  }

  void test_getSignalAtCoord_pointData() {
    // Create a test workspace
    const auto ws = createTestWorkspace(4, 5, 5);
    auto normType = Mantid::API::NoNormalization;

    // Get signal at coordinates
    std::vector<coord_t> coords = {-1.0, 1.0};
    coords[0] = -0.75;
    TS_ASSERT(std::isnan(ws.getSignalAtCoord(coords.data(), normType)));
    coords[0] = -0.25;
    TS_ASSERT_DELTA(ws.getSignalAtCoord(coords.data(), normType), 0.0, 1e-5);
    coords[0] = 0.0;
    TS_ASSERT_DELTA(ws.getSignalAtCoord(coords.data(), normType), 0.0, 1e-5);
    coords[0] = 0.25;
    TS_ASSERT_DELTA(ws.getSignalAtCoord(coords.data(), normType), 0.0, 1e-5);
    coords[0] = 0.75;
    TS_ASSERT_DELTA(ws.getSignalAtCoord(coords.data(), normType), 1.0, 1e-5);
    coords[0] = 1.0;
    TS_ASSERT_DELTA(ws.getSignalAtCoord(coords.data(), normType), 1.0, 1e-5);
    coords[0] = 4.25;
    TS_ASSERT_DELTA(ws.getSignalAtCoord(coords.data(), normType), 4.0, 1e-5);
    coords[0] = 4.75;
    TS_ASSERT(std::isnan(ws.getSignalAtCoord(coords.data(), normType)));
  }

  void test_getCoordAtSignal_regression() {
    /*
    Having more spectrum numbers (acutally vertical axis increments) than x bins
    in VolumeNormalisation mode
    should not cause any issues.
    */
    WorkspaceTester ws;
    const int nVertical = 4;

    const int nBins = 2;
    const int nYValues = 1;
    ws.initialize(nVertical, nBins, nYValues);
    auto verticalAxis = std::make_unique<NumericAxis>(nVertical);
    for (int i = 0; i < nVertical; ++i) {
      for (int j = 0; j < nBins; ++j) {
        if (j < nYValues) {
          ws.dataY(i)[j] = 1.0; // All y values are 1.
          ws.dataE(i)[j] = j;
        }
        ws.dataX(i)[j] = j; // x increments by 1
      }
      verticalAxis->setValue(i, double(i)); // Vertical axis increments by 1.
    }
    ws.replaceAxis(1, std::move(verticalAxis));
    // Signal is always 1 and volume of each box is 1. Therefore normalized
    // signal values by volume should always be 1.

    // Test at the top right.
    coord_t coord_top_right[2] = {static_cast<float>(ws.readX(0).back()),
                                  float(0)};
    signal_t value = 0;
    TS_ASSERT_THROWS_NOTHING(
        value = ws.getSignalAtCoord(coord_top_right, VolumeNormalization));
    TS_ASSERT_EQUALS(1.0, value);

    // Test at another location just to be sure.
    coord_t coord_bottom_left[2] = {
        static_cast<float>(ws.readX(nVertical - 1)[1]), float(nVertical - 1)};
    TS_ASSERT_THROWS_NOTHING(
        value = ws.getSignalAtCoord(coord_bottom_left, VolumeNormalization));
    TS_ASSERT_EQUALS(1.0, value);
  }

  void test_setMDMasking() {
    WorkspaceTester ws;
    TSM_ASSERT_THROWS("Characterisation test. This is not implemented.",
                      ws.setMDMasking(nullptr), const std::runtime_error &);
  }

  void test_clearMDMasking() {
    WorkspaceTester ws;
    TSM_ASSERT_THROWS("Characterisation test. This is not implemented.",
                      ws.clearMDMasking(), const std::runtime_error &);
  }

  void test_getSpecialCoordinateSystem_default() {
    WorkspaceTester ws;
    TSM_ASSERT_EQUALS("Should default to no special coordinate system.",
                      Mantid::Kernel::None, ws.getSpecialCoordinateSystem());
  }

  void test_getFirstPulseTime_getLastPulseTime() {
    WorkspaceTester ws;
    auto proton_charge = new TimeSeriesProperty<double>("proton_charge");
    DateAndTime startTime("2013-04-21T10:40:00");
    proton_charge->addValue(startTime, 1.0E-7);
    proton_charge->addValue(startTime + 1.0, 2.0E-7);
    proton_charge->addValue(startTime + 2.0, 3.0E-7);
    proton_charge->addValue(startTime + 3.0, 4.0E-7);
    ws.mutableRun().addLogData(proton_charge);

    TS_ASSERT_EQUALS(ws.getFirstPulseTime(), startTime);
    TS_ASSERT_EQUALS(ws.getLastPulseTime(), startTime + 3.0);
  }

  void test_getFirstPulseTime_getLastPulseTime_SNS1990bug() {
    WorkspaceTester ws;
    auto proton_charge = new TimeSeriesProperty<double>("proton_charge");
    DateAndTime startTime("1990-12-31T23:59:00");
    proton_charge->addValue(startTime, 1.0E-7);
    proton_charge->addValue(startTime + 1.0, 2.0E-7);
    ws.mutableRun().addLogData(proton_charge);

    // If fewer than 100 entries (unlikely to happen in reality), you just get
    // back the last one
    TS_ASSERT_EQUALS(ws.getFirstPulseTime(), startTime + 1.0);

    for (int i = 2; i < 62; ++i) {
      proton_charge->addValue(startTime + static_cast<double>(i), 1.0E-7);
    }
    TS_ASSERT_EQUALS(ws.getFirstPulseTime(),
                     DateAndTime("1991-01-01T00:00:00"));
  }

  void
  test_getFirstPulseTime_getLastPulseTime_throws_if_protoncharge_missing_or_empty() {
    WorkspaceTester ws;
    TS_ASSERT_THROWS(ws.getFirstPulseTime(), const std::runtime_error &);
    TS_ASSERT_THROWS(ws.getLastPulseTime(), const std::runtime_error &);
    ws.mutableRun().addLogData(new TimeSeriesProperty<double>("proton_charge"));
    TS_ASSERT_THROWS(ws.getFirstPulseTime(), const std::runtime_error &);
    TS_ASSERT_THROWS(ws.getLastPulseTime(), const std::runtime_error &);
  }

  void
  test_getFirstPulseTime_getLastPulseTime_throws_if_protoncharge_wrong_type() {
    WorkspaceTester ws;
    auto proton_charge = new TimeSeriesProperty<int>("proton_charge");
    proton_charge->addValue("2013-04-21T10:19:10", 1);
    proton_charge->addValue("2013-04-21T10:19:12", 2);
    ws.mutableRun().addLogData(proton_charge);
    TS_ASSERT_THROWS(ws.getFirstPulseTime(), const std::invalid_argument &);
    TS_ASSERT_THROWS(ws.getLastPulseTime(), const std::invalid_argument &);

    ws.mutableRun().addProperty(
        new PropertyWithValue<double>("proton_charge", 99.0), true);
    TS_ASSERT_THROWS(ws.getFirstPulseTime(), const std::invalid_argument &);
    TS_ASSERT_THROWS(ws.getLastPulseTime(), const std::invalid_argument &);
  }

  void test_getXMinMax() {
    double xmin, xmax;
    ws->getXMinMax(xmin, xmax);
    TS_ASSERT_EQUALS(xmin, 1.0);
    TS_ASSERT_EQUALS(xmax, 1.0);
    TS_ASSERT_EQUALS(ws->getXMin(), 1.0);
    TS_ASSERT_EQUALS(ws->getXMax(), 1.0);
  }

  void test_monitorWorkspace() {
    auto ws = std::make_shared<WorkspaceTester>();
    TSM_ASSERT("There should be no monitor workspace by default",
               !ws->monitorWorkspace())

    auto ws2 = std::make_shared<WorkspaceTester>();
    ws->setMonitorWorkspace(ws2);
    TSM_ASSERT_EQUALS("Monitor workspace not successfully set",
                      ws->monitorWorkspace(), ws2)

    ws->setMonitorWorkspace(std::shared_ptr<MatrixWorkspace>());
    TSM_ASSERT("Monitor workspace not successfully reset",
               !ws->monitorWorkspace())
  }

  void test_getXIndex() {
    WorkspaceTester ws;
    ws.initialize(1, 4, 3);
    auto &X = ws.dataX(0);
    X[0] = 1.0;
    X[1] = 2.0;
    X[2] = 3.0;
    X[3] = 4.0;

    auto ip = ws.getXIndex(0, 0.0, true);
    TS_ASSERT_EQUALS(ip.first, 0);
    TS_ASSERT_DELTA(ip.second, 0.0, 1e-15);

    ip = ws.getXIndex(0, 0.0, false);
    TS_ASSERT_EQUALS(ip.first, 4);
    TS_ASSERT_DELTA(ip.second, 0.0, 1e-15);

    ip = ws.getXIndex(0, 1.0, true);
    TS_ASSERT_EQUALS(ip.first, 0);
    TS_ASSERT_DELTA(ip.second, 0.0, 1e-15);

    ip = ws.getXIndex(0, 1.0, false);
    TS_ASSERT_EQUALS(ip.first, 4);
    TS_ASSERT_DELTA(ip.second, 0.0, 1e-15);

    ip = ws.getXIndex(0, 5.0, true);
    TS_ASSERT_EQUALS(ip.first, 4);
    TS_ASSERT_DELTA(ip.second, 0.0, 1e-15);

    ip = ws.getXIndex(0, 5.0, false);
    TS_ASSERT_EQUALS(ip.first, 3);
    TS_ASSERT_DELTA(ip.second, 0.0, 1e-15);

    ip = ws.getXIndex(0, 4.0, true);
    TS_ASSERT_EQUALS(ip.first, 4);
    TS_ASSERT_DELTA(ip.second, 0.0, 1e-15);

    ip = ws.getXIndex(0, 4.0, false);
    TS_ASSERT_EQUALS(ip.first, 3);
    TS_ASSERT_DELTA(ip.second, 0.0, 1e-15);

    ip = ws.getXIndex(0, 5.0, true, 5);
    TS_ASSERT_EQUALS(ip.first, 4);
    TS_ASSERT_DELTA(ip.second, 0.0, 1e-15);

    ip = ws.getXIndex(0, 5.0, false, 5);
    TS_ASSERT_EQUALS(ip.first, 4);
    TS_ASSERT_DELTA(ip.second, 0.0, 1e-15);

    ip = ws.getXIndex(0, 3.0, true, 5);
    TS_ASSERT_EQUALS(ip.first, 4);
    TS_ASSERT_DELTA(ip.second, 0.0, 1e-15);

    ip = ws.getXIndex(0, 3.0, false, 5);
    TS_ASSERT_EQUALS(ip.first, 4);
    TS_ASSERT_DELTA(ip.second, 0.0, 1e-15);

    ip = ws.getXIndex(0, 4.0, true, 5);
    TS_ASSERT_EQUALS(ip.first, 4);
    TS_ASSERT_DELTA(ip.second, 0.0, 1e-15);

    ip = ws.getXIndex(0, 4.0, false, 5);
    TS_ASSERT_EQUALS(ip.first, 4);
    TS_ASSERT_DELTA(ip.second, 0.0, 1e-15);

    ip = ws.getXIndex(0, 4.0, true, 4);
    TS_ASSERT_EQUALS(ip.first, 4);
    TS_ASSERT_DELTA(ip.second, 0.0, 1e-15);

    ip = ws.getXIndex(0, 4.0, false, 4);
    TS_ASSERT_EQUALS(ip.first, 4);
    TS_ASSERT_DELTA(ip.second, 0.0, 1e-15);

    ip = ws.getXIndex(0, 4.0, true, 3);
    TS_ASSERT_EQUALS(ip.first, 4);
    TS_ASSERT_DELTA(ip.second, 0.0, 1e-15);

    ip = ws.getXIndex(0, 4.0, false, 3);
    TS_ASSERT_EQUALS(ip.first, 3);
    TS_ASSERT_DELTA(ip.second, 0.0, 1e-15);

    ip = ws.getXIndex(0, 4.0, true);
    TS_ASSERT_EQUALS(ip.first, 4);
    TS_ASSERT_DELTA(ip.second, 0.0, 1e-15);

    ip = ws.getXIndex(0, 4.0, false);
    TS_ASSERT_EQUALS(ip.first, 3);
    TS_ASSERT_DELTA(ip.second, 0.0, 1e-15);

    ip = ws.getXIndex(0, 2.0, true, 3);
    TS_ASSERT_EQUALS(ip.first, 4);
    TS_ASSERT_DELTA(ip.second, 0.0, 1e-15);

    ip = ws.getXIndex(0, 2.0, false, 3);
    TS_ASSERT_EQUALS(ip.first, 3);
    TS_ASSERT_DELTA(ip.second, 0.0, 1e-15);

    ip = ws.getXIndex(0, 1.0, true, 3);
    TS_ASSERT_EQUALS(ip.first, 4);
    TS_ASSERT_DELTA(ip.second, 0.0, 1e-15);

    ip = ws.getXIndex(0, 1.0, false, 3);
    TS_ASSERT_EQUALS(ip.first, 3);
    TS_ASSERT_DELTA(ip.second, 0.0, 1e-15);

    ip = ws.getXIndex(0, 2.1, true);
    TS_ASSERT_EQUALS(ip.first, 1);
    TS_ASSERT_DELTA(ip.second, 0.1, 1e-15);

    ip = ws.getXIndex(0, 2.1, false);
    TS_ASSERT_EQUALS(ip.first, 2);
    TS_ASSERT_DELTA(ip.second, 0.9, 1e-15);
  }

  void test_getImage_0_width() {
    WorkspaceTester ws;
    ws.initialize(9, 2, 1);
    auto &X = ws.dataX(0);
    X[0] = 1.0;
    X[1] = 2.0;
    const size_t start = 0;
    const size_t stop = 8;
    size_t width = 0;
    TS_ASSERT_THROWS(ws.getImageY(start, stop, width),
                     const std::runtime_error &);
    width = 3;
    TS_ASSERT_THROWS_NOTHING(ws.getImageY(start, stop, width));
  }

  void test_getImage_wrong_start() {
    WorkspaceTester ws;
    ws.initialize(9, 2, 1);
    auto &X = ws.dataX(0);
    X[0] = 1.0;
    X[1] = 2.0;
    size_t start = 10;
    size_t stop = 8;
    size_t width = 3;
    TS_ASSERT_THROWS(ws.getImageY(start, stop, width),
                     const std::runtime_error &);
    start = 9;
    TS_ASSERT_THROWS(ws.getImageY(start, stop, width),
                     const std::runtime_error &);
    start = 0;
    TS_ASSERT_THROWS_NOTHING(ws.getImageY(start, stop, width));
  }

  void test_getImage_wrong_stop() {
    WorkspaceTester ws;
    ws.initialize(9, 2, 1);
    auto &X = ws.dataX(0);
    X[0] = 1.0;
    X[1] = 2.0;
    size_t start = 0;
    size_t stop = 18;
    size_t width = 3;
    TS_ASSERT_THROWS(ws.getImageY(start, stop, width),
                     const std::runtime_error &);
    stop = 9;
    TS_ASSERT_THROWS(ws.getImageY(start, stop, width),
                     const std::runtime_error &);
    stop = 8;
    TS_ASSERT_THROWS_NOTHING(ws.getImageY(start, stop, width));
  }

  void test_getImage_empty_set() {
    WorkspaceTester ws;
    ws.initialize(9, 2, 1);
    auto &X = ws.dataX(0);
    X[0] = 1.0;
    X[1] = 2.0;
    size_t start = 1;
    size_t stop = 0;
    size_t width = 1;
    TS_ASSERT_THROWS(ws.getImageY(start, stop, width),
                     const std::runtime_error &);
    stop = 1;
    TS_ASSERT_THROWS_NOTHING(ws.getImageY(start, stop, width));
  }

  void test_getImage_non_rectangular() {
    WorkspaceTester ws;
    ws.initialize(9, 2, 1);
    auto &X = ws.dataX(0);
    X[0] = 1.0;
    X[1] = 2.0;
    size_t start = 0;
    size_t stop = 7;
    size_t width = 3;
    TS_ASSERT_THROWS(ws.getImageY(start, stop, width),
                     const std::runtime_error &);
  }

  void test_getImage_wrong_indexStart() {
    WorkspaceTester ws;
    ws.initialize(9, 2, 1);
    auto &X = ws.dataX(0);
    X[0] = 1.0;
    X[1] = 2.0;
    const size_t start = 0;
    const size_t stop = 8;
    const size_t width = 3;
    double startX = 3;
    double endX = 4;
    TS_ASSERT_THROWS(ws.getImageY(start, stop, width, startX, endX),
                     const std::runtime_error &);

    WorkspaceTester wsh;
    wsh.initialize(9, 1, 1);
    startX = 2;
    endX = 2;
    TS_ASSERT_THROWS(wsh.getImageY(start, stop, width, startX, endX),
                     const std::runtime_error &);
  }

  void test_getImage_wrong_indexEnd() {
    WorkspaceTester ws;
    ws.initialize(9, 2, 1);
    auto &X = ws.dataX(0);
    X[0] = 1.0;
    X[1] = 2.0;
    const size_t start = 0;
    const size_t stop = 8;
    const size_t width = 3;
    double startX = 1.0;
    double endX = 0.0;
    TS_ASSERT_THROWS(ws.getImageY(start, stop, width, startX, endX),
                     const std::runtime_error &);

    WorkspaceTester wsh;
    wsh.initialize(9, 2, 2);
    auto &X1 = ws.dataX(0);
    X1[0] = 1.0;
    X1[1] = 2.0;
    startX = 1.0;
    endX = 0.0;
    TS_ASSERT_THROWS(wsh.getImageY(start, stop, width, startX, endX),
                     const std::runtime_error &);
  }

  void test_getImage_single_bin_histo() {
    WorkspaceTester ws;
    ws.initialize(9, 2, 1);
    auto &X = ws.dataX(0);
    X[0] = 1.0;
    X[1] = 2.0;
    for (size_t i = 0; i < ws.getNumberHistograms(); ++i) {
      ws.dataY(i)[0] = static_cast<double>(i + 1);
    }
    const size_t start = 0;
    const size_t stop = 8;
    const size_t width = 3;
    double startX = 0;
    double endX = 3;
    Mantid::API::MantidImage_sptr image;
    TS_ASSERT_THROWS_NOTHING(
        image = ws.getImageY(start, stop, width, startX, endX));
    if (!image)
      return;
    TS_ASSERT_EQUALS(image->size(), 3);
    TS_ASSERT_EQUALS((*image)[0].size(), 3);
    TS_ASSERT_EQUALS((*image)[1].size(), 3);
    TS_ASSERT_EQUALS((*image)[2].size(), 3);

    TS_ASSERT_EQUALS((*image)[0][0], 1);
    TS_ASSERT_EQUALS((*image)[0][1], 2);
    TS_ASSERT_EQUALS((*image)[0][2], 3);
    TS_ASSERT_EQUALS((*image)[1][0], 4);
    TS_ASSERT_EQUALS((*image)[1][1], 5);
    TS_ASSERT_EQUALS((*image)[1][2], 6);
    TS_ASSERT_EQUALS((*image)[2][0], 7);
    TS_ASSERT_EQUALS((*image)[2][1], 8);
    TS_ASSERT_EQUALS((*image)[2][2], 9);
  }

  void test_getImage_single_bin_points() {
    WorkspaceTester ws;
    ws.initialize(9, 1, 1);
    auto &X = ws.dataX(0);
    X[0] = 1.0;
    for (size_t i = 0; i < ws.getNumberHistograms(); ++i) {
      ws.dataY(i)[0] = static_cast<double>(i + 1);
    }
    const size_t start = 0;
    const size_t stop = 8;
    const size_t width = 3;
    double startX = 1;
    double endX = 1;
    Mantid::API::MantidImage_sptr image;
    TS_ASSERT_THROWS_NOTHING(
        image = ws.getImageY(start, stop, width, startX, endX));
    if (!image)
      return;
    TS_ASSERT_EQUALS(image->size(), 3);
    TS_ASSERT_EQUALS((*image)[0].size(), 3);
    TS_ASSERT_EQUALS((*image)[1].size(), 3);
    TS_ASSERT_EQUALS((*image)[2].size(), 3);

    TS_ASSERT_EQUALS((*image)[0][0], 1);
    TS_ASSERT_EQUALS((*image)[0][1], 2);
    TS_ASSERT_EQUALS((*image)[0][2], 3);
    TS_ASSERT_EQUALS((*image)[1][0], 4);
    TS_ASSERT_EQUALS((*image)[1][1], 5);
    TS_ASSERT_EQUALS((*image)[1][2], 6);
    TS_ASSERT_EQUALS((*image)[2][0], 7);
    TS_ASSERT_EQUALS((*image)[2][1], 8);
    TS_ASSERT_EQUALS((*image)[2][2], 9);
  }

  void test_getImage_multi_bin_histo() {
    WorkspaceTester ws;
    ws.initialize(9, 4, 3);
    auto &X = ws.dataX(0);
    X[0] = 1.0;
    X[1] = 2.0;
    X[2] = 3.0;
    X[3] = 4.0;
    for (size_t i = 0; i < ws.getNumberHistograms(); ++i) {
      ws.dataY(i)[0] = static_cast<double>(i + 1);
      ws.dataY(i)[1] = static_cast<double>(i + 2);
      ws.dataY(i)[2] = static_cast<double>(i + 3);
    }
    const size_t start = 0;
    const size_t stop = 8;
    const size_t width = 3;
    Mantid::API::MantidImage_sptr image;
    TS_ASSERT_THROWS_NOTHING(image = ws.getImageY(start, stop, width));
    if (!image)
      return;
    TS_ASSERT_EQUALS(image->size(), 3);
    TS_ASSERT_EQUALS((*image)[0].size(), 3);
    TS_ASSERT_EQUALS((*image)[1].size(), 3);
    TS_ASSERT_EQUALS((*image)[2].size(), 3);

    TS_ASSERT_EQUALS((*image)[0][0], 6);
    TS_ASSERT_EQUALS((*image)[0][1], 9);
    TS_ASSERT_EQUALS((*image)[0][2], 12);
    TS_ASSERT_EQUALS((*image)[1][0], 15);
    TS_ASSERT_EQUALS((*image)[1][1], 18);
    TS_ASSERT_EQUALS((*image)[1][2], 21);
    TS_ASSERT_EQUALS((*image)[2][0], 24);
    TS_ASSERT_EQUALS((*image)[2][1], 27);
    TS_ASSERT_EQUALS((*image)[2][2], 30);
  }

  void test_getImage_multi_bin_points() {
    WorkspaceTester ws;
    ws.initialize(9, 3, 3);
    auto &X = ws.dataX(0);
    X[0] = 1.0;
    X[1] = 2.0;
    X[2] = 3.0;
    for (size_t i = 0; i < ws.getNumberHistograms(); ++i) {
      ws.dataY(i)[0] = static_cast<double>(i + 1);
      ws.dataY(i)[1] = static_cast<double>(i + 2);
      ws.dataY(i)[2] = static_cast<double>(i + 3);
    }
    const size_t start = 0;
    const size_t stop = 8;
    const size_t width = 3;
    Mantid::API::MantidImage_sptr image;
    TS_ASSERT_THROWS_NOTHING(image = ws.getImageY(start, stop, width));
    if (!image)
      return;
    TS_ASSERT_EQUALS(image->size(), 3);
    TS_ASSERT_EQUALS((*image)[0].size(), 3);
    TS_ASSERT_EQUALS((*image)[1].size(), 3);
    TS_ASSERT_EQUALS((*image)[2].size(), 3);

    TS_ASSERT_EQUALS((*image)[0][0], 6);
    TS_ASSERT_EQUALS((*image)[0][1], 9);
    TS_ASSERT_EQUALS((*image)[0][2], 12);
    TS_ASSERT_EQUALS((*image)[1][0], 15);
    TS_ASSERT_EQUALS((*image)[1][1], 18);
    TS_ASSERT_EQUALS((*image)[1][2], 21);
    TS_ASSERT_EQUALS((*image)[2][0], 24);
    TS_ASSERT_EQUALS((*image)[2][1], 27);
    TS_ASSERT_EQUALS((*image)[2][2], 30);
  }

  void test_setImage_too_large() {
    auto image = createImage(2, 3);
    WorkspaceTester ws;
    ws.initialize(2, 2, 1);
    TS_ASSERT_THROWS(ws.setImageY(*image), const std::runtime_error &);
  }

  void test_setImage_not_single_bin() {
    auto image = createImage(2, 3);
    WorkspaceTester ws;
    ws.initialize(20, 3, 2);
    TS_ASSERT_THROWS(ws.setImageY(*image), const std::runtime_error &);
  }

  void test_setImageY() {
    auto image = createImage(2, 3);
    WorkspaceTester ws;
    ws.initialize(6, 2, 1);
    TS_ASSERT_THROWS_NOTHING(ws.setImageY(*image));
    TS_ASSERT_EQUALS(ws.readY(0)[0], 1);
    TS_ASSERT_EQUALS(ws.readY(1)[0], 2);
    TS_ASSERT_EQUALS(ws.readY(2)[0], 3);
    TS_ASSERT_EQUALS(ws.readY(3)[0], 4);
    TS_ASSERT_EQUALS(ws.readY(4)[0], 5);
    TS_ASSERT_EQUALS(ws.readY(5)[0], 6);
  }

  void test_setImageE() {
    auto image = createImage(2, 3);
    WorkspaceTester ws;
    ws.initialize(6, 2, 1);
    TS_ASSERT_THROWS_NOTHING(ws.setImageE(*image));
    TS_ASSERT_EQUALS(ws.readE(0)[0], 1);
    TS_ASSERT_EQUALS(ws.readE(1)[0], 2);
    TS_ASSERT_EQUALS(ws.readE(2)[0], 3);
    TS_ASSERT_EQUALS(ws.readE(3)[0], 4);
    TS_ASSERT_EQUALS(ws.readE(4)[0], 5);
    TS_ASSERT_EQUALS(ws.readE(5)[0], 6);
  }

  void test_setImageY_start() {
    auto image = createImage(2, 3);
    WorkspaceTester ws;
    ws.initialize(9, 2, 1);
    TS_ASSERT_THROWS_NOTHING(ws.setImageY(*image, 3));
    TS_ASSERT_EQUALS(ws.readY(3)[0], 1);
    TS_ASSERT_EQUALS(ws.readY(4)[0], 2);
    TS_ASSERT_EQUALS(ws.readY(5)[0], 3);
    TS_ASSERT_EQUALS(ws.readY(6)[0], 4);
    TS_ASSERT_EQUALS(ws.readY(7)[0], 5);
    TS_ASSERT_EQUALS(ws.readY(8)[0], 6);
  }

  void test_setImageE_start() {
    auto image = createImage(2, 3);
    WorkspaceTester ws;
    ws.initialize(9, 2, 1);
    TS_ASSERT_THROWS_NOTHING(ws.setImageE(*image, 2));
    TS_ASSERT_EQUALS(ws.readE(2)[0], 1);
    TS_ASSERT_EQUALS(ws.readE(3)[0], 2);
    TS_ASSERT_EQUALS(ws.readE(4)[0], 3);
    TS_ASSERT_EQUALS(ws.readE(5)[0], 4);
    TS_ASSERT_EQUALS(ws.readE(6)[0], 5);
    TS_ASSERT_EQUALS(ws.readE(7)[0], 6);
  }

  /**
   * Test declaring an input workspace and retrieving as const_sptr or sptr
   */
  void testGetProperty_const_sptr() {
    const std::string wsName = "InputWorkspace";
    MatrixWorkspace_sptr wsInput = std::make_shared<WorkspaceTester>();
    PropertyManagerHelper manager;
    manager.declareProperty(wsName, wsInput, Direction::Input);

    // Check property can be obtained as const_sptr or sptr
    MatrixWorkspace_const_sptr wsConst;
    MatrixWorkspace_sptr wsNonConst;
    TS_ASSERT_THROWS_NOTHING(
        wsConst = manager.getValue<MatrixWorkspace_const_sptr>(wsName));
    TS_ASSERT(wsConst != nullptr);
    TS_ASSERT_THROWS_NOTHING(
        wsNonConst = manager.getValue<MatrixWorkspace_sptr>(wsName));
    TS_ASSERT(wsNonConst != nullptr);
    TS_ASSERT_EQUALS(wsConst, wsNonConst);

    // Check TypedValue can be cast to const_sptr or to sptr
    PropertyManagerHelper::TypedValue val(manager, wsName);
    MatrixWorkspace_const_sptr wsCastConst;
    MatrixWorkspace_sptr wsCastNonConst;
    TS_ASSERT_THROWS_NOTHING(wsCastConst = (MatrixWorkspace_const_sptr)val);
    TS_ASSERT(wsCastConst != nullptr);
    TS_ASSERT_THROWS_NOTHING(wsCastNonConst = (MatrixWorkspace_sptr)val);
    TS_ASSERT(wsCastNonConst != nullptr);
    TS_ASSERT_EQUALS(wsCastConst, wsCastNonConst);
  }

  void test_x_uncertainty_can_be_set() {
    // Arrange
    WorkspaceTester ws;
    const size_t numspec = 4;
    const size_t j = 3;
    const size_t k = j;
    ws.initialize(numspec, j, k);

    double values[3] = {10, 11, 17};
    size_t workspaceIndexWithDx[3] = {0, 1, 2};

    Mantid::MantidVec dxSpec0(j, values[0]);
    auto dxSpec1 =
        Kernel::make_cow<Mantid::HistogramData::HistogramDx>(j, values[1]);
    auto dxSpec2 = std::make_shared<Mantid::HistogramData::HistogramDx>(
        Mantid::MantidVec(j, values[2]));

    // Act
    for (size_t spec = 0; spec < numspec; ++spec) {
      TSM_ASSERT("Should not have any x resolution values", !ws.hasDx(spec));
    }
    ws.dataDx(workspaceIndexWithDx[0]) = dxSpec0;
    ws.setSharedDx(workspaceIndexWithDx[1], dxSpec1);
    ws.setSharedDx(workspaceIndexWithDx[2], dxSpec2);

    // Assert
    auto compareValue = [&values](double data, size_t index) {
      return data == values[index];
    };
    for (auto &index : workspaceIndexWithDx) {
      TSM_ASSERT("Should have x resolution values", ws.hasDx(index));
      TSM_ASSERT_EQUALS("Should have a length of 3", ws.dataDx(index).size(),
                        j);
      auto compareValueForSpecificWorkspaceIndex =
          std::bind(compareValue, std::placeholders::_1, index);

      auto &dataDx = ws.dataDx(index);
      TSM_ASSERT("dataDx should allow access to the spectrum",
                 std::all_of(std::begin(dataDx), std::end(dataDx),
                             compareValueForSpecificWorkspaceIndex));

      auto &readDx = ws.readDx(index);
      TSM_ASSERT("readDx should allow access to the spectrum",
                 std::all_of(std::begin(readDx), std::end(readDx),
                             compareValueForSpecificWorkspaceIndex));

      auto refDx = ws.sharedDx(index);
      TSM_ASSERT("readDx should allow access to the spectrum",
                 std::all_of(std::begin(*refDx), std::end(*refDx),
                             compareValueForSpecificWorkspaceIndex));
    }

    TSM_ASSERT("Should not have any x resolution values", !ws.hasDx(3));
  }

  void test_scanning() {
    // Set up 2 workspaces to be merged
    auto ws1 = makeWorkspaceWithDetectors(1, 1);
    auto ws2 = makeWorkspaceWithDetectors(1, 1);
    auto &detInfo1 = ws1->mutableDetectorInfo();
    auto &detInfo2 = ws2->mutableDetectorInfo();
    auto &compInfo1 = ws1->mutableComponentInfo();
    auto &compInfo2 = ws2->mutableComponentInfo();

    detInfo1.setPosition(0, {1, 0, 0});
    detInfo2.setPosition(0, {2, 0, 0});
    compInfo1.setScanInterval({10, 20});
    compInfo2.setScanInterval({20, 30});

    // Merge
    auto merged = WorkspaceFactory::Instance().create(ws1, 2);
    merged->mutableComponentInfo().merge(ws2->componentInfo());

    // Setting IndexInfo without spectrum definitions will set up a 1:1 mapping
    // such that each spectrum corresponds to 1 time index of a detector.
    merged->setIndexInfo(IndexInfo(merged->getNumberHistograms()));

    const auto &specInfo = merged->spectrumInfo();
    TS_ASSERT(specInfo.hasDetectors(0));
    TS_ASSERT(specInfo.hasDetectors(1));
    // This is the order we get currently from the default mapping, but it is
    // not guaranteed by the interface and might change.
    TS_ASSERT_EQUALS(specInfo.position(0), V3D(1, 0, 0));
    TS_ASSERT_EQUALS(specInfo.position(1), V3D(2, 0, 0));

    TS_ASSERT_THROWS_NOTHING(specInfo.detector(0));
    const auto &det = specInfo.detector(0);
    // Failing legacy methods (use DetectorInfo/SpectrumInfo instead):
    TS_ASSERT_THROWS(det.getPos(), const std::runtime_error &);
    TS_ASSERT_THROWS(det.getRelativePos(), const std::runtime_error &);
    TS_ASSERT_THROWS(det.getRotation(), const std::runtime_error &);
    TS_ASSERT_THROWS(det.getRelativeRot(), const std::runtime_error &);
    TS_ASSERT_THROWS(det.getPhi(), const std::runtime_error &);
    // Failing methods, currently without replacement:
    TS_ASSERT_THROWS(det.solidAngle(V3D(0, 0, 0)), const std::runtime_error &);
    BoundingBox bb;
    TS_ASSERT_THROWS(det.getBoundingBox(bb), const std::runtime_error &);
    // Moving parent not possible since non-detector components do not have time
    // indices and thus DetectorInfo cannot tell which set of detector positions
    // to adjust.

    auto &compInfo = merged->mutableComponentInfo();

    // Try to move the parent
    TS_ASSERT_THROWS(compInfo.setPosition(compInfo.parent(compInfo.indexOf(
                                              det.getComponentID())),
                                          V3D(1, 2, 3)),
                     const std::runtime_error &);
    // Try to rotate the parent
    TS_ASSERT_THROWS(compInfo.setRotation(compInfo.parent(compInfo.indexOf(
                                              det.getComponentID())),
                                          Quat(1, 2, 3, 4)),
                     const std::runtime_error &);
  }

  void test_legacy_setting_spectrum_numbers_with_MPI() {
    ParallelTestHelpers::runParallel(
        run_legacy_setting_spectrum_numbers_with_MPI);
  }

  void test_detectorSignedTwoTheta() {
    checkDetectorSignedTwoTheta(Geometry::Y, {{1., 1., -1., -1.}});
    checkDetectorSignedTwoTheta(Geometry::X, {{1., -1., -1., 1.}});
  }

  void
  test_that_yIndexOfX_returns_the_ascending_index_for_a_histogram_workspace_when_provided_an_x_value_corresponding_to_the_first_bin() {
    std::vector<double> const xValues{1.0, 2.0, 3.0, 4.0};
    auto const workspace = getWorkspaceWithPopulatedX(1, 4, 3, xValues);

    TS_ASSERT_EQUALS(workspace.yIndexOfX(1.3), 0);
  }

  void
  test_that_yIndexOfX_returns_the_ascending_index_for_a_histogram_workspace_when_provided_an_x_value_which_is_a_bin_boundary() {
    std::vector<double> const xValues{1.0, 2.0, 3.0, 4.0};
    auto const workspace = getWorkspaceWithPopulatedX(1, 4, 3, xValues);

    TS_ASSERT_EQUALS(workspace.yIndexOfX(2.0), 0);
  }

  void
  test_that_yIndexOfX_returns_the_ascending_index_for_a_histogram_workspace_when_provided_an_x_value_which_is_mid_range() {
    std::vector<double> const xValues{1.0, 2.0, 3.0, 4.0};
    auto const workspace = getWorkspaceWithPopulatedX(1, 4, 3, xValues);

    TS_ASSERT_EQUALS(workspace.yIndexOfX(2.5), 1);
  }

  void
  test_that_yIndexOfX_returns_the_ascending_index_for_a_histogram_workspace_when_provided_an_x_value_which_is_only_just_in_the_second_bin() {
    std::vector<double> const xValues{1.0, 2.0, 3.0, 4.0};
    auto const workspace = getWorkspaceWithPopulatedX(1, 4, 3, xValues);

    TS_ASSERT_EQUALS(workspace.yIndexOfX(2.001), 1);
  }

  void
  test_that_yIndexOfX_returns_the_ascending_index_for_a_histogram_workspace_when_provided_an_x_value_which_is_in_the_last_bin() {
    std::vector<double> const xValues{1.0, 2.0, 3.0, 4.0};
    auto const workspace = getWorkspaceWithPopulatedX(1, 4, 3, xValues);

    TS_ASSERT_EQUALS(workspace.yIndexOfX(3.1), 2);
  }

  void
  test_that_yIndexOfX_returns_the_ascending_index_for_a_histogram_workspace_when_provided_an_x_value_which_is_the_last_bin_boundary() {
    std::vector<double> const xValues{1.0, 2.0, 3.0, 4.0};
    auto const workspace = getWorkspaceWithPopulatedX(1, 4, 3, xValues);

    TS_ASSERT_EQUALS(workspace.yIndexOfX(4.0), 2);
  }

  void
  test_that_yIndexOfX_throws_when_provided_an_index_which_is_out_of_range_for_ascending_x_values() {
    std::vector<double> const xValues{1.0, 2.0, 3.0, 4.0};
    auto const workspace = getWorkspaceWithPopulatedX(1, 4, 3, xValues);

    TS_ASSERT_EQUALS(workspace.getNumberHistograms(), 1);
    TS_ASSERT_THROWS(workspace.yIndexOfX(2.5, 1), const std::out_of_range &);
    TS_ASSERT_THROWS(workspace.yIndexOfX(2.5, -1), const std::out_of_range &);
  }

  void
  test_that_yIndexOfX_throws_when_provided_x_values_which_are_out_of_range_for_ascending_x_values() {
    std::vector<double> const xValues{1.0, 2.0, 3.0, 4.0};
    auto const workspace = getWorkspaceWithPopulatedX(1, 4, 3, xValues);

    TS_ASSERT_THROWS(workspace.yIndexOfX(5.), const std::out_of_range &);
    TS_ASSERT_THROWS(workspace.yIndexOfX(0.), const std::out_of_range &);
  }

  void
  test_that_yIndexOfX_returns_the_descending_index_for_a_histogram_workspace_when_provided_an_x_value_corresponding_to_the_first_bin() {
    std::vector<double> const xValues{5.3, 4.3, 3.3, 2.3};
    auto const workspace = getWorkspaceWithPopulatedX(1, 4, 3, xValues);

    TS_ASSERT_EQUALS(workspace.yIndexOfX(5.2), 0);
  }

  void
  test_that_yIndexOfX_returns_the_descending_index_for_a_histogram_workspace_when_provided_an_x_value_corresponding_to_the_first_boundary() {
    std::vector<double> const xValues{5.3, 4.3, 3.3, 2.3};
    auto const workspace = getWorkspaceWithPopulatedX(1, 4, 3, xValues);

    TS_ASSERT_EQUALS(workspace.yIndexOfX(5.3), 0)
  }

  void
  test_that_yIndexOfX_returns_the_descending_index_for_a_histogram_workspace_when_provided_an_x_value_which_is_a_bin_boundary() {
    std::vector<double> const xValues{5.3, 4.3, 3.3, 2.3};
    auto const workspace = getWorkspaceWithPopulatedX(1, 4, 3, xValues);

    TS_ASSERT_EQUALS(workspace.yIndexOfX(4.3), 0);
  }

  void
  test_that_yIndexOfX_returns_the_descending_index_for_a_histogram_workspace_when_provided_an_x_value_which_is_mid_range() {
    std::vector<double> const xValues{5.3, 4.3, 3.3, 2.3};
    auto const workspace = getWorkspaceWithPopulatedX(1, 4, 3, xValues);

    TS_ASSERT_EQUALS(workspace.yIndexOfX(3.8), 1);
  }

  void
  test_that_yIndexOfX_returns_the_descending_index_for_a_histogram_workspace_when_provided_an_x_value_which_is_only_just_in_the_second_bin() {
    std::vector<double> const xValues{5.3, 4.3, 3.3, 2.3};
    auto const workspace = getWorkspaceWithPopulatedX(1, 4, 3, xValues);

    TS_ASSERT_EQUALS(workspace.yIndexOfX(std::nextafter(3.3, 10.0)), 1);
  }

  void
  test_that_yIndexOfX_returns_the_descending_index_for_a_histogram_workspace_when_provided_an_x_value_which_is_in_the_last_bin() {
    std::vector<double> const xValues{5.3, 4.3, 3.3, 2.3};
    auto const workspace = getWorkspaceWithPopulatedX(1, 4, 3, xValues);

    TS_ASSERT_EQUALS(workspace.yIndexOfX(3.1), 2);
  }

  void
  test_that_yIndexOfX_returns_the_descending_index_for_a_histogram_workspace_when_provided_an_x_value_which_is_the_last_bin_boundary() {
    std::vector<double> const xValues{5.3, 4.3, 3.3, 2.3};
    auto const workspace = getWorkspaceWithPopulatedX(1, 4, 3, xValues);

    TS_ASSERT_EQUALS(workspace.yIndexOfX(2.3), 2);
  }

  void
  test_that_yIndexOfX_throws_when_provided_an_index_which_is_out_of_range_for_descending_x_values() {
    std::vector<double> const xValues{5.3, 4.3, 3.3, 2.3};
    auto const workspace = getWorkspaceWithPopulatedX(1, 4, 3, xValues);

    TS_ASSERT_EQUALS(workspace.getNumberHistograms(), 1);
    TS_ASSERT_THROWS(workspace.yIndexOfX(2.5, 1), const std::out_of_range &);
    TS_ASSERT_THROWS(workspace.yIndexOfX(2.5, -1), const std::out_of_range &);
  }

  void
  test_that_yIndexOfX_throws_when_provided_x_values_which_are_out_of_range_for_descending_x_values() {
    std::vector<double> const xValues{5.3, 4.3, 3.3, 2.3};
    auto const workspace = getWorkspaceWithPopulatedX(1, 4, 3, xValues);

    TS_ASSERT_THROWS(workspace.yIndexOfX(std::nextafter(5.3, 10.0)),
                     const std::out_of_range &);
    TS_ASSERT_THROWS(workspace.yIndexOfX(5.4), const std::out_of_range &);
    TS_ASSERT_THROWS(workspace.yIndexOfX(std::nextafter(2.3, 0.0)),
                     const std::out_of_range &);
    TS_ASSERT_THROWS(workspace.yIndexOfX(0.), const std::out_of_range &);
  }

  void
  test_that_yIndexOfX_returns_the_correct_index_for_a_histogram_workspace_when_the_x_bins_are_in_ascending_order() {
    std::vector<double> const xValues{1.0, 2.0, 3.0, 4.0};
    auto const workspace = getWorkspaceWithPopulatedX(1, 4, 3, xValues);

    TS_ASSERT_EQUALS(workspace.yIndexOfX(3.0), 1);
  }

  void
  test_that_yIndexOfX_returns_the_correct_index_for_a_histogram_workspace_when_the_x_bins_are_in_descending_order() {
    std::vector<double> const xValues{4.0, 3.0, 2.0, 1.0};
    auto const workspace = getWorkspaceWithPopulatedX(1, 4, 3, xValues);

    TS_ASSERT_EQUALS(workspace.yIndexOfX(3.0), 0);
  }

  void
  test_that_yIndexOfX_returns_the_correct_index_for_a_nonHistogram_workspace_when_the_x_bins_are_in_ascending_order() {
    std::vector<double> const xValues{1.0, 2.0, 3.0, 4.0};
    auto const workspace = getWorkspaceWithPopulatedX(1, 4, 4, xValues);

    TS_ASSERT_EQUALS(workspace.yIndexOfX(3.0), 2);
  }

  void
  test_that_yIndexOfX_returns_the_correct_index_for_a_nonHistogram_workspace_when_the_x_bins_are_in_descending_order() {
    std::vector<double> const xValues{4.0, 3.0, 2.0, 1.0};
    auto const workspace = getWorkspaceWithPopulatedX(1, 4, 4, xValues);

    TS_ASSERT_EQUALS(workspace.yIndexOfX(3.0), 1);
  }

  void
  test_that_yIndexOfX_returns_the_ascending_index_for_a_nonHistogram_workspace_when_passed_an_x_value_within_a_tolerance() {
    std::vector<double> const xValues{1.0, 1.9997, 3.0, 4.0};
    auto const workspace = getWorkspaceWithPopulatedX(1, 4, 4, xValues);

    TS_ASSERT_EQUALS(workspace.yIndexOfX(2.0, 0, 0.0005), 1);
  }

  void
  test_that_yIndexOfX_returns_the_ascending_index_for_a_nonHistogram_workspace_when_passed_an_x_value_within_a_tolerance2() {
    std::vector<double> const xValues{1.0, 1.9997, 3.0, 4.0};
    auto const workspace = getWorkspaceWithPopulatedX(1, 4, 4, xValues);

    TS_ASSERT_EQUALS(workspace.yIndexOfX(1.9995, 0, 0.0005), 1);
  }

  void
  test_that_yIndexOfX_returns_the_ascending_index_for_a_nonHistogram_workspace_when_passed_an_x_value_with_no_tolerance() {
    std::vector<double> const xValues{1.0, 1.9997, 3.0, 4.0};
    auto const workspace = getWorkspaceWithPopulatedX(1, 4, 4, xValues);

    TS_ASSERT_EQUALS(workspace.yIndexOfX(1.9997), 1);
  }

  void
  test_that_yIndexOfX_throws_for_a_nonHistogram_workspace_when_passed_an_x_value_just_outside_a_tolerance() {
    std::vector<double> const xValues{1.0, 1.9997, 3.0, 4.0};
    auto const workspace = getWorkspaceWithPopulatedX(1, 4, 4, xValues);

    TS_ASSERT_THROWS(workspace.yIndexOfX(2.0, 0, 0.0002),
                     const std::invalid_argument &);
  }

  void
  test_that_yIndexOfX_throws_for_a_nonHistogram_workspace_when_passed_an_x_value_just_outside_a_tolerance2() {
    std::vector<double> const xValues{1.0, 1.9997, 3.0, 4.0};
    auto const workspace = getWorkspaceWithPopulatedX(1, 4, 4, xValues);

    TS_ASSERT_THROWS(workspace.yIndexOfX(1.9992, 0, 0.0002),
                     const std::invalid_argument &);
  }

  void
  test_that_yIndexOfX_throw_when_passed_an_invalid_x_value_with_no_tolerance_for_an_ascending_vector_of_x_values() {
    std::vector<double> const xValues{1.0, 1.9997, 3.0, 4.0};
    auto const workspace = getWorkspaceWithPopulatedX(1, 4, 4, xValues);

    TS_ASSERT_THROWS(workspace.yIndexOfX(3.5), const std::invalid_argument &);
  }

  void
  test_that_yIndexOfX_returns_the_descending_index_for_a_nonHistogram_workspace_when_passed_an_x_value_within_a_tolerance() {
    std::vector<double> const xValues{4.0, 3.0, 1.9997, 1.0};
    auto const workspace = getWorkspaceWithPopulatedX(1, 4, 4, xValues);

    TS_ASSERT_EQUALS(workspace.yIndexOfX(1.9994, 0, 0.0005), 2);
  }

  void
  test_that_yIndexOfX_returns_the_descending_index_for_a_nonHistogram_workspace_when_passed_an_x_value_within_a_tolerance2() {
    std::vector<double> const xValues{4.0, 3.0, 1.9997, 1.0};
    auto const workspace = getWorkspaceWithPopulatedX(1, 4, 4, xValues);

    TS_ASSERT_EQUALS(workspace.yIndexOfX(1.9999, 0, 0.0005), 2);
  }

  void
  test_that_yIndexOfX_returns_the_descending_index_for_a_nonHistogram_workspace_when_passed_an_x_value_with_no_tolerance() {
    std::vector<double> const xValues{4.0, 3.0, 1.9997, 1.0};
    auto const workspace = getWorkspaceWithPopulatedX(1, 4, 4, xValues);

    TS_ASSERT_EQUALS(workspace.yIndexOfX(1.9997), 2);
  }

  void
  test_that_yIndexOfX_returns_the_descending_index_for_a_nonHistogram_workspace_when_passed_an_x_value_just_outside_a_tolerance() {
    std::vector<double> const xValues{4.0, 3.0, 1.9997, 1.0};
    auto const workspace = getWorkspaceWithPopulatedX(1, 4, 4, xValues);

    TS_ASSERT_THROWS(workspace.yIndexOfX(1.9994, 0, 0.0002),
                     const std::invalid_argument &);
  }

  void
  test_that_yIndexOfX_returns_the_descending_index_for_a_nonHistogram_workspace_when_passed_an_x_value_just_outside_a_tolerance2() {
    std::vector<double> const xValues{4.0, 3.0, 1.9997, 1.0};
    auto const workspace = getWorkspaceWithPopulatedX(1, 4, 4, xValues);

    TS_ASSERT_THROWS(workspace.yIndexOfX(2.0, 0, 0.0002),
                     const std::invalid_argument &);
  }

  void
  test_that_yIndexOfX_throw_when_passed_an_invalid_x_value_with_no_tolerance_for_a_descending_vector_of_x_values() {
    std::vector<double> const xValues{4.0, 3.0, 1.9997, 1.0};
    auto const workspace = getWorkspaceWithPopulatedX(1, 4, 4, xValues);

    TS_ASSERT_THROWS(workspace.yIndexOfX(3.5), const std::invalid_argument &);
  }

  void
  test_YUnitLabel_Correct_For_Distribution_Workspace_Custom_m_YUnitLabel_Not_Set() {
    auto testWS = generateTestWorkspaceWithDistributionAndLabelSet(true, "");
    TS_ASSERT_EQUALS(testWS->YUnitLabel(false, false),
                     "Counts per microsecond");
    TS_ASSERT_EQUALS(testWS->YUnitLabel(false, true), "Counts per microsecond");
    TS_ASSERT_EQUALS(testWS->YUnitLabel(true, false),
                     "Counts ($\\mu s$)$^{-1}$");
    TS_ASSERT_EQUALS(testWS->YUnitLabel(true, true),
                     "Counts ($\\mu s$)$^{-1}$");
  }

  void
  test_YUnitLabel_Correct_For_Distribution_Workspace_Custom_m_YUnitLabel_Set() {
    auto testWS =
        generateTestWorkspaceWithDistributionAndLabelSet(true, "Custom Label");
    TS_ASSERT_EQUALS(testWS->YUnitLabel(false, false), "Custom Label");
    TS_ASSERT_EQUALS(testWS->YUnitLabel(false, true), "Custom Label");
    TS_ASSERT_EQUALS(testWS->YUnitLabel(true, false), "Custom Label");
    TS_ASSERT_EQUALS(testWS->YUnitLabel(true, true), "Custom Label");
  }

  void
  test_YUnitLabel_Correct_For_Non_Distribution_Workspace_Custom_m_YUnitLabel_Not_Set() {
    auto testWS = generateTestWorkspaceWithDistributionAndLabelSet(false, "");
    TS_ASSERT_EQUALS(testWS->YUnitLabel(false, false), "Counts");
    TS_ASSERT_EQUALS(testWS->YUnitLabel(false, true), "Counts per microsecond");
    TS_ASSERT_EQUALS(testWS->YUnitLabel(true, false), "Counts");
    TS_ASSERT_EQUALS(testWS->YUnitLabel(true, true),
                     "Counts ($\\mu s$)$^{-1}$");
  }

  void
  test_YUnitLabel_Correct_For_Non_Distribution_Workspace_Custom_m_YUnitLabel_Set() {
    auto testWS =
        generateTestWorkspaceWithDistributionAndLabelSet(false, "Custom Label");
    TS_ASSERT_EQUALS(testWS->YUnitLabel(false, false), "Custom Label");
    TS_ASSERT_EQUALS(testWS->YUnitLabel(false, true),
                     "Custom Label per microsecond");
    TS_ASSERT_EQUALS(testWS->YUnitLabel(true, false), "Custom Label");
    TS_ASSERT_EQUALS(testWS->YUnitLabel(true, true),
                     "Custom Label ($\\mu s$)$^{-1}$");
  }

  void test_YUnitLabel_Correct_For_Empty_Y_Labels() {
    auto testWS = std::make_shared<WorkspaceTester>();
    testWS->initialize(1, 2, 1);
    testWS->setDistribution(false);
    testWS->getAxis(0)->setUnit("TOF");
    TS_ASSERT_EQUALS(testWS->YUnitLabel(false, false), "");
    TS_ASSERT_EQUALS(testWS->YUnitLabel(false, true), " per microsecond");
    TS_ASSERT_EQUALS(testWS->YUnitLabel(true, false), "");
    TS_ASSERT_EQUALS(testWS->YUnitLabel(true, true), "($\\mu s$)$^{-1}$");
  }

  void test_YUnitLabel_Correct_For_Empty_X_And_Y_Labels() {
    auto testWS = std::make_shared<WorkspaceTester>();
    testWS->initialize(1, 2, 1);
    testWS->setDistribution(false);
    TS_ASSERT_EQUALS(testWS->YUnitLabel(false, false), "");
    TS_ASSERT_EQUALS(testWS->YUnitLabel(false, true), "");
    TS_ASSERT_EQUALS(testWS->YUnitLabel(true, false), "");
    TS_ASSERT_EQUALS(testWS->YUnitLabel(true, true), "");
  }

  void test_findY() {
    auto ws = std::make_shared<WorkspaceTester>();
    ws->initialize(2, 2, 2);
    ws->mutableY(0) = 1.;
    ws->mutableY(1) = 2.;
    auto idx = ws->findY(0., {0, 0});
    TS_ASSERT_EQUALS(idx.first, -1);
    TS_ASSERT_EQUALS(idx.second, -1);
    idx = ws->findY(1., {0, 0});
    TS_ASSERT_EQUALS(idx.first, 0);
    TS_ASSERT_EQUALS(idx.second, 0);
    idx = ws->findY(1., {0, 1});
    TS_ASSERT_EQUALS(idx.first, 0);
    TS_ASSERT_EQUALS(idx.second, 1);
    idx = ws->findY(2., {1, 0});
    TS_ASSERT_EQUALS(idx.first, 1);
    TS_ASSERT_EQUALS(idx.second, 0);
    idx = ws->findY(2., {1, 1});
    TS_ASSERT_EQUALS(idx.first, 1);
    TS_ASSERT_EQUALS(idx.second, 1);
    ws->mutableY(1) = NAN;
    idx = ws->findY(NAN, {0, 0});
    TS_ASSERT_EQUALS(idx.first, 1);
    TS_ASSERT_EQUALS(idx.second, 0);
  }
  void testGetIntegratedSpectra() {
    WorkspaceTester workspace;
    workspace.initialize(5, 5, 4);
    MantidVec xValues = {1., 2., 3., 4., 5.};
    //set some values
    for (size_t wsIndex = 0; wsIndex < workspace.getNumberHistograms(); wsIndex++) {
      for (size_t binIndex = 0; binIndex < workspace.blocksize(); binIndex++) {
        // incrementing numbers
        double fillValue = static_cast<double>(binIndex + 1);
        // apart from some exceptions
        if (wsIndex == 1) {
          // all NaN
          fillValue = std::numeric_limits<double>::quiet_NaN();
        } else if (wsIndex == 2) {
          // all infinite
          fillValue = std::numeric_limits<double>::infinity();
        } else if ((wsIndex == 3) && (binIndex % 2 == 1)) {
          // alternate value NaN
          fillValue = std::numeric_limits<double>::quiet_NaN();
        } else if ((wsIndex == 4) && (binIndex % 2 == 0)) {
          // other alternate value inf
          fillValue = std::numeric_limits<double>::infinity();
        }
        workspace.mutableY(wsIndex)[binIndex] = fillValue;
      }
      //set the x values
      std::copy(xValues.begin(), xValues.end(), workspace.mutableX(wsIndex).begin());
    }
    MantidVec integratedValues;
    //the enitre range
    workspace.getIntegratedSpectra(integratedValues, 0, 0, true);
    MantidVec expected = {10., 0., 0., 4., 6.};
    TS_ASSERT_EQUALS(integratedValues, expected);
    // just the first two values
    workspace.getIntegratedSpectra(integratedValues, 0.0, 2.0, false);
    expected = {3., 0., 0., 1., 2.};
    TS_ASSERT_EQUALS(integratedValues, expected);
    // just the middle two values
    workspace.getIntegratedSpectra(integratedValues, 2.0, 3.9, false);
    expected = {5., 0., 0., 3., 2.};
    TS_ASSERT_EQUALS(integratedValues, expected);
    // just the last two values
    workspace.getIntegratedSpectra(integratedValues, 3.0, 5.0, false);
    expected = {7., 0., 0., 3., 4.};
    TS_ASSERT_EQUALS(integratedValues, expected);
  }

private:
  std::shared_ptr<WorkspaceTester>
  generateTestWorkspaceWithDistributionAndLabelSet(const bool distribution,
                                                   const std::string &yLabel) {
    auto testWS = std::make_shared<WorkspaceTester>();
    testWS->initialize(1, 2, 1);
    testWS->setDistribution(distribution);
    testWS->setYUnit("Counts");
    if (!yLabel.empty())
      testWS->setYUnitLabel(yLabel);
    testWS->getAxis(0)->setUnit("TOF");
    return testWS;
  }

  WorkspaceTester getWorkspaceWithPopulatedX(
      std::size_t const &nVectors, std::size_t const &xLength,
      std::size_t const &yLength, std::vector<double> const &xValues) {
    WorkspaceTester workspace;
    workspace.initialize(nVectors, xLength, yLength);
    auto &X = workspace.dataX(0);

    std::copy(xValues.begin(), xValues.end(), X.begin());
    return workspace;
  }

  void checkDetectorSignedTwoTheta(const Geometry::PointingAlong thetaSignAxis,
                                   const std::array<double, 4> &signs) {
    constexpr size_t numDets{4};
    constexpr size_t numBins{1};
    const auto frameUp = Geometry::Y;
    const auto frameAlongBeam = Geometry::Z;
    const auto frameSideways = Geometry::X;
    const auto frameThetaSign = thetaSignAxis;
    const auto frameHandedness = Geometry::Right;
    const std::string frameOrigin{"source"};
    auto refFrame = std::make_shared<ReferenceFrame>(
        frameUp, frameAlongBeam, frameThetaSign, frameHandedness, frameOrigin);
    std::shared_ptr<MatrixWorkspace> ws = std::make_shared<WorkspaceTester>();
    ws->initialize(numDets, numBins, numBins);
    // Create instrument with four detectors to play with.
    auto instrument = std::make_shared<Instrument>("TestInstrument");
    instrument->setReferenceFrame(refFrame);
    constexpr double twoTheta{4.2 / 180. * M_PI};
    for (size_t i = 0; i < numDets; ++i) {
      Detector *det =
          new Detector("pixel", static_cast<detid_t>(i), instrument.get());
      constexpr double r{1.};
      const double rotation =
          (45. + 90. * static_cast<double>(i)) / 180. * M_PI;
      const double x = r * std::sin(twoTheta) * std::cos(rotation);
      const double y = r * std::sin(twoTheta) * std::sin(rotation);
      const double z = r * std::cos(twoTheta);
      V3D pos;
      pos[frameUp] = y;
      pos[frameAlongBeam] = z;
      pos[frameSideways] = x;
      det->setShape(ComponentCreationHelper::createSphere(0.01, pos, "1"));
      det->setPos(pos);
      instrument->add(det);
      instrument->markAsDetector(det);
      ws->getSpectrum(i).addDetectorID(static_cast<detid_t>(i));
    }
    V3D pos(0., 0., 0.);
    ComponentCreationHelper::addSampleToInstrument(instrument, pos);
    pos[frameAlongBeam] = -1.;
    ComponentCreationHelper::addSourceToInstrument(instrument, pos);
    ws->setInstrument(instrument);
    for (detid_t detid = 0; static_cast<size_t>(detid) < numDets; ++detid) {
      auto det = instrument->getDetector(detid);
      const auto signedTwoTheta = ws->detectorSignedTwoTheta(*det);
      TS_ASSERT_DELTA(signedTwoTheta, signs[detid] * twoTheta, 1e-12)
    }
  }

  Mantid::API::MantidImage_sptr createImage(const size_t width,
                                            const size_t height) {
    auto image =
        std::make_shared<Mantid::API::MantidImage>(height, MantidVec(width));
    double startingValue = 1.0;
    for (auto &row : *image) {
      std::iota(row.begin(), row.end(), startingValue);
      startingValue += static_cast<double>(width);
    }
    return image;
  }

  /**
   * Create a test workspace. Can be histo or points depending on x/yLength.
   * @param nVectors :: [input] Number of vectors
   * @param xLength :: [input] Length of X vector
   * @param yLength :: [input] Length of Y, E vectors
   * @returns :: workspace
   */
  WorkspaceTester createTestWorkspace(size_t nVectors, size_t xLength,
                                      size_t yLength) {
    WorkspaceTester ws;
    ws.initialize(nVectors, xLength, yLength);
    // X data
    std::vector<double> xData(xLength);
    std::iota(xData.begin(), xData.end(), 0.0);

    // Y data
    const auto yCounts = [&yLength](size_t wi) {
      std::vector<double> v(yLength);
      std::iota(v.begin(), v.end(), static_cast<double>(wi) * 10.0);
      return v;
    };

    // E data
    const auto errors = [&yLength](size_t wi) {
      std::vector<double> v(yLength);
      std::generate(v.begin(), v.end(), [&wi]() {
        return std::sqrt(static_cast<double>(wi) * 10.0);
      });
      return v;
    };

    for (size_t wi = 0; wi < nVectors; ++wi) {
      if (xLength == yLength) {
        ws.setPoints(wi, xData);
      } else if (xLength == yLength + 1) {
        ws.setBinEdges(wi, xData);
      } else {
        throw std::invalid_argument(
            "yLength must either be equal to xLength or xLength - 1");
      }
      ws.setCounts(wi, yCounts(wi));
      ws.setCountStandardDeviations(wi, errors(wi));
    }
    return ws;
  }

  std::shared_ptr<MatrixWorkspace> ws;
};

class MatrixWorkspaceTestPerformance : public CxxTest::TestSuite {

public:
  static MatrixWorkspaceTestPerformance *createSuite() {
    return new MatrixWorkspaceTestPerformance();
  }
  static void destroySuite(MatrixWorkspaceTestPerformance *suite) {
    delete suite;
  }

  MatrixWorkspaceTestPerformance() : m_workspace() {
    using namespace Mantid::Geometry;

    size_t numberOfHistograms = 10000;
    size_t numberOfBins = 1;
    m_workspace.initialize(numberOfHistograms, numberOfBins + 1, numberOfBins);
    bool includeMonitors = false;
    bool startYNegative = true;
    const std::string instrumentName("SimpleFakeInstrument");
    InstrumentCreationHelper::addFullInstrumentToWorkspace(
        m_workspace, includeMonitors, startYNegative, instrumentName);

    Mantid::Kernel::V3D sourcePos(0, 0, 0);
    Mantid::Kernel::V3D samplePos(0, 0, 1);
    Mantid::Kernel::V3D trolley1Pos(0, 0, 3);
    Mantid::Kernel::V3D trolley2Pos(0, 0, 6);
    m_paramMap = std::make_shared<Mantid::Geometry::ParameterMap>();

    auto baseInstrument = ComponentCreationHelper::sansInstrument(
        sourcePos, samplePos, trolley1Pos, trolley2Pos);

    auto sansInstrument =
        std::make_shared<Instrument>(baseInstrument, m_paramMap);

    // See component creation helper for instrument definition
    m_sansBank = sansInstrument->getComponentByName("Bank1");

    numberOfHistograms = sansInstrument->getNumberDetectors();
    m_workspaceSans.initialize(numberOfHistograms, numberOfBins + 1,
                               numberOfBins);
    m_workspaceSans.setInstrument(sansInstrument);
    m_workspaceSans.getAxis(0)->setUnit("TOF");
    m_workspaceSans.rebuildSpectraMapping();

    m_zRotation =
        Mantid::Kernel::Quat(180, V3D(0, 0, 1)); // rotate 180 degrees around z

    m_pos = Mantid::Kernel::V3D(1, 1, 1);
  }
  /// This test is equivalent to GeometryInfoFactoryTestPerformance, see there.
  void test_typical() {
    auto instrument = m_workspace.getInstrument();
    auto source = instrument->getSource();
    auto sample = instrument->getSample();
    auto L1 = source->getDistance(*sample);
    double result = 0.0;
    for (size_t i = 0; i < 10000; ++i) {
      auto detector = m_workspace.getDetector(i);
      result += L1;
      result += detector->getDistance(*sample);
      result += m_workspace.detectorTwoTheta(*detector);
    }
    // We are computing and using the result to fool the optimizer.
    TS_ASSERT_DELTA(result, 5214709.740869, 1e-6);
  }

  void test_calculateL2() {

    /*
     * Simulate the L2 calculation performed via the Workspace/Instrument
     * interface.
     */
    auto instrument = m_workspaceSans.getInstrument();
    auto sample = instrument->getSample();
    double l2 = 0;
    for (size_t i = 0; i < m_workspaceSans.getNumberHistograms(); ++i) {
      auto detector = m_workspaceSans.getDetector(i);
      l2 += detector->getDistance(*sample);
    }
    // Prevent optimization
    TS_ASSERT(l2 > 0);
  }

  void test_calculateL2_x10() {

    /*
     * Simulate the L2 calculation performed via the Workspace/Instrument
     * interface. Repeat several times to benchmark any caching/optmisation that
     * might be taken place in parameter maps.
     */
    auto instrument = m_workspaceSans.getInstrument();
    auto sample = instrument->getSample();
    double l2 = 0;
    int count = 0;
    while (count < 10) {
      for (size_t i = 0; i < m_workspaceSans.getNumberHistograms(); ++i) {
        auto detector = m_workspaceSans.getDetector(i);
        l2 += detector->getDistance(*sample);
      }
      ++count;
    }
    // Prevent optimization
    TS_ASSERT(l2 > 0);
  }

  /*
   * Rotate a bank in the workspace and read the positions out again. Very
   * typical.
   */
  void test_rotate_bank_and_read_positions_x10() {

    using namespace Mantid::Geometry;
    using namespace Mantid::Kernel;

    int count = 0;
    // Repeated execution to improve statistics and for comparison purposes with
    // future updates
    while (count < 10) {
      // Rotate the bank
      auto &compInfo = m_workspaceSans.mutableComponentInfo();
      compInfo.setRotation(compInfo.indexOf(m_sansBank->getComponentID()),
                           m_zRotation);

      V3D pos;
      for (size_t i = 1; i < m_workspaceSans.getNumberHistograms(); ++i) {
        pos += m_workspaceSans.getDetector(i)->getPos();
      }
      ++count;
    }
  }

  /*
   * Move a bank in the workspace and read the positions out again. Very
   * typical.
   */
  void test_move_bank_and_read_positions_x10() {

    using namespace Mantid::Geometry;
    using namespace Mantid::Kernel;

    int count = 0;
    // Repeated execution to improve statistics and for comparison purposes with
    // future updates
    while (count < 10) {
      // move the bank
      auto &compInfo = m_workspaceSans.mutableComponentInfo();
      compInfo.setPosition(compInfo.indexOf(m_sansBank->getComponentID()),
                           m_pos);

      V3D pos;
      for (size_t i = 1; i < m_workspaceSans.getNumberHistograms(); ++i) {
        pos += m_workspaceSans.getDetector(i)->getPos();
      }
      ++count;
    }
  }

  // As test_rotate_bank_and_read_positions_x10 but based on SpectrumInfo.
  void test_rotate_bank_and_read_positions_SpectrumInfo_x10() {
    int count = 0;
    while (count < 10) {
      // Rotate the bank
      auto &compInfo = m_workspaceSans.mutableComponentInfo();
      compInfo.setRotation(compInfo.indexOf(m_sansBank->getComponentID()),
                           m_zRotation);

      V3D pos;
      const auto &spectrumInfo = m_workspaceSans.spectrumInfo();
      for (size_t i = 1; i < m_workspaceSans.getNumberHistograms(); ++i) {
        pos += spectrumInfo.position(i);
      }
      ++count;
    }
  }

  // As test_move_bank_and_read_positions_x10 but based on SpectrumInfo.
  void test_move_bank_and_read_positions_SpectrumInfo_x10() {
    int count = 0;
    while (count < 10) {
      // move the bank
      auto &compInfo = m_workspaceSans.mutableComponentInfo();
      compInfo.setPosition(compInfo.indexOf(m_sansBank->getComponentID()),
                           m_pos);

      V3D pos;
      const auto &spectrumInfo = m_workspaceSans.spectrumInfo();
      for (size_t i = 1; i < m_workspaceSans.getNumberHistograms(); ++i) {
        pos += spectrumInfo.position(i);
      }
      ++count;
    }
  }

  void test_hasOrientedLattice() {
    // create a workspace without an oriented lattice (or sample)
    std::shared_ptr<MatrixWorkspace> ws(makeWorkspaceWithDetectors(3, 1));
    TSM_ASSERT_EQUALS(
        "A newly created workspace should not have an oriented lattice",
        ws->hasOrientedLattice(), false);

    // add an oriented lattice
    ws->mutableSample().setOrientedLattice(
        std::make_unique<OrientedLattice>(1.0, 2.0, 3.0, 90, 90, 90));
    TSM_ASSERT_EQUALS("A workspace with an oriented lattice should report true",
                      ws->hasOrientedLattice(), true);

    // remove it again
    ws->mutableSample().clearOrientedLattice();
    TSM_ASSERT_EQUALS(
        "workspace with it's oriented lattice cleared should report false",
        ws->hasOrientedLattice(), false);
  }

  void test_isGroup() {
    std::shared_ptr<MatrixWorkspace> ws(makeWorkspaceWithDetectors(3, 1));
    TS_ASSERT_EQUALS(ws->isGroup(), false);
  }

private:
  WorkspaceTester m_workspace;
  WorkspaceTester m_workspaceSans;
  Mantid::Kernel::Quat m_zRotation;
  Mantid::Kernel::V3D m_pos;
  Mantid::Geometry::IComponent_const_sptr m_sansBank;
  std::shared_ptr<Mantid::Geometry::ParameterMap> m_paramMap;
};
