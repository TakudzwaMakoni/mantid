#ifndef MANTID_MDALGORITHMS_ACCUMULATEMDTEST_H_
#define MANTID_MDALGORITHMS_ACCUMULATEMDTEST_H_

#include <cxxtest/TestSuite.h>
#include "MantidMDAlgorithms/AccumulateMD.h"
#include "MantidTestHelpers/WorkspaceCreationHelper.h"
#include "MantidTestHelpers/MDEventsTestHelper.h"
#include "MantidKernel/ConfigService.h"
#include "MantidAPI/AlgorithmManager.h"
#include <Poco/Path.h>
#include <Poco/File.h>

using Mantid::MDAlgorithms::AccumulateMD;
using namespace Mantid::API;
using namespace Mantid::DataObjects;

class AccumulateMDTest : public CxxTest::TestSuite {
public:
  // This pair of boilerplate methods prevent the suite being created statically
  // This means the constructor isn't called when running other tests
  static AccumulateMDTest *createSuite() { return new AccumulateMDTest(); }

  static void destroySuite(AccumulateMDTest *suite) { delete suite; }

  void test_Init() {
    AccumulateMD alg;
    TS_ASSERT_THROWS_NOTHING(alg.initialize());
  }

  void test_pad_parameter_vector_empty() {
    std::vector<double> test_param_vector;
    unsigned long grow_to = 8;
    Mantid::MDAlgorithms::padParameterVector(test_param_vector, grow_to);

    TS_ASSERT_EQUALS(test_param_vector.size(), 8);
    TS_ASSERT_EQUALS(test_param_vector[4], 0.0);
  }

  void test_pad_parameter_vector_values() {
    std::vector<double> test_param_vector(1, 3.7);
    unsigned long grow_to = 8;
    Mantid::MDAlgorithms::padParameterVector(test_param_vector, grow_to);

    TS_ASSERT_EQUALS(test_param_vector.size(), 8);
    TS_ASSERT_EQUALS(test_param_vector[4], 3.7);
  }

  void test_filter_to_existing_sources_file_nonexist() {
    // Create vector of data_sources to filter
    std::vector<std::string> data_sources;

    // Create vector for other parameters
    std::vector<double> psi(1, 0.0);
    std::vector<double> gl(1, 0.0);
    std::vector<double> gs(1, 0.0);
    std::vector<double> efix(1, 0.0);

    // Add absolute path to a file which doesn't exist
    Poco::Path filepath =
        Poco::Path(Mantid::Kernel::ConfigService::Instance().getTempDir(),
                   "ACCUMULATEMDTEST_NONEXISTENTFILE");
    data_sources.push_back(filepath.toString());

    Mantid::MDAlgorithms::filterToExistingSources(data_sources, psi, gl, gs,
                                                  efix);

    TS_ASSERT(data_sources.empty());
  }

  void test_filter_to_existing_sources_workspace_nonexist() {
    // Create vector of data_sources to filter
    std::vector<std::string> data_sources;

    // Create vector for other parameters
    std::vector<double> psi(1, 0.0);
    std::vector<double> gl(1, 0.0);
    std::vector<double> gs(1, 0.0);
    std::vector<double> efix(1, 0.0);

    data_sources.push_back("ACCUMULATEMDTEST_NONEXISTENTWORKSPACE");

    Mantid::MDAlgorithms::filterToExistingSources(data_sources, psi, gl, gs,
                                                  efix);

    TS_ASSERT(data_sources.empty());
  }

  void test_filter_to_existing_sources_workspace_exist() {
    // Create vector of data_sources to filter
    std::vector<std::string> data_sources;

    // Create vector for other parameters
    std::vector<double> psi(1, 0.0);
    std::vector<double> gl(1, 0.0);
    std::vector<double> gs(1, 0.0);
    std::vector<double> efix(1, 0.0);

    // Create a cheap workspace
    std::string ws_name = "ACCUMULATEMDTEST_EXISTENTWORKSPACE";
    auto bkg_ws = WorkspaceCreationHelper::Create1DWorkspaceRand(1);
    // add to ADS (no choice but to use ADS here)
    AnalysisDataService::Instance().add(ws_name, bkg_ws);

    data_sources.push_back(ws_name);

    Mantid::MDAlgorithms::filterToExistingSources(data_sources, psi, gl, gs,
                                                  efix);

    TS_ASSERT(!data_sources.empty());

    // Remove workspace from the data service.
    AnalysisDataService::Instance().remove(ws_name);
  }

  void test_filter_to_existing_sources_file_exist() {
    // Create vector of data_sources to filter
    std::vector<std::string> data_sources;

    // Create vector for other parameters
    std::vector<double> psi(1, 0.0);
    std::vector<double> gl(1, 0.0);
    std::vector<double> gs(1, 0.0);
    std::vector<double> efix(1, 0.0);

    // Create a temporary file to find
    Poco::Path filepath =
        Poco::Path(Mantid::Kernel::ConfigService::Instance().getTempDir(),
                   "ACCUMULATEMDTEST_EXISTENTFILE");
    Poco::File existent_file(filepath);
    existent_file.createFile();
    data_sources.push_back(filepath.toString());

    Mantid::MDAlgorithms::filterToExistingSources(data_sources, psi, gl, gs,
                                                  efix);

    TS_ASSERT(!data_sources.empty());

    // Remove the temporary file
    existent_file.remove();
  }

  void test_filter_to_new_none_new() {
    std::vector<std::string> input_data;
    input_data.push_back("test1");
    input_data.push_back("test2");
    input_data.push_back("test3");
    std::vector<std::string> current_data = input_data;

    // Create vector for other parameters
    std::vector<double> psi(3, 0.0);
    std::vector<double> gl(3, 0.0);
    std::vector<double> gs(3, 0.0);
    std::vector<double> efix(3, 0.0);

    Mantid::MDAlgorithms::filterToNew(input_data, current_data, psi, gl, gs,
                                      efix);

    // Two input vectors were identical, so we should get an empty vector back
    TS_ASSERT(input_data.empty());

    // Parameter vectors should also have been emptied
    TS_ASSERT(psi.empty());
    TS_ASSERT(gl.empty());
    TS_ASSERT(gs.empty());
    TS_ASSERT(efix.empty());
  }

  void test_filter_to_new() {
    std::vector<std::string> input_data;
    input_data.push_back("test1");
    input_data.push_back("test2");
    input_data.push_back("test3");
    input_data.push_back("test4");
    input_data.push_back("test5");

    std::vector<std::string> current_data;
    current_data.push_back("test1");
    current_data.push_back("test3");
    current_data.push_back("test4");

    // Create vector for other parameters
    std::vector<double> psi(5, 0.0);
    std::vector<double> gl(5, 0.0);
    std::vector<double> gs(5, 0.0);
    std::vector<double> efix(5, 0.0);

    Mantid::MDAlgorithms::filterToNew(input_data, current_data, psi, gl, gs,
                                      efix);

    // test2 and test5 is new data (it is in input_data but not current_data)
    // and so should be returned in the vector
    TS_ASSERT_EQUALS(input_data[0], "test2");
    TS_ASSERT_EQUALS(input_data[1], "test5");

    // Parameter vectors should have been reduced to the same size
    TS_ASSERT_EQUALS(psi.size(), input_data.size());
    TS_ASSERT_EQUALS(gl.size(), input_data.size());
    TS_ASSERT_EQUALS(gs.size(), input_data.size());
    TS_ASSERT_EQUALS(efix.size(), input_data.size());
  }

  void test_insert_data_sources() {
    std::string data_sources = "test1,test2,test3";
    std::set<std::string> data_sources_set;
    Mantid::MDAlgorithms::insertDataSources(data_sources, data_sources_set);

    // Check set contains "test1", "test2" and "test3"
    std::set<std::string>::iterator iter;
    iter = data_sources_set.find("test1");
    TS_ASSERT(iter != data_sources_set.end());

    iter = data_sources_set.find("test2");
    TS_ASSERT(iter != data_sources_set.end());

    iter = data_sources_set.find("test3");
    TS_ASSERT(iter != data_sources_set.end());
  }

  void test_insert_data_sources_with_whitespace() {
    std::string data_sources = " test1,test2 , test3";
    std::set<std::string> data_sources_set;
    Mantid::MDAlgorithms::insertDataSources(data_sources, data_sources_set);

    // Check set contains "test1", "test2" and "test3"
    std::set<std::string>::iterator iter;
    iter = data_sources_set.find("test1");
    TS_ASSERT(iter != data_sources_set.end());

    iter = data_sources_set.find("test2");
    TS_ASSERT(iter != data_sources_set.end());

    iter = data_sources_set.find("test3");
    TS_ASSERT(iter != data_sources_set.end());
  }

  void test_algorithm_append_data() {

    const std::string out_ws_name = "mdew_output";
    const std::string in_ws_name = "mdew_input";
    const std::string sample_ws_name = "sample_data_1";

    auto alg = Mantid::API::AlgorithmManager::Instance().create(
        "CreateSimulationWorkspace");
    alg->initialize();
    alg->setPropertyValue("Instrument", "MAR");
    alg->setPropertyValue("BinParams", "-3,1,3");
    alg->setPropertyValue("UnitX", "DeltaE");
    alg->setProperty("OutputWorkspace", sample_ws_name);
    alg->execute();

    alg = Mantid::API::AlgorithmManager::Instance().create("AddSampleLog");
    alg->initialize();
    alg->setPropertyValue("Workspace", sample_ws_name);
    alg->setPropertyValue("LogName", "Ei");
    alg->setPropertyValue("LogText", "3.0");
    alg->setPropertyValue("LogType", "Number");
    alg->execute();

    Mantid::API::IMDEventWorkspace_sptr in_ws =
        MDEventsTestHelper::makeFakeMDEventWorkspace(in_ws_name, 500);

    AccumulateMD acc_alg;
    acc_alg.initialize();
    acc_alg.setPropertyValue("InputWorkspace", in_ws_name);
    acc_alg.setPropertyValue("OutputWorkspace", out_ws_name);
    acc_alg.setPropertyValue("DataSources", sample_ws_name);
    acc_alg.setPropertyValue("Alatt", "1.4165,1.4165,1.4165");
    acc_alg.setPropertyValue("Angdeg", "90,90,90");
    acc_alg.setPropertyValue("u", "1,0,0");
    acc_alg.setPropertyValue("v", "0,1,0");
    TS_ASSERT_THROWS_NOTHING(acc_alg.execute());

    // TODO Get output workspace and check it has the sum of the number of
    // events in sample_data_1 and mdew_input
    auto sample_ws = alg->getProperty("Workspace");
    auto out_ws = acc_alg.getProperty("OutputWorkspace");

    // Clean up
    // Remove workspaces from the data service.
    AnalysisDataService::Instance().clear();
  }

  void test_algorithm_clean_option() {
    // Same as previous test but with Clean option
  }
};

#endif /* MANTID_MDALGORITHMS_ACCUMULATEMDTEST_H_ */
