/**
 * @file PRSketchTest.h
 * @author XierLabber (you@domain.com)
 * @brief Test PR Sketch
 *
 * @copyright Copyright (c) 2022
 *
 */
#pragma once

#include <common/test.h>
#include <sketch/PRSketch.h>

#define PR_PARA_PATH "PR.para"
#define PR_TEST_PATH "PR.test"
#define PR_DATA_PATH "PR.data"

namespace OmniSketch::Test {

/**
 * @brief Testing class for PR Sketch
 *
 */
template <int32_t key_len, typename T, typename hash_t = Hash::AwareHash>
class PRSketchTest : public TestBase<key_len, T> {
  using TestBase<key_len, T>::config_file;

public:
  /**
   * @brief Constructor
   * @details Names from left to right are
   * - show name
   * - config file
   * - path to the node that contains metrics of interest (concatenated with
   * '.')
   */
  PRSketchTest(const std::string_view config_file)
      : TestBase<key_len, T>("PR Sketch", config_file, PR_TEST_PATH) {}

  /**
   * @brief Test Bloom Filter
   * @details An overriden method
   */
  void runTest() override;
};

} // namespace OmniSketch::Test

//-----------------------------------------------------------------------------
//
///                        Implementation of templated methods
//
//-----------------------------------------------------------------------------

namespace OmniSketch::Test {

template <int32_t key_len, typename T, typename hash_t>
void PRSketchTest<key_len, T, hash_t>::runTest() {
  /**
   * @brief shorthand for convenience
   *
   */
  using StreamData = Data::StreamData<key_len>;

  /// Part I.
  ///   Parse the config file
  ///
  /// Step i.  First we list the variables to parse, namely:
  ///
  int32_t counter_length, counter_hash_num,
          filter_length, filter_hash_num;   
  T phi;                                       // sketch config
  std::string data_file; // data config
  toml::array arr;       // shortly we will convert it to format
  /// Step ii. Open the config file
  Util::ConfigParser parser(config_file);
  if (!parser.succeed()) {
    return;
  }
  /// Step iii. Set the working node of the parser.
  parser.setWorkingNode(
      PR_PARA_PATH); // do not forget to to enclose it with braces
  /// Step iv. Parse num_bits and num_hash
  if (!parser.parseConfig(counter_length, "counter_length"))
    return;
  if (!parser.parseConfig(counter_hash_num, "counter_hash_num"))
    return;
  if (!parser.parseConfig(filter_length, "filter_length"))
    return;
  if (!parser.parseConfig(filter_hash_num, "filter_hash_num"))
    return;
  if (!parser.parseConfig(phi, "phi"))
    return;
  /// Step v. Move to the data node
  parser.setWorkingNode(PR_DATA_PATH);
  /// Step vi. Parse data and format
  if (!parser.parseConfig(data_file, "data"))
    return;
  if (!parser.parseConfig(arr, "format"))
    return;
  Data::DataFormat format(arr); // conver from toml::array to Data::DataFormat
  /// [Optional] User-defined rules
  ///
  /// Step vii. Parse Cnt Method.
  std::string method;
  Data::CntMethod cnt_method = Data::InLength;
  if (!parser.parseConfig(method, "cnt_method"))
    return;
  if (!method.compare("InPacket")) {
    cnt_method = Data::InPacket;
  }

  /// Part II.
  ///   Prepare sketch and data
  ///
  /// Step i. Initialize a sketch
  std::unique_ptr<Sketch::SketchBase<key_len, T>> ptr(
      new Sketch::PRSketch<key_len, T, hash_t>(counter_length, counter_hash_num, 
                                               filter_length, filter_hash_num, phi));
  /// remember that the left ptr must point to the base class in order to call
  /// the methods in it

  /// Step ii. Get ground truth
  ///
  ///       1. read data
  StreamData data(data_file, format); // specify both data file and data format
  if (!data.succeed())
    return;
  Data::GndTruth<key_len, T> gnd_truth;
  gnd_truth.getGroundTruth(data.begin(), data.end(), cnt_method);
  ///       2. [optional] show data info
  fmt::print("DataSet: {:d} records with {:d} keys ({})\n", data.size(),
             gnd_truth.size(), data_file);
  /// Step iii. Insert the samples and then look up all the flows
  ///
  ///        1. update records into the sketch
  this->testUpdate(ptr, data.begin(), data.end(),
                   cnt_method); // metrics of interest are in config file
  ///        2. query for all the flowkeys
  this->testQuery(ptr, gnd_truth); // metrics of interest are in config file
  ///        3. size
  this->testSize(ptr);
  ///        3. show metrics
  this->show();

  return;
}

} // namespace OmniSketch::Test

#undef PR_PARA_PATH
#undef PR_TEST_PATH
#undef PR_DATA_PATH

// Driver instance:
//      AUTHOR: dromniscience
//      CONFIG: sketch_config.toml  # with respect to the `src/` directory
//    TEMPLATE: <13, int32_t, Hash::AwareHash>