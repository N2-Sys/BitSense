
/**
 * @file NZESketchTest.h
 * @author XierLabber (you@domain.com)
 * @brief Test NZE
 *
 * @copyright Copyright (c) 2022
 *
 */
#pragma once

#include <common/test.h>
#include <sketch/NZESketch.h>

#define NZE_PARA_PATH "NZE.para"
#define NZE_TEST_PATH "NZE.test"
#define NZE_DATA_PATH "NZE.data"

namespace OmniSketch::Test {

/**
 * @brief Testing class for NZESketch
 *
 */
template <int32_t key_len, typename T, typename hash_t = Hash::AwareHash>
class NZESketchTest : public TestBase<key_len, T> {
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
  NZESketchTest(const std::string_view config_file)
      : TestBase<key_len, T>("NZE Sketch", config_file, NZE_TEST_PATH) {}

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
void NZESketchTest<key_len, T, hash_t>::runTest() {
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
  int32_t HTLength, BFBitsNum, BFHashNum, FSdepth, FSwidth;  // sketch config
  std::string data_file; // data config
  toml::array arr;       // shortly we will convert it to format
  /// Step ii. Open the config file
  Util::ConfigParser parser(config_file);
  if (!parser.succeed()) {
    return;
  }
  /// Step iii. Set the working node of the parser.
  parser.setWorkingNode(
      NZE_PARA_PATH); // do not forget to to enclose it with braces
  /// Step iv. Parse num_bits and num_hash
  if (!parser.parseConfig(HTLength, "HTLength"))
    return;
  if (!parser.parseConfig(BFBitsNum, "BFBitsNum"))
    return;
  if (!parser.parseConfig(BFHashNum, "BFHashNum"))
    return;
  if (!parser.parseConfig(FSdepth, "FSdepth"))
    return;
  if (!parser.parseConfig(FSwidth, "FSwidth"))
    return;
  /// Step v. Move to the data node
  parser.setWorkingNode(NZE_DATA_PATH);
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
      new Sketch::NZESketch<key_len, T, hash_t>(HTLength, BFBitsNum, BFHashNum, FSdepth, FSwidth));
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
  double result = this->testQuery(ptr, gnd_truth); // metrics of interest are in config file
  ///        3. size
  this->testSize(ptr);
  ///        3. show metrics
  this->show();

  printf("RESULT: %lf\n",result);

  return;
}

} // namespace OmniSketch::Test

#undef NZE_PARA_PATH
#undef NZE_TEST_PATH
#undef NZE_DATA_PATH

// Driver instance:
//      AUTHOR: XierLabber
//      CONFIG: config/sketch_config.toml  # with respect to the `simulator/` directory
//    TEMPLATE: <13, int32_t, Hash::AwareHash>
