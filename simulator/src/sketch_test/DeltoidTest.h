/**
 * @file DeltoidTest.h
 * @author deadlycat (lsmfttb@gmail.com)
 * @brief Test Deltoid
 *
 * @copyright Copyright (c) 2022
 *
 */
#pragma once

#include <common/test.h>
#include <sketch/Deltoid.h>

#define DT_PARA_PATH "DT.para"
#define DT_TEST_PATH "DT.test"
#define DT_DATA_PATH "DT.data"

namespace OmniSketch::Test {

/**
 * @brief Testing class for Deltoid
 *
 */
template <int32_t key_len, typename T, typename hash_t = Hash::AwareHash>
class DeltoidTest : public TestBase<key_len, T> {
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
  DeltoidTest(const std::string_view config_file)
      : TestBase<key_len, T>("Deltoid", config_file, DT_TEST_PATH) {}

  /**
   * @brief Test Deltoid
   * @details An overriden method
   */
  void runTest() override;
};

} // namespace OmniSketch::Test

//-----------------------------------------------------------------------------
//
///                        Implementation of template methods
//
//-----------------------------------------------------------------------------

namespace OmniSketch::Test {

template <int32_t key_len, typename T, typename hash_t>
void DeltoidTest<key_len, T, hash_t>::runTest() {
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
  int32_t num_hash, num_group; // sketch config
  double num_heavy_hitter;
  std::string data_file; // data config
  toml::array arr;       // shortly we will convert it to format
  /// Step ii. Open the config file
  Util::ConfigParser parser(config_file);
  if (!parser.succeed()) {
    return;
  }
  /// Step iii. Set the working node of the parser. In this case, [Deltoid].
  parser.setWorkingNode(
      DT_PARA_PATH); // do not forget to to enclose it with braces
  /// Step iv. Parse num_hash, num_group and threshold
  if (!parser.parseConfig(num_hash, "num_hash"))
    return;
  if (!parser.parseConfig(num_group, "num_group"))
    return;
  /// Step v. To know about the data, we move to the [Deltoid.data] node.
  parser.setWorkingNode(DT_DATA_PATH);
  /// Step vi. Parse data and format
  if (!parser.parseConfig(num_heavy_hitter, "threshold_heavy_hitter"))
    return;
  if (!parser.parseConfig(data_file, "data"))
    return;
  if (!parser.parseConfig(arr, "format"))
    return;
  Data::DataFormat format(arr); // conver from toml::array to Data::DataFormat
  /// [Optional] User-defined rules
  std::string method;
  Data::HXMethod hx_method = Data::TopK;
  if (!parser.parseConfig(method, "hx_method"))
    return;
  if (!method.compare("Percentile")) {
    hx_method = Data::Percentile;
  }
  Data::CntMethod cnt_method = Data::InLength;
  if (!parser.parseConfig(method, "cnt_method"))
    return;
  if (!method.compare("InPacket")) {
    cnt_method = Data::InPacket;
  }

  /// Step ii. Get ground truth
  ///
  ///       1. read data
  StreamData data(data_file, format); // specify both data file and data format
  if (!data.succeed())
    return;
  Data::GndTruth<key_len, T> gnd_truth, gnd_truth_heavy_hitters;
  gnd_truth.getGroundTruth(data.begin(), data.end(), cnt_method);
  gnd_truth_heavy_hitters.getHeavyHitter(gnd_truth, num_heavy_hitter,
                                         hx_method);
  ///       2. [optional] show data info
  fmt::print("DataSet: {:d} records with {:d} keys ({})\n", data.size(),
             gnd_truth.size(), data_file);

  OmniSketch::Hash::AwareHash(1);

  /// Part II.
  ///   Prepare sketch and data
  ///
  /// Step i. Initialize a sketch
  std::unique_ptr<Sketch::SketchBase<key_len, T>> ptr(
      new Sketch::Deltoid<key_len, T, hash_t>(num_hash, num_group));
  /// remember that the left ptr must point to the base class in order to call
  /// the methods in it
  /// Step iii. Insert the samples and then look up all the flows
  ///
  ///        1. update records into the sketch
  this->testUpdate(ptr, data.begin(), data.end(),
                   cnt_method); // metrics of interest are in config file
  ///        2. query for all the flowkeys
  this->testQuery(ptr, gnd_truth); // metrics of interest are in config file
  ///        3. query for heavy hitters
  if (hx_method == Data::TopK) {
    this->testHeavyHitter(
        ptr, gnd_truth_heavy_hitters.min(),
        gnd_truth_heavy_hitters); // metrics of interest are in config file
  } else {
    this->testHeavyHitter(
        ptr, std::floor(gnd_truth.totalValue() * num_heavy_hitter + 1),
        gnd_truth_heavy_hitters); // gnd_truth_heavy_hitter: >, yet HashPipe: >=
  }
  ///        4. size
  this->testSize(ptr);
  ///        5. show metrics
  this->show();

  return;
}

} // namespace OmniSketch::Test

#undef DT_PARA_PATH
#undef DT_TEST_PATH
#undef DT_DATA_PATH

// Driver instance:
//      AUTHOR:
//      CONFIG: sketch_config.toml  # with respect to the `src/` directory
//    TEMPLATE: <13, int32_t, Hash::AwareHash>