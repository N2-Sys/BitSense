/**
 * @file HeavyKeeper.h
 * @author XierLabber (you@domain.com)
 * @brief Heavy Keeper
 *
 * @copyright Copyright (c) 2022
 *
 */
#pragma once

#include <common/test.h>
#include <sketch/HeavyKeeper.h>

#define HK_PARA_PATH "HK.para"
#define HK_TEST_PATH "HK.test"
#define HK_DATA_PATH "HK.data"

namespace OmniSketch::Test {
/**
 * @brief Testing class for Heavy Keeper
 *
 */
 template <int32_t key_len, typename T, typename hash_t = Hash::AwareHash>
  class HeavyKeeperTest : public TestBase<key_len, T> {
    using TestBase<key_len, T>::config_file;

  public:
  /**
   * @brief Constructor
   * @details Names from left to right are
   * - name for showing
   * - config file [This one should be passed from the .cpp file]
   * - table name in the config file
   * - sub-table name (also a node name) under the table, containing metrics
   * of interest
   * - sub-table name (also a node name) under the table, containing data info
   */
  HeavyKeeperTest(const std::string_view config_file)
      : TestBase<key_len, T>("Heavy Keeper", config_file, HK_TEST_PATH) {}

  /**
   * @brief Test Heavy Keeper
   * @details An overriden method
   */
  void runTest() override;

 };

}// namespace OmniSketch::Test

//-----------------------------------------------------------------------------
//
///                        Implementation of template methods
//
//-----------------------------------------------------------------------------

namespace OmniSketch::Test {

 template <int32_t key_len, typename T, typename hash_t>
 void HeavyKeeperTest<key_len, T, hash_t> ::runTest() {
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
  int32_t depth_, width_, num_threshold_;
  double b_, hash_table_alpha;
  // sketch config
  double num_heavy_hitter_;  
  std::string data_file; // data config
  toml::array arr;       // shortly we will convert it to format
  /// Step ii. Open the config file
  Util::ConfigParser parser(config_file);
  if (!parser.succeed()) {
    return;
  }
  /// Step iii. Set the working node of the parser. In this case, [HK].
  parser.setWorkingNode(
      HK_PARA_PATH); // do not forget to to enclose it with braces
  /// Step iv. Parse num_bits and num_hash
  if (!parser.parseConfig(depth_, "depth"))
    return;
  if (!parser.parseConfig(width_, "width"))
    return;
  if (!parser.parseConfig(num_threshold_, "num_threshold"))
    return;
  if (!parser.parseConfig(b_, "b"))
    return;
  if (!parser.parseConfig(hash_table_alpha, "hash_table_alpha"))
    return;
  /// Step v. To know about the data, we move to the [HK.data] node.
  parser.setWorkingNode(HK_DATA_PATH);
  /// Step vi. Parse data and format
  if (!parser.parseConfig(num_heavy_hitter_, "threshold_heavy_hitter"))
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

  /// Part II.
  ///   Prepare sketch and data
  ///
  /// Step i. Initialize a sketch
  std::unique_ptr<Sketch::SketchBase<key_len, T>> ptr(
      new Sketch::HeavyKeeper<key_len, T, hash_t>
                             (depth_, width_, num_threshold_, b_, hash_table_alpha));
  /// remember that the left ptr must point to the base class in order to call
  /// the methods in it

  /// Step ii. Get ground truth
  ///
  ///       1. read data
  StreamData data(data_file, format); // specify both data file and data format
  if (!data.succeed())
    return;
  Data::GndTruth<key_len, T> gnd_truth, gnd_truth_heavy_hitters;
  gnd_truth.getGroundTruth(data.begin(), data.end(), cnt_method);
  gnd_truth_heavy_hitters.getHeavyHitter(gnd_truth, num_heavy_hitter_,
                                         hx_method);
  ///       2. [optional] show data info
  fmt::print("DataSet: {:d} records with {:d} keys ({})\n", data.size(),
             gnd_truth.size(), data_file);
  /// Step iii. Insert the samples and then look up all the flows
  ///
  ///        1. update records into the sketch
  this->testUpdate(ptr, data.begin(), data.end(),
                   cnt_method); // metrics of interest are in config file
  ///        2. query for all the flowkeys
  if (hx_method == Data::TopK) {
    this->testHeavyHitter(
        ptr, gnd_truth_heavy_hitters.min(),
        gnd_truth_heavy_hitters); // metrics of interest are in config file
  } else {
    this->testHeavyHitter(
        ptr, std::floor(gnd_truth.totalValue() * num_heavy_hitter_ + 1),
        gnd_truth_heavy_hitters); // gnd_truth_heavy_hitter: >, yet HashPipe: >=
  }
  ///        3. size
  this->testSize(ptr);
  ///        3. show metrics
  this->show();

  return;
}

} // namespace OmniSketch::Test

#undef HK_PARA_PATH
#undef HK_TEST_PATH
#undef HK_DATA_PATH

// Driver instance:
//      AUTHOR: XierLabber
//      CONFIG: sketch_config.toml  # with respect to the `src/` directory
//    TEMPLATE: <13, int32_t, Hash::AwareHash>
