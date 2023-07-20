/**
 * @file CounterBraidsTest.h
 * @author dromniscience (you@domain.com)
 * @brief Test CounterBraids
 * @copyright Copyright (c) 2022
 *
 */
#include <common/test.h>
#include <sketch/CounterBraids.h>

#define CB_PARA_PATH "CB.para"
#define CB_TEST_PATH "CB.test"
#define CB_DATA_PATH "CB.data"

namespace OmniSketch::Test {

/**
 * @brief Testing class for Count Min Sketch
 *
 */
template <int32_t key_len, int32_t no_layer, typename T,
          typename hash_t = Hash::AwareHash>
class CounterBraidsTest : public TestBase<key_len, T> {
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
  CounterBraidsTest(const std::string_view config_file)
      : TestBase<key_len, T>("Counter Braids", config_file, CB_TEST_PATH) {}

  /**
   * @brief Test CH-optimized Count Min Sketch
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

template <int32_t key_len, int32_t no_layer, typename T, typename hash_t>
void CounterBraidsTest<key_len, no_layer, T, hash_t>::runTest() {
  /**
   * @brief shorthand for convenience
   *
   */
  using StreamData = Data::StreamData<key_len>;

  /// Part I.
  ///   Parse the config file
  std::vector<size_t> no_cnt, width_cnt, no_hash;

  std::string data_file; // data config
  toml::array arr;       // shortly we will convert it to format

  /// Step ii. Open the config file
  Util::ConfigParser parser(config_file);
  if (!parser.succeed()) {
    return;
  }
  /// Step iii. Set the working node of the parser.
  parser.setWorkingNode(CB_PARA_PATH);
  if (!parser.parseConfig(no_cnt, "no_cnt"))
    return;
  if (!parser.parseConfig(width_cnt, "width_cnt"))
    return;
  if (!parser.parseConfig(no_hash, "no_hash"))
    return;
  /// Step v. Move to the data node
  parser.setWorkingNode(CB_DATA_PATH);
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
  parser.setWorkingNode(CB_DATA_PATH);
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
      new Sketch::CounterBraids<key_len, no_layer, T, hash_t>(no_cnt, width_cnt,
                                                              no_hash));
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

  Data::GndTruth<key_len, T> heavy_part;
  heavy_part.getHeavyHitter(gnd_truth, gnd_truth.size() * 0.3, Data::TopK);
  printf("SIZE: %ld\n", heavy_part.size());

  ///       2. [optional] show data info
  // heavy_part.getHeavyHitter(gnd_truth, 0.3, Data::Percentile);
  fmt::print("DataSet: {:d} records with {:d} keys ({})\n", data.size(),
             gnd_truth.size(), data_file);
  /// Step iii. Insert the samples and then look up all the flows
  ///
  ///        1. update records into the sketch
  this->testUpdate(ptr, data.begin(), data.end(),
                   cnt_method); // metrics of interest are in config file
  ///        2. decode
  this->testDecode(ptr, heavy_part); // metrics of interest are in config file
  // this->testDecode(ptr, gnd_truth); // metrics of interest are in config file
  ///        3. size
  this->testSize(ptr);
  ///        3. show metrics
  this->show();

  return;
}

} // namespace OmniSketch::Test

#undef CB_PARA_PATH
#undef CB_TEST_PATH
#undef CB_DATA_PATH

// Driver instance:
//      AUTHOR: dromniscience
//      CONFIG: sketch_config.toml  # with respect to the `src/` directory
//    TEMPLATE: <13, 2, int32_t, Hash::AwareHash>
