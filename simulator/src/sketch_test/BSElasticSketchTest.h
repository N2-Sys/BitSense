/**
 * @file BSElasticSketchTest.h
 * @author XierLabber (you@domain.com)
 * @brief Test BS-optimized Elastic Sketch Sketch
 *
 * @copyright Copyright (c) 2022
 *
 */
#include <common/test.h>
#include <sketch/BSElasticSketch.h>

#define BSES_PARA_PATH "ES.para"
#define BSES_TEST_PATH "ES.test"
#define BSES_DATA_PATH "ES.data"
#define BSES_BS_PATH "ES.bs"

namespace OmniSketch::Test {

/**
 * @brief Testing class for Elastic Sketch Sketch
 *
 */
template <int32_t key_len, int32_t no_layer, typename T,
          typename hash_t = Hash::AwareHash>
class BSElasticSketchTest : public TestBase<key_len, T> {
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
  BSElasticSketchTest(const std::string_view config_file)
      : TestBase<key_len, T>("Elastic Sketch with BS", config_file, BSES_TEST_PATH) {
  }

  /**
   * @brief Test BS-optimized Elastic Sketch Sketch
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
void BSElasticSketchTest<key_len, no_layer, T, hash_t>::runTest() {
  /**
   * @brief shorthand for convenience
   *
   */
  using StreamData = Data::StreamData<key_len>;

  /// Part I.
  ///   Parse the config file
  int32_t num_buckets, num_per_bucket, l_depth, l_width;
  double cnt_no_ratio;
  std::vector<size_t> width_cnt, no_hash;
  int32_t heavy_cm_r, heavy_cm_w;
  double cm_cnt_no_ratio;
  std::vector<size_t> cm_width_cnt, cm_no_hash;

  std::string data_file; // data config
  toml::array arr;       // shortly we will convert it to format

  /// Step ii. Open the config file
  Util::ConfigParser parser(config_file);
  if (!parser.succeed()) {
    return;
  }
  /// Step iii. Set the working node of the parser.
  parser.setWorkingNode(
      BSES_PARA_PATH); // do not forget to to enclose it with braces
  /// Step iv. Parse num_bits and num_hash
  if (!parser.parseConfig(num_buckets, "num_buckets"))
    return;
  if (!parser.parseConfig(num_per_bucket, "num_per_bucket"))
    return;
  if (!parser.parseConfig(l_depth, "l_depth"))
    return;
  if (!parser.parseConfig(l_width, "l_width"))
    return;
  /// Step v. Move to the data node
  parser.setWorkingNode(BSES_DATA_PATH);
  /// Step vi. Parse data and format
  if (!parser.parseConfig(data_file, "data"))
    return;
  if (!parser.parseConfig(arr, "format"))
    return;
  /// Step vi. Move to the BS node
  parser.setWorkingNode(BSES_BS_PATH);
  if (!parser.parseConfig(heavy_cm_r, "heavy_cm_r"))
    return;
  if (!parser.parseConfig(heavy_cm_w, "heavy_cm_w"))
    return;
  if (!parser.parseConfig(cnt_no_ratio, "cnt_no_ratio"))
    return;
  if (!parser.parseConfig(width_cnt, "width_cnt"))
    return;
  if (!parser.parseConfig(no_hash, "no_hash"))
    return;
  if (!parser.parseConfig(cm_cnt_no_ratio, "cm_cnt_no_ratio"))
    return;
  if (!parser.parseConfig(cm_width_cnt, "cm_width_cnt"))
    return;
  if (!parser.parseConfig(cm_no_hash, "cm_no_hash"))
    return;

  Data::DataFormat format(arr); // conver from toml::array to Data::DataFormat
  /// [Optional] User-defined rules
  ///
  /// Step vii. Parse Cnt Method.
  parser.setWorkingNode(BSES_DATA_PATH);
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
      new Sketch::BSElasticSketch<key_len, no_layer, T, hash_t>(
          num_buckets, num_per_bucket, l_depth, l_width, cnt_no_ratio,
          width_cnt, no_hash, heavy_cm_r, heavy_cm_w, cm_cnt_no_ratio,
          cm_width_cnt, cm_no_hash));
  /// remember that the left ptr must point to the base class in order to call
  /// the methods in it

  this->testSize(ptr);
  this->show();

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

  printf("\n  BSES NUM BUCKETS: %d\n  BSES NUM PER BUCKET: %d\n  BSES LIGHT DEPTH: %d\n  BSES LIGHT WIDTH: %d\n", num_buckets, num_per_bucket, l_depth, l_width);
  printf("  WIDTH_CNT: [");
  for(int i = 0; i < width_cnt.size();i++)
  {
    printf("%ld", width_cnt[i]);
    if(i != width_cnt.size() - 1)
    {
      printf(", ");
    }
  }
  printf("]\n");
  printf("  RATIO: %lf\n\n", cnt_no_ratio);
  printf("  CM_WIDTH_CNT: [");
  for(int i = 0; i < cm_width_cnt.size();i++)
  {
    printf("%ld", cm_width_cnt[i]);
    if(i != cm_width_cnt.size() - 1)
    {
      printf(", ");
    }
  }
  printf("]\n");
  printf("  CM_RATIO: %lf\n\n", cm_cnt_no_ratio);
  printf("============================================\n");

  return;
}

} // namespace OmniSketch::Test

#undef BSES_PARA_PATH
#undef BSES_TEST_PATH
#undef BSES_DATA_PATH
#undef BSES_BS_PATH

// Driver instance:
//      AUTHOR: XierLabber
//      CONFIG: config/sketch_config.toml  # with respect to the `simulator/` directory
//    TEMPLATE: <13, 2, int32_t, Hash::AwareHash>