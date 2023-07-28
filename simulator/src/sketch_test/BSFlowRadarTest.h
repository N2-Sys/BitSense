/**
 * @file BSFlowRadarTest.h
 * @author XierLabber
 * @brief Test BS Flow Radar
 *
 * @copyright Copyright (c) 2022
 *
 */
#pragma once

#include <common/test.h>
#include <common/hash.h>
#include <common/sketch.h>
#include <sketch/BSFlowRadar.h>

#define BSFR_PARA_PATH "FlowRadar.para"
#define BSFR_TEST_PATH "FlowRadar.test"
#define BSFR_DATA_PATH "FlowRadar.data"
#define BSFR_BS_PATH "FlowRadar.bs"

namespace OmniSketch::Test {

/**
 * @brief Testing class for BS Flow Radar
 *
 */
template <int32_t key_len, int32_t no_layer, typename T,
          typename hash_t>
class BSFlowRadarTest : public TestBase<key_len, T> {
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
  BSFlowRadarTest(const std::string_view config_file)
      : TestBase<key_len, T>("Flow Radar with BS", config_file, BSFR_TEST_PATH) {}

  /**
   * @brief Test BS Flow Radar
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

template <int32_t key_len, int32_t no_layer, typename T,
          typename hash_t>
void BSFlowRadarTest<key_len, no_layer, T, hash_t>::runTest(){
  // for convenience only
  using StreamData = Data::StreamData<key_len>;

  // parse config
  int32_t flow_filter_bit, flow_filter_hash, count_table_num,
      count_table_hash;  // sketch config
  double flow_cnt_no_ratio, packet_cnt_no_ratio;
  std::vector<size_t> flow_width_cnt, flow_no_hash, packet_width_cnt, packet_no_hash;
  std::string data_file; // data config
  toml::array arr;       // shortly we will convert it to format

  Util::ConfigParser parser(config_file);
  if (!parser.succeed()) {
    return;
  }

  parser.setWorkingNode(BSFR_PARA_PATH);
  if (!parser.parseConfig(flow_filter_bit, "flow_filter_bit"))
    return;
  if (!parser.parseConfig(flow_filter_hash, "flow_filter_hash"))
    return;
  if (!parser.parseConfig(count_table_num, "count_table_num"))
    return;
  if (!parser.parseConfig(count_table_hash, "count_table_hash"))
    return;

  // prepare data
  parser.setWorkingNode(BSFR_DATA_PATH);
  if (!parser.parseConfig(data_file, "data"))
    return;
  if (!parser.parseConfig(arr, "format"))
    return;

  parser.setWorkingNode(BSFR_BS_PATH);
  if (!parser.parseConfig(flow_cnt_no_ratio, "flow_cnt_no_ratio"))
    return;
  if (!parser.parseConfig(flow_width_cnt, "flow_width_cnt"))
    return;
  if (!parser.parseConfig(flow_no_hash, "flow_no_hash"))
    return;
  if (!parser.parseConfig(packet_cnt_no_ratio, "packet_cnt_no_ratio"))
    return;
  if (!parser.parseConfig(packet_width_cnt, "packet_width_cnt"))
    return;
  if (!parser.parseConfig(packet_no_hash, "packet_no_hash"))
    return;

  Data::DataFormat format(arr);
  StreamData data(data_file, format);
  if (!data.succeed())
    return;
  Data::GndTruth<key_len, T> gnd_truth;
  gnd_truth.getGroundTruth(data.begin(), data.end(), Data::InPacket);
  fmt::print("DataSet: {:d} records with {:d} keys ({})\n", data.size(),
             gnd_truth.size(), data_file);

  std::unique_ptr<Sketch::SketchBase<key_len, T>> ptr(
      new Sketch::BSFlowRadar<key_len, no_layer, T, hash_t>(
          flow_filter_bit, flow_filter_hash, count_table_num,
          count_table_hash, flow_cnt_no_ratio, flow_width_cnt,
          flow_no_hash, packet_cnt_no_ratio, packet_width_cnt,
          packet_no_hash));

  this->testSize(ptr);
  this->show();

  this->testUpdate(ptr, data.begin(), data.end(), Data::InPacket);
  this->testDecode(ptr, gnd_truth);
  this->testSize(ptr);
  // show
  this->show();

  printf("\n  FLOW FILTER BIT SIZE: %d\n  FLOW FILTER HASH NUM: %d\n  COUNT TABLE NUM: %d\n  COUNT HASH NUM: %d\n", flow_filter_bit, flow_filter_hash, count_table_num, count_table_hash);
  printf("  FLOW WIDTH_CNT: [");
  for(int i = 0; i < flow_width_cnt.size();i++)
  {
    printf("%ld", flow_width_cnt[i]);
    if(i != flow_width_cnt.size() - 1)
    {
      printf(", ");
    }
  }
  printf("]\n");
  printf("  FLOW BS RATIO: %lf\n", flow_cnt_no_ratio);
  printf("  PACKET WIDTH_CNT: [");
  for(int i = 0; i < packet_width_cnt.size();i++)
  {
    printf("%ld", packet_width_cnt[i]);
    if(i != packet_width_cnt.size() - 1)
    {
      printf(", ");
    }
  }
  printf("]\n");
  printf("  PACKET BS RATIO: %lf\n\n", packet_cnt_no_ratio);    
  printf("============================================\n");


  return;
}

} // namespace OmniSketch::Test

#undef BSFR_PARA_PATH
#undef BSFR_TEST_PATH
#undef BSFR_DATA_PATH
#undef BSFR_BS_PATH

// Driver instance:
//      AUTHOR: dromniscience
//      CONFIG: config/sketch_config.toml  # with respect to the `simulator/` directory
//    TEMPLATE: <13, 2, int32_t, Hash::AwareHash>