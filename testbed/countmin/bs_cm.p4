/*******************************************************************************
 *  INTEL CONFIDENTIAL
 *
 *  Copyright (c) 2021 Intel Corporation
 *  All Rights Reserved.
 *
 *  This software and the related documents are Intel copyrighted materials,
 *  and your use of them is governed by the express license under which they
 *  were provided to you ("License"). Unless the License provides otherwise,
 *  you may not use, modify, copy, publish, distribute, disclose or transmit
 *  this software or the related documents without Intel's prior written
 *  permission.
 *
 *  This software and the related documents are provided as is, with no express
 *  or implied warranties, other than those that are expressly stated in the
 *  License.
 ******************************************************************************/


#include <core.p4>
#if __TARGET_TOFINO__ == 3
#include <t3na.p4>
#elif __TARGET_TOFINO__ == 2
#include <t2na.p4>
#else
#include <tna.p4>
#endif

typedef bit<8> value_t;
typedef bit<17> index_t;
typedef bit<14> bs_index_t;


#define CM_REGISTER_AND_HASH(num, seed) \
    Register<value_t, index_t>(131072) sketch##num;                          \
                                                                            \
    RegisterAction<value_t, index_t, value_t> (sketch##num) update##num = { \
        void apply(inout value_t value) {                                   \
            value = value + 1;                                              \
        }                                                                   \
    };                                                                      \
                                                                            \
    CRCPolynomial<bit<32>>(seed,                                            \
                           true,                                            \
                           false,                                           \
                           false,                                           \
                           32w0xFFFFFFFF,                                   \
                           32w0xFFFFFFFF                                    \
                           ) poly##num;                                     \
    Hash<index_t>(HashAlgorithm_t.CUSTOM, poly##num) hash##num;             \
                                                                            \
    index_t tmp_index##num;                                                 \
    value_t value##num

#define APPLY_HASH(num)                   \
    tmp_index##num = hash##num.get({          \
      hdr.ipv4.src_addr,                  \
      hdr.ipv4.dst_addr,                  \
      hdr.tcp.src_port,                   \
      hdr.tcp.dst_port,                   \
      hdr.ipv4.protocol                   \
    });                                   \
    value##num = update##num.execute(tmp_index##num); \
    ig_md.index##num = tmp_index##num

#define BITSENSE(num, seed)                      \
    Register<value_t, bs_index_t>(40000) bs_sketch##num; \
      \
    RegisterAction<value_t, bs_index_t, value_t>(bs_sketch##num) bs_update##num = { \
        void apply(inout value_t value) {           \
            value = value + 1;                                              \
        }                                                                   \
    };  \
      \
    CRCPolynomial<bit<32>>(seed,                                            \
                           true,                                            \
                           false,                                           \
                           false,                                           \
                           32w0xFFFFFFFF,                                   \
                           32w0xFFFFFFFF                                    \
                           ) bs_poly##num;                                  \
    Hash<bs_index_t>(HashAlgorithm_t.CUSTOM, bs_poly##num) bs_hash##num;    \
                                                                            \
    bs_index_t bs_index##num

#define UPDATE_BITSENSE(num)                   \
    bs_index##num = bs_hash##num.get(          \
      {ig_md.index##num} \
    );                                   \
    bs_update##num.execute(bs_index##num)


struct metadata_t {
  index_t index0;
  index_t index1;
  index_t index2;
}

#include "parser.p4"


control SwitchIngress(
        inout header_t hdr,
        inout metadata_t ig_md,
        in ingress_intrinsic_metadata_t ig_intr_md,
        in ingress_intrinsic_metadata_from_parser_t ig_prsr_md,
        inout ingress_intrinsic_metadata_for_deparser_t ig_dprsr_md,
        inout ingress_intrinsic_metadata_for_tm_t ig_tm_md) {
    
    ///
    /// port-based forwarding
    ///
    action hit(PortId_t port) {
        ig_tm_md.ucast_egress_port = port;
    }

    action miss(bit<3> drop) {
        ig_dprsr_md.drop_ctl = drop; // Drop packet.
    }

    table forward {
        key = {
            ig_intr_md.ingress_port : exact;
        }

        actions = {
            hit;
            @defaultonly miss;
        }
        const default_action = miss(0x1);
        size = 16;
    }

    ///
    /// count-min sketch
    ///
    CM_REGISTER_AND_HASH(0, 32w0x04C11DB7);
    CM_REGISTER_AND_HASH(1, 32w0x34FD110C);
    CM_REGISTER_AND_HASH(2, 32w0x203AD4E3);

    // BitSense
    BITSENSE(0, 32w0x20230117);
    BITSENSE(1, 32w0x20230116);
    BITSENSE(2, 32w0x20230115);

    apply {
        // No need for egress processing, skip it and use empty controls for egress.
        ig_tm_md.bypass_egress = 1w1;
        // port-based forward
        forward.apply();
        // count-min
        if(hdr.tcp.isValid()){
            APPLY_HASH(0);
            if(value0 == (value_t)0xFF){
                UPDATE_BITSENSE(0);
            }

            APPLY_HASH(1);
            if(value1 == (value_t)0xFF){
                UPDATE_BITSENSE(1);
            }

            APPLY_HASH(2);
            if(value2 == (value_t)0xFF){
                UPDATE_BITSENSE(2);
            }
        }
    }
}

Pipeline(SwitchIngressParser(),
         SwitchIngress(),
         SwitchIngressDeparser(),
         EmptyEgressParser(),
         EmptyEgress(),
         EmptyEgressDeparser()) pipe;

Switch(pipe) main;
