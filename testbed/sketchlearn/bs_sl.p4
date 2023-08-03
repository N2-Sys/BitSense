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

typedef bit<32> value_t;
typedef bit<16> index_t;
typedef bit<8> bs_index_t;


#define CM_REGISTER(num) \
    Register<value_t, index_t>(131072) sketch##num;                          \
                                                                            \
    RegisterAction<value_t, index_t, value_t> (sketch##num) update##num = { \
        void apply(inout value_t value) {                                   \
            value = value + 1;                                              \
        }                                                                   \
    }                                                                       \

struct metadata_t {}

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
    /// sketchlearn
    ///
    Hash<index_t>(HashAlgorithm_t.CRC32) hash;

    CM_REGISTER(0);
    CM_REGISTER(1);
    CM_REGISTER(2);
    CM_REGISTER(3);
    CM_REGISTER(4);
    CM_REGISTER(5);
    CM_REGISTER(6);
    CM_REGISTER(7);
    // CM_REGISTER(8);
    // CM_REGISTER(9);
    // CM_REGISTER(10);
    // CM_REGISTER(11);
    // CM_REGISTER(12);
    // CM_REGISTER(13);
    // CM_REGISTER(14);
    // CM_REGISTER(15);
    // CM_REGISTER(16);
    // CM_REGISTER(17);
    // CM_REGISTER(18);
    // CM_REGISTER(19);
    // CM_REGISTER(20);
    // CM_REGISTER(21);
    // CM_REGISTER(22);
    // CM_REGISTER(23);
    // CM_REGISTER(24);
    // CM_REGISTER(25);
    // CM_REGISTER(26);
    // CM_REGISTER(27);
    // CM_REGISTER(28);
    // CM_REGISTER(29);
    // CM_REGISTER(30);
    // CM_REGISTER(31);
    CM_REGISTER(32);

    index_t index;

    apply {
        // No need for egress processing, skip it and use empty controls for egress.
        ig_tm_md.bypass_egress = 1w1;
        // port-based forward
        forward.apply();
        // sketchlearn
        index = hash.get({hdr.ethernet.src_addr});

        update32.execute(index);
        if(hdr.ethernet.src_addr[0:0] & 0x1 == 0x1) {
            update0.execute(index);
        }

        if(hdr.ethernet.src_addr[1:1] & 0x1 == 0x1) {
            update1.execute(index);
        }

        if(hdr.ethernet.src_addr[2:2] & 0x1 == 0x1) {
            update2.execute(index);
        }

        if(hdr.ethernet.src_addr[3:3] & 0x1 == 0x1) {
            update3.execute(index);
        }

        if(hdr.ethernet.src_addr[4:4] & 0x1 == 0x1) {
            update4.execute(index);
        }

        if(hdr.ethernet.src_addr[5:5] & 0x1 == 0x1) {
            update5.execute(index);
        }

        if(hdr.ethernet.src_addr[6:6] & 0x1 == 0x1) {
            update6.execute(index);
        }

        if(hdr.ethernet.src_addr[7:7] & 0x1 == 0x1) {
            update7.execute(index);
        }
/*
        if(hdr.ethernet.src_addr[8:8] & 0x1 == 0x1) {
            update8.execute(index);
        }

        if(hdr.ethernet.src_addr[9:9] & 0x1 == 0x1) {
            update9.execute(index);
        }

        if(hdr.ethernet.src_addr[10:10] & 0x1 == 0x1) {
            update10.execute(index);
        }

        if(hdr.ethernet.src_addr[11:11] & 0x1 == 0x1) {
            update11.execute(index);
        }

        if(hdr.ethernet.src_addr[12:12] & 0x1 == 0x1) {
            update12.execute(index);
        }

        if(hdr.ethernet.src_addr[13:13] & 0x1 == 0x1) {
            update13.execute(index);
        }

        if(hdr.ethernet.src_addr[14:14] & 0x1 == 0x1) {
            update14.execute(index);
        }

        if(hdr.ethernet.src_addr[15:15] & 0x1 == 0x1) {
            update15.execute(index);
        }
*/
/*
        if(hdr.ethernet.src_addr[16:16] & 0x1 == 0x1) {
            update16.execute(index);
        }

        if(hdr.ethernet.src_addr[17:17] & 0x1 == 0x1) {
            update17.execute(index);
        }

        if(hdr.ethernet.src_addr[18:18] & 0x1 == 0x1) {
            update18.execute(index);
        }

        if(hdr.ethernet.src_addr[19:19] & 0x1 == 0x1) {
            update19.execute(index);
        }

        if(hdr.ethernet.src_addr[20:20] & 0x1 == 0x1) {
            update20.execute(index);
        }

        if(hdr.ethernet.src_addr[21:21] & 0x1 == 0x1) {
            update21.execute(index);
        }

        if(hdr.ethernet.src_addr[22:22] & 0x1 == 0x1) {
            update22.execute(index);
        }

        if(hdr.ethernet.src_addr[23:23] & 0x1 == 0x1) {
            update23.execute(index);
        }

        if(hdr.ethernet.src_addr[24:24] & 0x1 == 0x1) {
            update24.execute(index);
        }

        if(hdr.ethernet.src_addr[25:25] & 0x1 == 0x1) {
            update25.execute(index);
        }

        if(hdr.ethernet.src_addr[26:26] & 0x1 == 0x1) {
            update26.execute(index);
        }

        if(hdr.ethernet.src_addr[27:27] & 0x1 == 0x1) {
            update27.execute(index);
        }

        if(hdr.ethernet.src_addr[28:28] & 0x1 == 0x1) {
            update28.execute(index);
        }

        if(hdr.ethernet.src_addr[29:29] & 0x1 == 0x1) {
            update29.execute(index);
        }

        if(hdr.ethernet.src_addr[30:30] & 0x1 == 0x1) {
            update30.execute(index);
        }

        if(hdr.ethernet.src_addr[31:31] & 0x1 == 0x1) {
            update31.execute(index);
        }*/
    }
}

Pipeline(SwitchIngressParser(),
         SwitchIngress(),
         SwitchIngressDeparser(),
         EmptyEgressParser(),
         EmptyEgress(),
         EmptyEgressDeparser()) pipe;

Switch(pipe) main;
