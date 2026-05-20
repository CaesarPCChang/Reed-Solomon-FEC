// SPDX-License-Identifier: MIT
`timescale 1ns/1ps

import rs127_121_pkg::*;

// Systematic RS(127,121) encoder over GF(2^7).
//
// Input : 121 symbols, packed as data_i[symbol_index*7 +: 7].
// Output: 127 symbols = 121 systematic symbols followed by 6 parity symbols.
//
// This is intended as the first hardware anchor for the MATLAB algorithm in
// this repo. For JPEG XS transport, instantiate this block per byte-lane or
// symbol-lane after an external packet/interleaver maps RTP payload bytes into
// 7-bit RS symbols.
module rs127_121_encoder_block (
  input  logic [K*SYM_W-1:0] data_i,
  output logic [N*SYM_W-1:0] code_o
);
  gf_t parity [0:PAR-1];
  gf_t next_p [0:PAR-1];
  gf_t fb;
  gf_t d;

  always_comb begin
    for (int p = 0; p < PAR; p++) begin
      parity[p] = '0;
    end

    // LFSR polynomial division by g(x) = x^6 + g5*x^5 + ... + g0.
    for (int i = 0; i < K; i++) begin
      d  = data_i[i*SYM_W +: SYM_W];
      fb = d ^ parity[PAR-1];

      next_p[PAR-1] = parity[PAR-2] ^ gf_mul(fb, GEN[PAR-1]);
      next_p[PAR-2] = parity[PAR-3] ^ gf_mul(fb, GEN[PAR-2]);
      next_p[PAR-3] = parity[PAR-4] ^ gf_mul(fb, GEN[PAR-3]);
      next_p[PAR-4] = parity[PAR-5] ^ gf_mul(fb, GEN[PAR-4]);
      next_p[PAR-5] = parity[PAR-6] ^ gf_mul(fb, GEN[PAR-5]);
      next_p[0]     =              gf_mul(fb, GEN[0]);

      for (int p = 0; p < PAR; p++) begin
        parity[p] = next_p[p];
      end
    end

    // Systematic codeword layout.
    code_o = '0;
    for (int i = 0; i < K; i++) begin
      code_o[i*SYM_W +: SYM_W] = data_i[i*SYM_W +: SYM_W];
    end
    for (int p = 0; p < PAR; p++) begin
      code_o[(K+p)*SYM_W +: SYM_W] = parity[p];
    end
  end
endmodule
