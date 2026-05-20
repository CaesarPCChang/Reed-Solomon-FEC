// SPDX-License-Identifier: MIT
`timescale 1ns/1ps

import rs127_121_pkg::*;

// Six-syndrome calculator for RS(127,121), t=3.
// Matches compute_syndromes.m conceptually: evaluate R(alpha^1)..R(alpha^6)
// using Horner form over 127 received symbols.
module rs127_121_syndrome (
  input  logic [N*SYM_W-1:0] code_i,
  output logic [PAR*SYM_W-1:0] syndrome_o,
  output logic                 syndrome_zero_o
);
  gf_t s [0:PAR-1];
  gf_t a [0:PAR-1];
  gf_t r;

  always_comb begin
    for (int j = 0; j < PAR; j++) begin
      s[j] = '0;
      a[j] = gf_alpha_pow(j+1);
    end

    // Treat code_i[126] as highest-degree term and code_i[0] as constant.
    for (int idx = N-1; idx >= 0; idx--) begin
      r = code_i[idx*SYM_W +: SYM_W];
      for (int j = 0; j < PAR; j++) begin
        s[j] = gf_mul(s[j], a[j]) ^ r;
      end
    end

    syndrome_zero_o = 1'b1;
    for (int j = 0; j < PAR; j++) begin
      syndrome_o[j*SYM_W +: SYM_W] = s[j];
      if (s[j] != '0) syndrome_zero_o = 1'b0;
    end
  end
endmodule
