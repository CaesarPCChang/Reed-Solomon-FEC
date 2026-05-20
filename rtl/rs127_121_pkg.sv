// SPDX-License-Identifier: MIT
// Reed-Solomon RS(127,121) helper package.
//
// This package is a SystemVerilog RTL translation of the MATLAB model in this
// repository:
//   - gf_arithmetic.m: GF(2^7) primitive polynomial x^7 + x^3 + 1
//   - compute_syndromes.m: six syndrome values for t=3
//   - berlekamp_massey.m / chien_search.m / forney_algorithm_parallel.m
//
// Symbol convention used by the RTL files:
//   codeword[0] is the earliest/systematic symbol.
//   codeword[126] is the latest symbol.
//   Packed buses use symbol i at bits [i*SYM_W +: SYM_W].

package rs127_121_pkg;
  parameter int SYM_W = 7;
  parameter int N     = 127;
  parameter int K     = 121;
  parameter int PAR   = 6;
  parameter int T     = 3;

  typedef logic [SYM_W-1:0] gf_t;

  // Primitive polynomial: x^7 + x^3 + 1.
  // MATLAB gf_arithmetic.m uses p_x_dec = hex2dec('89').  During modular
  // reduction of x^7, this becomes the lower-term mask x^3 + 1 = 0x09.
  localparam gf_t PRIM_RED = 7'h09;

  // Generator polynomial coefficients g0..g5 from gf_arithmetic.m comments.
  // g(x) = g0 + g1*x + ... + g5*x^5 + x^6.
  localparam gf_t GEN [0:PAR-1] = '{7'h6d, 7'h22, 7'h64, 7'h44, 7'h40, 7'h7e};

  function automatic gf_t gf_add(input gf_t a, input gf_t b);
    gf_add = a ^ b;
  endfunction

  function automatic gf_t gf_mul(input gf_t a, input gf_t b);
    logic [13:0] p;
    begin
      p = '0;
      for (int i = 0; i < SYM_W; i++) begin
        if (b[i]) p[i +: SYM_W] ^= a;
      end
      // Reduce degrees 12 down to 7 using x^7 = x^3 + 1.
      for (int i = 12; i >= 7; i--) begin
        if (p[i]) begin
          p[i]   = 1'b0;
          p[i-4] = p[i-4] ^ 1'b1; // x^(i-7+3)
          p[i-7] = p[i-7] ^ 1'b1; // x^(i-7)
        end
      end
      gf_mul = p[SYM_W-1:0];
    end
  endfunction

  function automatic gf_t gf_pow_elem(input gf_t a, input int unsigned e);
    gf_t result;
    gf_t base;
    int unsigned exp;
    begin
      result = 7'h01;
      base   = a;
      exp    = e;
      while (exp != 0) begin
        if (exp[0]) result = gf_mul(result, base);
        base = gf_mul(base, base);
        exp  = exp >> 1;
      end
      gf_pow_elem = result;
    end
  endfunction

  function automatic gf_t gf_alpha_pow(input int signed e);
    int signed ee;
    gf_t result;
    begin
      ee = e % N;
      if (ee < 0) ee += N;
      result = 7'h01;
      for (int i = 0; i < ee; i++) begin
        result = gf_mul(result, 7'h02);
      end
      gf_alpha_pow = result;
    end
  endfunction

  function automatic gf_t gf_inv(input gf_t a);
    begin
      // In GF(2^7), a^126 = 1 for nonzero a, so a^-1 = a^125.
      gf_inv = (a == '0) ? '0 : gf_pow_elem(a, 125);
    end
  endfunction

  function automatic gf_t gf_div(input gf_t a, input gf_t b);
    begin
      gf_div = (a == '0) ? '0 : gf_mul(a, gf_inv(b));
    end
  endfunction

  function automatic gf_t gf_poly_eval(
    input gf_t poly [0:PAR],
    input int unsigned degree,
    input gf_t x
  );
    gf_t acc;
    begin
      acc = '0;
      for (int i = degree; i >= 0; i--) begin
        acc = gf_mul(acc, x) ^ poly[i];
      end
      gf_poly_eval = acc;
    end
  endfunction
endpackage
