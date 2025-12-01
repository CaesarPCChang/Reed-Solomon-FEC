% function [decoded_bits, stats] = rs_decoder_main(encoded_words, tables, crc_enable)
    % rs_decoder_main - Block-level RS-FEC decoder orchestrator
    % encoded_words: 128x1 uint32, each word is 19-bit {is_k, din[17:0]}
    
    % poly_utils = gf_poly_utils(tables);
    utils = gf_poly_utils(tables);

    if numel(encoded_words) ~= 128
        error('encoded_words must have 128 entries');
    end

    % 1) Unpack three slices and is_k vector
    [slice_A, slice_B, slice_C, is_k_vec] = unpack_slices(encoded_wordsCopy);

    % 2) Per-slice decode
    [decoded_A, stats_A] = rs_decode_slice(slice_A, tables);
    [decoded_B, stats_B] = rs_decode_slice(slice_B, tables);
    [decoded_C, stats_C] = rs_decode_slice(slice_C, tables);

    % 3) Repack corrected 121 symbols
    words_19bit = repack_data(decoded_A, decoded_B, decoded_C, is_k_vec);

    % 4) CRC check (optional)
    if crc_enable
        stats.crc_pass = crc_check(words_19bit);
    else
        stats.crc_pass = true;
    end

    % 5) Convert words to bits (LSB-first per word)
    decoded_bits = words_to_bits(words_19bit, crc_enable);

    % 6) Aggregate stats
    stats.corrected_errors = double(stats_A.num_errors) + double(stats_B.num_errors) + double(stats_C.num_errors);
    stats.uncorrectable = logical(stats_A.uncorrectable || stats_B.uncorrectable || stats_C.uncorrectable);
% end

function [sliceA, sliceB, sliceC, is_k] = unpack_slices(words)
    BIT_MASK_18 = uint32(hex2dec('3FFFF'));
    sliceA = zeros(127,1,'uint8');
    sliceB = zeros(127,1,'uint8');
    sliceC = zeros(127,1,'uint8');
    is_k   = zeros(121,1,'uint8');

    % ---------------- Data/CRC (N=0~120) ----------------
    for n = 1:121
        w = uint32(words(n));
        is_k(n) = uint8(bitget(w,19));
        d = bitand(w, BIT_MASK_18);
        msb = bitshift(is_k(n), 6);
        sliceA(n) = bitor(uint8(bitand(bitshift(d,-12), 63)), msb);
        sliceB(n) = bitor(uint8(bitand(bitshift(d,-6),  63)), msb);
        sliceC(n) = bitor(uint8(bitand(d,               63)), msb);
    end

    % 2) Parity 비트스트림 읽기
    parity_bits = zeros(126,1,'uint8');
    bit_pos = 1;
    for word_idx = 122:128
        w = uint32(words(word_idx));
        for bit_in_word = 0:17
            parity_bits(bit_pos) = uint8(bitget(w, bit_in_word+1));
            bit_pos = bit_pos + 1;
        end
    end

    % 7비트씩 묶어서 다시 pA/pB/pC로 분해
    idx = 1;
    for k = 0:5
        % pA
        val = 0;
        for b = 0:6
            val = bitset(val, b+1, parity_bits(idx)); idx=idx+1;
        end
        sliceA(122+k) = uint8(val);

        % pB
        val = 0;
        for b = 0:6
            val = bitset(val, b+1, parity_bits(idx)); idx=idx+1;
        end
        sliceB(122+k) = uint8(val);

        % pC
        val = 0;
        for b = 0:6
            val = bitset(val, b+1, parity_bits(idx)); idx=idx+1;
        end
        sliceC(122+k) = uint8(val);
    end

end

function words = repack_data(A, B, C, is_k)
    % Pack 121 corrected data symbols back to 19-bit words
    % words = zeros(121,1,'uint32');
    %d18 = bitor(bitor(bitshift(uint32(A(1:121)),12), bitshift(uint32(B(1:121)),6)), uint32(C(1:121)));
    %words(1:121) = bitor(d18, bitshift(uint32(is_k),18));
    A_6bit = bitand(uint32(A(1:121)), 63);
    B_6bit = bitand(uint32(B(1:121)), 63);
    C_6bit = bitand(uint32(C(1:121)), 63);

    % 3개의 6-bit 데이터를 18-bit로 재조립
    d18 = bitor(bitor(bitshift(A_6bit, 12), bitshift(B_6bit, 6)), C_6bit);

    % 18-bit 데이터와 1-bit is_k를 합쳐 19-bit 워드 생성
    words = bitor(d18, bitshift(uint32(is_k), 18));
end

function pass = crc_check(words)
    % Compute CRC-18 over first 120 words and compare to word 121 (lower 18 bits)
    if numel(words) < 121, pass=false; return; end
    crc_calc = compute_crc_lfsr(words(1:120));
    crc_word = bitand(uint32(words(121)), uint32(hex2dec('3FFFF')));
    pass = (crc_calc == crc_word);
end

function crc_val = compute_crc_lfsr(data_words_19bit)
    % CRC-18 계산: 18-bit 데이터만 처리 (is_k 비트 제외)
    % 입력: 19-bit 워드 배열 (is_k 포함)
    % 디코더에서는 is_k를 제외한 18-bit만 CRC 계산에 사용
    
    lfsr = uint32(hex2dec('3FFFF'));           % 18-bit all ones
    poly_mask = uint32(hex2dec('BEA7'));       % g(x)=x^18+x^15+...+1
    BIT_MASK_18 = uint32(hex2dec('3FFFF'));
    
    for i = 1:length(data_words_19bit)
        word = uint32(data_words_19bit(i));
        % is_k 제외한 18-bit 데이터만 추출
        data_18bit = bitand(word, BIT_MASK_18);
        
        % Process 18 bits only (din[0]...din[17]), LSB-first
        for j = 0:17  % 18비트만 처리 (is_k 제외)
            inb = bitget(data_18bit, j+1);
            fb  = bitget(lfsr, 18);
            lfsr = bitshift(lfsr,1);
            lfsr = bitset(lfsr,1, bitxor(inb, fb));
            if fb
                lfsr = bitxor(lfsr, poly_mask);
            end
        end
    end
    crc_val = bitand(lfsr, uint32(hex2dec('3FFFF')));
end

function bits = words_to_bits(words, crc_enable)
    % Flatten 19-bit words to a bit vector (LSB-first per word)
    BIT_WIDTH = 19;

    if crc_enable
        num_words = 120; % CRC 모드: 120개 데이터 워드만 변환
    else
        num_words = 121; % CRC 비활성: 121개 데이터 워드 변환
    end

    % 변환할 워드만 선택
    words_to_convert = words(1:num_words);
    N = numel(words_to_convert);

    bits = zeros(N * BIT_WIDTH, 1, 'uint8');
    idx = 1;
    
    for i=1:N
        w = uint32(words(i));
        % for b = 0:18
        for b = BIT_WIDTH-1:-1:0
            bits(idx) = uint8(bitget(w, b+1));
            idx = idx + 1;
        end
    end
end