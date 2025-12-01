% 입력/손상/정정 결과 (열벡터라고 가정)
ref_bits      = input_bits;           % 원본
corrupted     = corrupted_bits;       % 에러 주입 후
decoded       = decoded_bits;         % RS 디코드 후

% 에러 마스크
err_injected  = bitxor(ref_bits, corrupted);   % 1이면 주입된 에러 위치
err_remaining = bitxor(ref_bits, decoded);     % 정정 후에도 남은 에러

idx = 1:length(ref_bits);

figure;
stem(idx, err_injected, 'r', 'filled', 'LineWidth', 1.2); hold on;
stem(idx, err_remaining, 'bo', 'LineWidth', 1.2);
ylim([-0.2 1.2]);
xlabel('Bit index'); ylabel('Error mask');
title('3-bit error injection and post-correction');
legend({'Injected errors','Remaining errors'}, 'Location','best');
grid on;
