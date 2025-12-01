function [sigma, omega] = berlekamp_massey(syndromes, tables)
    % berlekamp_massey - 표준 Berlekamp-Massey 알고리즘 구현
    %
    %   입력:
    %       syndromes: 1x6 벡터 [S1, S2, S3, S4, S5, S6] (GF 심볼 값)
    %       tables:    GF 테이블 구조체 (exp, log, m)
    %
    %   출력:
    %       sigma:     오류 위치 다항식 σ(x) 계수 [1, σ1, σ2, ...]
    %       omega:     오류 평가 다항식 ω(x) 계수 [ω0, ω1, ...] (Forney용)

    % --- 1. 내부 GF 연산 함수 정의 ---
    gf_add = @(a, b) bitxor(a, b);
    
    function res = gf_mul(a, b)
        if a == 0 || b == 0, res = 0; return; end
        res = tables.exp(mod(double(tables.log(a+1)) + double(tables.log(b+1)), 2^tables.m - 1) + 1);
    end

    function res = gf_div(a, b)
        if b == 0, error('Division by zero in GF'); end
        if a == 0, res = 0; return; end
        % log(a/b) = log(a) - log(b)
        diff = double(tables.log(a+1)) - double(tables.log(b+1));
        if diff < 0, diff = diff + (2^tables.m - 1); end
        res = tables.exp(diff + 1);
    end

    % 다항식 덧셈 (길이가 다르면 0으로 패딩)
    function res = poly_add(p1, p2)
        len = max(length(p1), length(p2));
        p1_pad = [p1, zeros(1, len - length(p1))];
        p2_pad = [p2, zeros(1, len - length(p2))];
        res = bitxor(p1_pad, p2_pad);
        % 끝부분의 0 제거 (필수는 아니지만 최적화용)
        % if length(res) > 1 && res(end) == 0, res(end) = []; end 
    end

    % 다항식 스칼라 곱
    function res = poly_scale(p, s)
        res = arrayfun(@(x) gf_mul(x, s), p);
    end

    % 다항식 시프트 (x^k 곱하기 -> 앞에 0 추가)
    function res = poly_shift(p, k)
        res = [zeros(1, k), p];
    end

    % --- 2. 알고리즘 초기화 ---
    % n = 2t = 6
    % syndromes S = [S1, ..., S6]
    
    sigma = 1;      % σ(x) 초기값 = 1
    omega = 1;      % ω(x) (사용 안 함, 나중에 신드롬으로 계산 가능하지만 여기선 BM 과정에서 업데이트)
                    % Note: 표준 BM에서는 주로 sigma만 구하고 omega는 S*sigma로 따로 구함.
                    % 여기서는 sigma만 정확히 구하는 것에 집중합니다.
    
    B = 1;          % 보조 다항식 B(x) 초기값 = 1
    L = 0;          % 현재 차수 L
    
    % 알고리즘 반복 (k=1 to 2t)
    % 주의: 신드롬 인덱스는 1-based로 S(1)=S1
    
    for n = 0 : 5  % n은 단계 (0부터 2t-1까지, 총 2t=6번 반복)
        
        % (a) Discrepancy (d) 계산
        % d = S_{n+1} + sigma_1*S_n + ... + sigma_L*S_{n+1-L}
        
        d = syndromes(n + 1); % S_{n+1} (첫 항)
        
        for i = 1:L
            % sigma는 [1, sigma_1, sigma_2, ...] 구조
            % sigma_i는 sigma(i+1)에 위치
            % S_{n+1-i}는 syndromes(n+1-i)에 위치
            if (n + 1 - i) > 0
                term = gf_mul(sigma(i + 1), syndromes(n + 1 - i));
                d = gf_add(d, term);
            end
        end
        
        if d == 0
            % Case 1: d = 0
            % B(x) <- x * B(x)
            B = poly_shift(B, 1);
        else
            % Case 2: d != 0
            % T(x) = sigma(x) - d * B(x) (여기서 B(x)는 이미 시프트 된 상태라고 가정하거나 아래에서 처리)
            
            % 이전 sigma 저장 (업데이트 전)
            sigma_prev = sigma;
            
            % sigma(x) 업데이트: sigma(x) - d * x * B(x) 가 아니라
            % 표준 알고리즘에서는 B가 업데이트되면서 x가 곱해져 있음.
            % 여기서는 B를 '이전 단계의 보정항'으로 관리.
            
            % 보정항: d * B(x)
            % 근데 여기서 B(x)는 이전 반복에서의 B(x) * x^1 상태여야 함.
            
            % 정확한 수식:
            % sigma_new = sigma - d * b^-1 * x * B_old
            % Inversionless가 아니므로 b(이전 d)를 저장해야 함.
            % 편의상 책의 알고리즘(Massey)을 따름:
            
            % T(x) = sigma(x) + d * B(x) (단, B는 x가 곱해진 상태)
            term = poly_scale(B, d);    % d * B(x)
            % x * (d*B(x)) -> B는 이미 이전 단계에서 shift 되었어야 함.
            % 로직을 단순화하기 위해 변수 하나를 더 씀.
            
            % --- 간소화된 로직 ---
            % sigma_new = sigma + d * B_shifted
            B_shifted = poly_shift(B, 1); % x * B(x)
            term = poly_scale(B_shifted, d);
            sigma_new = poly_add(sigma, term);
            
            if 2 * L <= n
                % L 업데이트 필요
                L = n + 1 - L;
                % B(x) 업데이트: B(x) = sigma_prev / d
                % 나눗셈 사용 (표준 BM)
                d_inv = gf_div(1, d);
                B = poly_scale(sigma_prev, d_inv);
            else
                % L 유지
                % B(x) = x * B(x)
                B = B_shifted;
                sigma_new = sigma_new; % 그대로
            end
            
            sigma = sigma_new;
        end
    end
    
    % --- 3. Omega(x) 계산 (Forney 알고리즘용) ---
    % Key Equation: Omega(x) = S(x) * Sigma(x) (mod x^2t)
    % S(x) = 1 + S1*x + S2*x^2 + ...
    % 여기서는 S1이 상수항이 아니라 x 계수인 정의를 따름?
    % 보통 Forney용 Omega는: Omega(x) = [S(x) * Sigma(x)] mod x^2t
    % S_poly = [S1, S2, ..., S6]
    
    % 다항식 곱셈 (GF)
    function res = poly_mul(p1, p2)
        len1 = length(p1); len2 = length(p2);
        res = zeros(1, len1 + len2 - 1);
        for ii = 1:len1
            for jj = 1:len2
                term = gf_mul(p1(ii), p2(jj));
                pos = ii + jj - 1;
                res(pos) = gf_add(res(pos), term);
            end
        end
    end

    % Omega 계산
    % S(x) = S1 + S2*x + ... + S6*x^5 (주의: 차수 정의에 따라 다름)
    % 보통 신드롬 다항식 S(x) = S_1 + S_2 x + ... + S_{2t} x^{2t-1}
    S_poly = syndromes; 
    
    full_product = poly_mul(S_poly, sigma);
    
    % mod x^2t (상위 차수 버림)
    % 우리가 필요한건 Omega의 차수가 (2t-1) 이하인 부분
    max_deg = 6; % 2t
    if length(full_product) > max_deg
        omega = full_product(1:max_deg);
    else
        omega = full_product;
    end
    
    % Omega 계산은 정의에 따라 S*Sigma의 상위항인지 하위항인지가 갈림.
    % 일반적인 Forney: Omega(x) = S(x)*Sigma(x) mod x^2t
    
end