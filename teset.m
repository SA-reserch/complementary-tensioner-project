clear;
clc;
%---------固定値の代入及び計算---------%
%固定値
R_1=9.55;
R_2=19.1;
r_1=5;
r_2=5;
max_count = 1000;
%右プーリの座標hita
p_left=[0;0];
%右プーリの座標
p_right=[230;0];
%曲線の始点座標
c_start=[170;-25];
%図心の始点座標
p_centre=[c_start(1);c_start(2)-r_1];
%上の曲線のパラメータ

%右プーリから左プーリへのベクトル
vec_a = p_left - p_right;
%正規直交基底の接線方向ベクトル
vec_e = [-1;0];
%正規直行基底の法線ベクトル
n = [0;-1];
%(-pi/2)の回転行列
Rot_90= [0 , 1;
         -1, 0];
L = 230;
%-----------逐次処理------------%
%右プーリから図心へのベクトル（初期値）
p_0 = p_centre - p_right;
%右プーリから曲線の始点へのベクトル（初期値）
c_0= c_start - p_right;
%下のタイベル左側の巻き出し初期値
l_length_left_0 = low_length_left(p_0, R_1, r_1, vec_a); 
%下のタイベル総長初期値
l_length_all_0 = low_length(p_0, R_1, R_2, r_1, vec_a);
%Δ幅を設定
ds = 0.15;
c=c_0;
p=p_0;




% 初期x位置
x0 = 10;

% f(x) 定義
f = @(x) 20*sin((x - 10)/200*2*pi).^2 + 10 + 20*sin((x - 10)/400*2*pi);

% 安全な関数
arccos_safe = @(x) real(acos(max(min(x, 1), -1)));
sqrt_safe = @(x) sqrt(max(x, 0));

% l1(x), l2(x) 定義（ここ注意！r_1, r_2に統一）
l1 = @(x) r_1 * (pi - atan2(f(x), x) - arccos_safe(r_1 ./ sqrt(x.^2 + f(x).^2))) + ...
          sqrt_safe(x.^2 + f(x).^2 - r_1^2);

l2 = @(x) r_2 * (pi - atan2(f(x), L - x) - arccos_safe(r_2 ./ sqrt((L - x).^2 + f(x).^2))) + ...
          sqrt_safe((L - x).^2 + f(x).^2 - r_2^2);

% 初期長さを計算
l1_0 = l1(x0);
l2_0 = l2(x0);




for count = 1:max_count
    % cを更新
    c = c + ds * vec_e;
    
    % --- 最適化関数（残差＆幾何拘束） ---
    epsilon = 1e-5;  % 微小値：方向を指定するため
    %全体の長さの変化を0に抑える．
    %下部の長さの変化量+上部の長さの変化量=0：の計算式を導入
    func = @(p_var) [
    low_length(p_var, R_1, R_2, r_1, vec_a) - l_length_all_0 ...
        + delta_total_from_delta_l2(low_length_left(p_var, R_1, r_1, vec_a) - l_length_left_0, R_1, R_2, L);
    norm(p_var - c) - r_2;                      % 半径制約
    p_var(1) - c(1) + epsilon;                  % x方向：左（p_x < c_x）
    %p_var(2) - c(2) + epsilon                   % y方向：下（p_y < c_y）
    ];


    
    % --- 初期推定値は前のnの方向にcをr_1だけ伸ばした点 ---
    p_init = c + r_1*n;
    
    % --- fsolve設定＆求解 ---
    options = optimoptions('fsolve','Display', 'none','Algorithm', 'levenberg-marquardt','TolX', 1e-8);
    p = fsolve(func, p_init, options);
    
    % cとpから法線ベクトルを求める
    n = (p-c)/norm(p-c);
    % 法線ベクトルnから接線ベクトルeを求める
    vec_e = Rot_90 * n;

    % 各ベクトルを記録
    p_list(:, count) = p;
    c_list(:, count) = c;
    n_list(:, count) = n;
    e_list(:, count) = vec_e;

     % === ここから追加する処理 ===
    % 今のpについて、Low_length（下部ベルト長）を計算
    low_len_now = low_length(p, R_1, R_2, r_1, vec_a);

    % 今のpについて、delta_total_from_delta_l2を使って上部ベルト長変化を計算
    delta_total_now = delta_total_from_delta_l2(...
        low_length_left(p, R_1, r_1, vec_a) - l_length_left_0...%下部の左変化量
        , r_1, r_2, L);

    % 上部ベルト長（絶対量）を取得（注意：delta_totalは変化量なので！）
    upper_belt_now = (l1_0 + l2_0) + delta_total_now;  % ←これがポイント

    % 下部ベルト長と上部ベルト長を足して、全体長にする
    total_belt_now = low_len_now + upper_belt_now;

    % 記録
    total_belt_length_list(count) = total_belt_now;

    % 下部と上部それぞれ記録
    lower_belt_length_list(count) = low_len_now;
    upper_belt_length_list(count) = upper_belt_now;
    
    % 全体ベルト長さも記録
    total_belt_length_list(count) = total_belt_now;
    % --- 変化量を計算 ---
    delta_lower_belt_length_list = lower_belt_length_list - lower_belt_length_list(1);
    delta_upper_belt_length_list = upper_belt_length_list - upper_belt_length_list(1);
    delta_total_belt_length_list = total_belt_length_list - total_belt_length_list(1);


    end
% 逐次の c - a ベクトルを保存
p_minus_a_list = c_list - vec_a;  % a は列ベクトル [2×1]、p_list は [2×N]

figure;
plot(p_minus_a_list(1, :), p_minus_a_list(2, :), '-o', 'LineWidth', 1.5);
axis equal;
xlabel('x (mm)');
ylabel('y (mm)');
title('p - a vector trajectory');
grid on;

figure;
plot(1:max_count, total_belt_length_list, 'LineWidth', 1.5);
xlabel('Step');
ylabel('Total Belt Length (mm)');
title('Evolution of Total Belt Length Over Time');
grid on;

%==== 下部タイベルの変化量グラフ ====
figure;
plot(1:max_count, delta_lower_belt_length_list, '-o', 'LineWidth', 1.5);
xlabel('Step');
ylabel('Lower Belt Length Change (mm)');
title('Change in Lower Belt Length Over Time');
grid on;

%==== 上部タイベルの変化量グラフ ====
figure;
plot(1:max_count, delta_upper_belt_length_list, '-o', 'LineWidth', 1.5);
xlabel('Step');
ylabel('Upper Belt Length Change (mm)');
title('Change in Upper Belt Length Over Time');
grid on;



%vec_pに対するその時の下部のタイベルの総長が分かる
function Low_length = low_length(vec_p, R_1, R_2, r_1, vec_a)
    theta_1 = acos((R_1 - r_1) / norm(vec_p));
    theta_2 = acos(dot(vec_a, vec_p) / (norm(vec_a) * norm(vec_p)));
    theta_3 = acos(dot(-vec_a, -vec_a + vec_p) / (norm(-vec_a) * norm(-vec_a + vec_p)));
    theta_4 = acos((R_2 - r_1) / norm(vec_a - vec_p));

    Low_length = R_1*(pi - theta_1 - theta_2) + ...
                 R_2*(pi - theta_3 - theta_4) + ...
                 r_1*(theta_1 + theta_2 + theta_3 + theta_4 - pi) + ...
                 norm(vec_p)*sin(theta_1) + ...
                 norm(-vec_a + vec_p)*sin(theta_4);
end


%vec_pに対するその時のタイベルの巻き出し量が分かる
function Low_length_left = low_length_left(vec_p, R_1, r_1, a)
    theta_1 = acos((R_1 - r_1) / norm(vec_p));
    theta_2 = acos(dot(a, vec_p) / (norm(a) * norm(vec_p)));
    Low_length_left = R_1*(pi - theta_1 - theta_2) + ...
                      r_1*(theta_1 + theta_2 - pi/2) + ...
                      norm(vec_p)*sin(theta_1);
end
%タイベルの巻き出し量から上のタイベルの総長変化量を求める
function delta_total = delta_total_from_delta_l2(delta_l2_target, r1, r2, L)
    f = @(x) 20*sin((x - 10)/200*2*pi).^2 + 10 + 20*sin((x - 10)/400*2*pi);
    arccos_safe = @(x) real(acos(max(min(x, 1), -1)));
    sqrt_safe = @(x) sqrt(max(x, 0));

    l1 = @(x) r1 * (pi - atan2(f(x), x) - arccos_safe(r1 ./ sqrt(x.^2 + f(x).^2))) + ...
              sqrt_safe(x.^2 + f(x).^2 - r1^2);
    l2 = @(x) r2 * (pi - atan2(f(x), L - x) - arccos_safe(r2 ./ sqrt((L - x).^2 + f(x).^2))) + ...
              sqrt_safe((L - x).^2 + f(x).^2 - r2^2);

    x0 = 10;
    l1_0 = l1(x0);
    l2_0 = l2(x0);
    l2_target = l2_0 - delta_l2_target;%下部の伸びを打ち消す必要があるので-delta

    cost = @(x) abs(l2(x) - l2_target);
    options = optimset('TolX', 1e-12);
    x_best = fminbnd(cost, 10, 210, options);

    delta_total = (l1(x_best) + l2(x_best)) - (l1_0 + l2_0);

    % 現在のワークスペースに 'x_best_list' が存在するか確認
    if evalin('base', 'exist(''x_best_list'', ''var'')')
    % 既に配列があるなら、末尾に追加
    assignin('base', 'x_best_list', [evalin('base', 'x_best_list'), x_best]);
    else
    % 初回なら、新しく作る
    assignin('base', 'x_best_list', x_best);
    end
    % ==== delta_total_listにも保存 ====
    if evalin('base', 'exist(''delta_total_list'', ''var'')')
    % 既に配列があれば末尾に追加
    assignin('base', 'delta_total_list', [evalin('base', 'delta_total_list'), delta_total]);
    else
    % なければ新しく作成
    assignin('base', 'delta_total_list', delta_total);
    end


end



