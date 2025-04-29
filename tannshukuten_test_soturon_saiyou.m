% 定数の設定
syms r1 r2 L

r1 = 9.55; % 半径(HTPA20)
r2 = 19.1; % 半径(HTPA40)
L = 230; % 全長

% 変数設定
% R1 = 14;  % 半径R1
Rr1=9.55;%(HTPA20)
Rr2=6;%(外径12mmベアリング)
k1=35;
k2=35.42937455;%(腕)


% x と f(x) を連続関数として定義
x_values = linspace(10, 210, 100); % x の範囲
f_values = 20*sin((x_values-10)/200*2*pi).^2 + 10 + 20*sin((x_values-10)/400*2*pi); % f(x) の値

% f をインデックス付き関数として定義
f = @(x_idx) f_values(x_idx);

% l1 と l2 の計算式
l1 = @(x_val, f_val) r1 * (pi - atan2(f_val, x_val) - acos(r1 ./ sqrt(x_val.^2 + f_val.^2))) + ...
                     sqrt(x_val.^2 + f_val.^2 - r1^2);

l2 = @(x_val, f_val) r2 * (pi - atan2(f_val, L - x_val) - acos(r2 ./ sqrt((L - x_val).^2 + f_val.^2))) + ...
                     sqrt((L - x_val).^2 + f_val.^2 - r2^2);

% 初期値の計算 (x = 10)
x_0 = 10;
f_0 = f(1);  % f_0 を前半と同様に定義
l1_0 = l1(x_0, f_0);
l2_0 = l2(x_0, f_0);

% delta_deta の計算（delta_dataとして使用）
l1_values = l1(x_values, f_values);
l2_values = l2(x_values, f_values);
delta_deta = (l1_values + l2_values) - (l1_0 + l2_0);  % delta_detaを計算

% delta_dataとしてdelta_detaを使用
delta_data = delta_deta;  % delta_dataに代入（もし別途設定する必要があれば）

% 以降の処理ではdelta_dataを使用します
a_values = zeros(1, length(delta_data));  % aの値を保存する配列
l_checks = zeros(1, length(delta_data));
a_0=pi/3;
% l_0 = 2*R1*(pi/2 - acos(2*R1/k2) + asin((k2*sin(a_0))/sqrt(k1^2 + k2^2 - 2*k1*k2*cos(a_0)))) ...
%     + 2*R1*(pi/2 - acos(2*R1/k2) + a_0) + k2*sin(acos(2*R1/k2)) + sqrt(k1^2 + k2^2 - 2*k1*k2*cos(a_0))*sin(acos(2*R1/k2));
l_0= ...
    k2*(1 - (Rr1 + Rr2)^2/k2^2)^(1/2) + ...
    (1 - (Rr1 + Rr2)^2/(k1^2 - 2*cos(a_0)*k1*k2 + k2^2))^(1/2)*(k1^2 - 2*cos(a_0)*k1*k2 + k2^2)^(1/2) + ...
    Rr1*(a_0 + pi - acos((Rr1 + Rr2)/(k1^2 - 2*cos(a_0)*k1*k2 + k2^2)^(1/2)) - acos((Rr1 + Rr2)/k2) + atan2((k2*sin(a_0)),(k1 - k2*cos(a_0)))) + ...
    Rr2*(a_0 + pi - acos((Rr1 + Rr2)/(k1^2 - 2*cos(a_0)*k1*k2 + k2^2)^(1/2)) - acos((Rr1 + Rr2)/k2) + atan2((k2*sin(a_0)),(k1 - k2*cos(a_0))));

l = l_0;  % lの初期値

% 逐次的にaを計算する
for i = 1:length(delta_data)
    % 各deltaに対応するlを計算
    l = l_0 - delta_data(i);  % l_0からdelta_data(i)を引く
    disp(l); 
    % 初期値aを数値として設定（この部分が問題）
    a = double(pi/3);  % aの初期値を数値に変換
    
    % lに基づいてaを求める
   
    opt_func = @(a_new)abs(l-(...
    k2*(1 - (Rr1 + Rr2)^2/k2^2)^(1/2) + ...
    (1 - (Rr1 + Rr2)^2/(k1^2 - 2*cos(a_new)*k1*k2 + k2^2))^(1/2)*(k1^2 - 2*cos(a_new)*k1*k2 + k2^2)^(1/2) + ...
    Rr1*(a_new + pi - acos((Rr1 + Rr2)/(k1^2 - 2*cos(a_new)*k1*k2 + k2^2)^(1/2)) - acos((Rr1 + Rr2)/k2) + atan2((k2*sin(a_new)),(k1 - k2*cos(a_new)))) + ...
    Rr2*(a_new + pi - acos((Rr1 + Rr2)/(k1^2 - 2*cos(a_new)*k1*k2 + k2^2)^(1/2)) - acos((Rr1 + Rr2)/k2) + atan2((k2*sin(a_new)),(k1 - k2*cos(a_new))))));
    % 最適化関数を最小化することでaを求める
    a = fminsearch(opt_func, a);  % 初期推定値a
    a_values(i) = a;  % 更新されたaを保存
    l_checks(i) = k2*(1 - (Rr1 + Rr2)^2/k2^2)^(1/2) + (1 - (Rr1 + Rr2)^2/(k1^2 - 2*cos(a)*k1*k2 + k2^2))^(1/2)*(k1^2 - 2*cos(a)*k1*k2 + k2^2)^(1/2) + Rr1*(a + pi - acos((Rr1 + Rr2)/(k1^2 - 2*cos(a)*k1*k2 + k2^2)^(1/2)) - acos((Rr1 + Rr2)/k2) + atan2((k2*sin(a)),(k1 - k2*cos(a)))) + Rr2*(a + pi - acos((Rr1 + Rr2)/(k1^2 - 2*cos(a)*k1*k2 + k2^2)^(1/2)) - acos((Rr1 + Rr2)/k2) + atan2((k2*sin(a)),(k1 - k2*cos(a))));
end

disp(l_checks);
% ラジアンの値も表示
% disp('各ステップで求められたaの値（ラジアン単位）:');
% disp(a_values);  % ラジアン単位での値
% 
% 度の値を計算して表示
a_values_deg = rad2deg(a_values);  % ラジアンから度に変換
% disp('各ステップで求められたaの値（度単位）:');
% disp(a_values_deg);  % 度単位での値


% k2 * cos(a_values) と k2 * sin(a_values) を計算
x_values_cos = k2 * cos(a_values);
% disp(x_values_cos);
y_values_sin = k2 * sin(a_values);

% グラフのプロット
% figure;
% plot(x_values_cos, y_values_sin, 'o-', 'LineWidth', 1.5, 'DisplayName', 'k2 * sin(a) vs k2 * cos(a)');
% xlabel('k2 * cos(a)', 'FontSize', 12);
% ylabel('k2 * sin(a)', 'FontSize', 12);
% title('k2 * cos(a) と k2 * sin(a) の関係', 'FontSize', 14);
% grid on;
% legend('show');

% 定数の設定
L_2 = 25; % 変更
u = 1/2; % 減速比

% L2 の計算
l2 = @(x_idx) r2 * (pi - atan2(f(x_idx), L - x_values(x_idx)) - acos(r2 ./ sqrt((L - x_values(x_idx)).^2 + f(x_idx).^2))) + ...
              sqrt((L - x_values(x_idx)).^2 + f(x_idx).^2 - r2^2);

l2_0 = r2 * (pi - atan2(f(1), L - x_values(1)) - acos(r2 ./ sqrt((L - x_values(1)).^2 + f(1).^2))) + ...
         sqrt((L - x_values(1)).^2 + f(1).^2 - r2^2);

% L2 の計算
L2 = @(x_idx) L_2 + u * (l2_0 - l2(x_idx));

% L2 の結果を表示
disp('L2 の値:');
L2_values = arrayfun(L2, 1:length(x_values));
% disp(L2_values);

% L2_values と k2 * cos(a) を計算する
x_values_L2 = L2_values;  % 横軸にL2_valuesを使用
x_values_combined = L2_values + x_values_cos;  % L2_values と k2 * cos(a) の合計
% disp(x_values_combined);  % 度単位での値
% k2 * sin(a) の縦軸データ
y_values_sin = k2 * sin(a_values)+Rr1;

% グラフを作成
figure;
plot(x_values_combined, y_values_sin, '-', 'LineWidth', 3, 'DisplayName', 'k2 * sin(a)');
xlabel('X-coordinate at cource [mm]', 'FontSize', 12);
xlim([0 250]);
xticks(0:20:250);
axis equal;
ylabel('Y-coordinate at cource [mm]', 'FontSize', 12);
ylim([0 60]); % yの範囲
yticks(0:10:60); % 10刻みで目盛りを設定
grid off;
legend off;



