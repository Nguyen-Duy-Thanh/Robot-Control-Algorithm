clear all
close all

% Khoi tao thong so
dt = 0.05;                  % Thoi gian lay mau;
runTime = 40;               % Tong thoi gian chay
step = floor(runTime / dt);
v = 0.1;                    % Van toc ban dau
alpha = 0.2;                % Goc lai ban dau
d = sqrt(5);                % Chieu dai truc co so

% Landmarks
x_A = 3;
y_A = 4;

x_B = 4;
y_B = 4;

% Khoi tao ma tran hiep phuong sai
P(:, :, 1)=[0 0 0; 0 0 0; 0 0 0]; 

% Khoi tao vi tri
% Vi tri tinh toan dung Kalman
x(1) = 0;  
y(1) = 0;
theta(1) = 0;

% Vi tri that theo thong so trong contrlsig
x_actual(1) = 0;  
y_actual(1) = 0;
theta_actual(1) = 0;

% Vi tri du doan (khong dung Kalman)
x_predict(1) = 0;  
y_predict(1) = 0;
theta_predict(1) = 0;

% Phuong sai nhieu dau vao
Q = [2.5*10^-3 0; 
     0         3.6*10^-3];   
 
% Phuong sai nhieu do luong
R = [(10^-6) 0       0            0;
      0      (10^-6) 0            0;
      0      0       (7.62*10^-6) 0;
      0      0       0            (7.62*10^-6)];                 

% Phep do
z(:, 1) = [0; 0; 0; 0];

% Nhan thong so 
load('contrlsig.mat');
load('measdat25.mat');
v_actual = inp(1,:);     % Van toc do duoc
alpha_actual = inp(2,:); % Goc lai do duoc
rA_actual = me(1,:);     % Khoang cach toi diem A do duoc
rB_actual = me(2,:);     % Khoang cach toi diem B do duoc
thetaA_actual = me(3,:); % Goc toi diem A do duoc
thetaB_actual = me(4,:); % Goc toi diem B do duoc


for k = 2 : step
    % Vi tri that cua robot
    x_actual(k) = x_actual(k-1) + v_actual(k-1)*cos(theta_actual(k-1))*dt;
    y_actual(k) = y_actual(k-1) + v_actual(k-1)*sin(theta_actual(k-1))*dt;
    theta_actual(k) = theta_actual(k-1) +v*tan(alpha_actual(k-1))/d*dt;
    
    
    % Du doan dua tren mo hinh he thong (khong su dung bo loc Kalman)
    x_predict(k) = x_predict(k-1) + v*cos(theta_predict(k-1))*dt;
    y_predict(k) = y_predict(k-1) + v*sin(theta_predict(k-1))*dt;
    theta_predict(k) = theta_predict(k-1) + v*tan(alpha)/d*dt;

    % Du doan dung bo loc Kalman
    % Pha du doan
    x(k) = x(k-1) + v*cos(theta(k-1))*dt;
    y(k) = y(k-1) + v*sin(theta(k-1))*dt;
    theta(k) = theta(k-1) + v*tan(alpha)/d*dt;
    
    A=[1 0 -v*sin(theta(k-1))*dt;
       0 1 v*cos(theta(k-1))*dt;
       0 0 1];
   
    W = [cos(theta(k-1))*dt 0;
         sin(theta(k-1))*dt 0;
         tan(alpha)*dt/d    v*dt/(d*(cos(alpha)^2))];

    H=[(x(k-1)-x_A)/sqrt((x(k-1)-x_A)^2 + (y(k-1)-y_A)^2) (y(k-1)-y_A)/sqrt((x(k-1)-x_A)^2 + (y(k-1)-y_A)^2) 0;
       (x(k-1)-x_B)/sqrt((x(k-1)-x_B)^2 + (y(k-1)-y_B)^2) (y(k-1)-y_B)/sqrt((x(k-1)-x_B)^2 + (y(k-1)-y_B)^2) 0;
       (y_A-x(k-1))/((x(k-1)-x_A)^2 + (y(k-1)-y_A)^2)     (x(k-1)-x_A)/((x(k-1)-x_A)^2 + (y(k-1)-y_A)^2)     -1;
       (y_B-x(k-1))/((x(k-1)-x_B)^2 + (y(k-1)-y_B)^2)     (x(k-1)-x_B)/((x(k-1)-x_B)^2 + (y(k-1)-y_B)^2)     -1];
    
    Pe(:,:,k)=A*P(:,:,k-1)*A' + W*Q*W';        % Cap nhat ma tran hiep phuong sai 
    K(:,:,k)=Pe(:,:,k)*H'/(H*Pe(:,:,k)*H'+R);  % hHe so Kalman
    
    % Du doan thong so khoang cach
    rA = sqrt((x(k)-x_A)^2 + (y(k)-y_A)^2);
    rB = sqrt((x(k)-x_B)^2 + (y(k)-y_B)^2);
    betaA = atan((y_A-y)/(x_A-x)) - theta(k);
    betaB = atan((y_B-y)/(x_B-x)) - theta(k);
    z_pre(:,k) = [rA; rB; betaA; betaB];

    % Thong so khoang cach that do duoc
    rA_a = sqrt((x_actual(k)-x_A)^2 + (y_actual(k)-y_A)^2);
    rB_b = sqrt((x_actual(k)-x_B)^2 + (y_actual(k)-y_B)^2);
    betaA_a = atan((y_A-y_actual)/(x_A-x_actual)) - theta_actual(k);
    betaB_b = atan((y_B-y_actual)/(x_B-x_actual)) - theta_actual(k);
    z(:,k) = [rA_a; rB_b; betaA_a; betaB_b];
    
    % Vi tri du doan
    X_pre(:,k) = [x(k);y(k);theta(k)];  
    
    % Pha cap nhat
    r(:,k)=z(:,k)-z_pre(:,k); % Inovation
    X = X_pre(:,k) + K(:,:,k)*r(:,k);  % Uoc tinh vi tri dung Kalman 
    x(k) = X(1);
    y(k) = X(2);
    z(k) = X(3);

    P(:,:,k)=(eye(3)-K(:,:,k)*H)*Pe(:,:,k); % Uoc tinh hiep phuong sai
end

% ve do thi ket qua
% do thi duong di trong 3 truong hop: that, du doan trong 2TH: theo mo hinh he thong va su dung bo loc Kalman
time=0:step-1;
figure(1)
plot(x,y,x_predict,y_predict,x_actual,y_actual,'--')
xlabel('X (m)');
ylabel('Y (m)');
legend('Kalman','Prediction','Actual','Location','SouthEast' )
title('Actual and Estimated Trajectories')

% do thi the hien sai so theo truc x khi du doan so voi ket qua that (ca 2 TH)
figure(2)
plot(time,x-x_actual,time,x_predict-x_actual,'--')
legend('Kalman error','Prediction error')
xlabel('Time (step)');
ylabel('Error (m)');
title('Error in X')
grid

% do thi the hien sai so theo truc y khi du doan so voi ket qua that (ca 2 TH)
figure(3)
plot(time,y-y_actual,time,y_predict-y_actual,'--')
xlabel('Time (step)');
ylabel('Error (m)');
legend('Kalman error','Prediction error')
title('Error in Y')
grid

% do thi the hien sai so theo goc di chuyen cua xe khi du doan so voi ket qua that (ca 2 TH)
figure(4)
plot(time,theta-theta_actual,time,theta_predict-theta_actual,'--')
xlabel('Time (step)');
ylabel('Error (rad)');
legend('Kalman error','Prediction error')
title('Error in Theta')
grid