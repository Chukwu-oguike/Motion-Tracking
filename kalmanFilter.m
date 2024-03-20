
%Function filters gravity components and magnetometer components
%gData are the components of gravity, mData are the components of the earth's magnetic field,
%wData are the angular velocity components and dw is the derivative of the angular velocity components

%This function returns the filtered gravity components, magnetometer
%components and the kalman gain
function [gFiltered, mFiltered,measResidual,kGain]=kalmanFilter(gData, mData, wData, dwData, Cov,version)


% covariance matrices
Q1 = cell(1,length(wData));
for i=1:length(wData)
    newcov_w = Cov.w;
    Q1{i} = [newcov_w zeros(3,3); zeros(3,3) newcov_w];%error Covariance for rate-gyro data
end


Q2 = [Cov.proc zeros(3,3); zeros(3,3) Cov.proc];
R = [Cov.g zeros(3,3); zeros(3,3) Cov.m]; %Measurement covariance

%in_out data sample time
dt = 0.01;

% measuremente model
H = eye(6);

%initialization
thlen = length(gData);
x_kk_k = zeros(6,thlen);%a_priori state vector
x_kk_kk = zeros(6,thlen);%a_posteriori state vector
measResidual = zeros(6,thlen);
Ktrace = zeros(6,thlen);
P_kk_k = cell(1,thlen);%a_priori Covariance Matrix 
P_kk_kk = cell(1,thlen);%a_posteriori
kGain = cell(1,thlen);%kalman gain

x_kk_kk(:,1) = [gData(:,1); mData(:,1)];
P_kk_kk{1} = Q1{1}+Q2;
kGain{1} = zeros(size(x_kk_kk,1),6);
RotMat = eye(3);

for j=2:thlen
    
    %measurements
    y = [gData(:,j); mData(:,j)];
    
    %sketwmat creates a skewmatrix with the components of the vectors
    Sw = skewmat(wData(:,j));
    Sdw = skewmat(dwData(:,j));
    Sxg = skewmat(x_kk_kk(1:3,j-1));
    Sxm = skewmat(x_kk_kk(4:6,j-1));
    
    %F is the process matrix and W is zero-mean white noise 
    if strcmp(version,'original')
        F =  [-Sw zeros(3,3); zeros(3,3) -Sw]*dt + eye(6);
        W = eye(6);
        
    elseif strcmp(version,'firstorder')
        %first order approx
        F =  [-Sw zeros(3,3); zeros(3,3) -Sw]*dt + eye(6);
        W = [-Sxg zeros(3,3); zeros(3,3) -Sxm]*dt;

    elseif strcmp(version,'secondorder')
        %second order approx
        
        
        F1 = (-Sw*dt-0.5*Sdw*dt^2+0.5*Sw*Sw*dt^2);
        F = [F1 zeros(3,3); zeros(3,3) F1] + eye(6);
        
        
        %neglecting delta in W (the complete expression includes terms in S(delta) and S(S(delta)*x)
        Wg = (Sw*dt^2-eye(3)*dt)*Sxg;
        Wm = (Sw*dt^2-eye(3)*dt)*Sxm;

        
        
        W = [Wg zeros(3,3); zeros(3,3) Wm];
              
    end

    %prediction Step
    x_kk_k(:,j) = F*x_kk_kk(:,j-1);
    P_kk_k{j} = F*P_kk_kk{j-1}*F' + W*Q1{j}*W' + Q2;

    %update step
    kGain{j} = P_kk_k{j}*H' / (H*P_kk_k{j}*H'+R);
    measResidual(:,j) = (y-H*x_kk_k(:,j));
    x_kk_kk(:,j) = x_kk_k(:,j) + kGain{j}*measResidual(:,j);
    P_kk_kk{j} = (eye(length(x_kk_kk(:,1))) - kGain{j}*H) * P_kk_k{j};
    

end

g_out = x_kk_kk(1:3,:);
m_out = x_kk_kk(4:6,:);