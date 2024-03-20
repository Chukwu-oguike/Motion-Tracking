

%v is the body frame and w is the inertial reference frame for the vectors
%g and m. (g and m represent the components of the earth's gravitational
%and magnetic fields respectively)

%This function returns the optimal quaternion representing the rotation between two reference frames
%, and the rotation matrix between two reference frames
function [q_opt,R,res,q] = questalgorithm(v,w,dataStandDev)

if nargin<1
v.g = [1 1 0]';
v.m = [0 -1 0]';
w.g = [0 1 0]';
w.m = [1 -1 0]';
end

%convert to unit vectors
v_g =v.g/norm(v.g);
v_m =v.m/norm(v.m);
w_g =w.g/norm(w.g);
w_m =w.m/norm(w.m);



if nargin<3
    dataStandDev.v_g = 0.1;
    dataStandDev.v_m = 0.1;
    dataStandDev.w_g = 0.1;
    dataStandDev.w_m = 0.1;
end

%sum of variance parameters in two reference frames
totDatVar_g = dataStandDev.v_g^2 + dataStandDev.w_g^2;
totDatVar_m = dataStandDev.v_m^2 + dataStandDev.w_m^2;

%total variance 
sq_sig_tot = 1/(1/totDatVar_g+1/totDatVar_m);

%weights for vectors (Wahba's problem)
a_g = sq_sig_tot/totDatVar_g;
a_m = sq_sig_tot/totDatVar_m;


%Attitude  profile Matrix 
B = a_g*w_g*v_g' + a_m*w_m*v_m';

%parameters of the K matrix
sLowerCase = trace(B);
S = B+B';
Z = a_g*cross(w_g,v_g) + a_m*cross(w_m,v_m);

K = [S-sLowerCase*eye(3), Z; Z', sLowerCase];


%Find eigen value for K matrix
try
    [q,lambda] = eig(K);
catch me
    q_opt = [0 0 0 0]';
    R = zeros(3,3);
    res = zeros(6,1);
    q = zeros(4,4);
    return
end

%largest eigen value optimizes the gain function
[remove_a,eig_idx] = max(max(lambda));

q_opt = q(:,eig_idx);

q_i = q_opt(1:3);
q_r = q_opt(4);
R = (q_r^2-q_i'*q_i)*eye(3)+2*(q_i*q_i')+2*q_r*[0 q_i(3) -q_i(2); -q_i(3) 0 q_i(1); q_i(2) -q_i(1) 0];
res.g = w_g-R*v_g;
res.m = w_m-R*v_m;

R = cell(4);

%Check to ensure that max eigen value maximizes gain function
for i=1:4
q_opt = q(:,i);
q_i = q_opt(1:3);
q_r = q_opt(4);
R{i} = (q_r^2-q_i'*q_i)*eye(3)+2*(q_i*q_i')+2*q_r*[0 q_i(3) -q_i(2); -q_i(3) 0 q_i(1); q_i(2) -q_i(1) 0];
res.g(:,i) = w_g-R{i}*v_g;
res.m(:,i) = w_m-R{i}*v_m;
myres(i) = norm([res.g(:,i); res.m(:,i)]);
end

[remove_b,bestidx] = min(myres);

q_opt = q(:,bestidx);
R = R{bestidx};
res.g = res.g(:,bestidx);
res.m = res.m(:,bestidx);

    
for i=1:4
    if q(4,i)<-0.01
       q(:,i)=-q(:,i);
    end
end
      
    
for i=1:size(q_opt,2)
    if q_opt(4,i)<-0.01
           q_opt=-q_opt;
    end
end
