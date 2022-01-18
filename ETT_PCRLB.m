function [posCRLB] = ETT_PCRLB(runTime,MCRun,model,measProcess)
% Implementation of PCRLB recursion for target state (both kinematic and extent states) in 
% E. Sar?ta? and U. Orguner, “Posterior Cramer-Rao Lower Bounds for Extended Target Tracking with Random Matrices,” 
% in Proceedings of 19th International Conference on Information Fusion (FUSION'16), Jul. 2016.
% INPUTS:
% T: 1x1, sampling time
% runTime: 1x1 duration of the scenario
% MCRun: 1x1, number of Monte Carlo runs
% model: a structure describing the motion model with fields
%   -model.T: 1x1, sampling time
%   -model.F: 4x4, state transition model
%   -model.B: 4x2, the gain matrix for process noise
%   -model.Q: 2x2, covariance of the process noise
%   -model.H: 2x4, measurement model
%   -model.R: 2x2, covariance of the measurement noise
%   -model.s: 1x1, scaling factor
%   -model.n0: 1x1, initial degree of freedom used for the initilization
%                  of the posterior CRLB
%   -model.n: 1x1, degree of freedom 
%   -model.angle: 1x1, fixed orientation angle in degrees
%   -model.semiMajor: 1x1, extent semi-major axis length
%   -model.semiMinor: 1x1, extent semi-minor axis length
% measProcess: a structure describing measurement process with the field
%   -measProcess.mbar: 1x1, expected number of measurements
% OUTPUTS:
% posCRLB: a structure containing posterior CRLB's for kinematic and extent
%          states together with the bound on semi-major and semi-minor axes 


% Random extent generation to calculate expected values
randomExtension = zeros(2,2,runTime,MCRun);
for k = 1:MCRun   
    randomExtension(:,:,:,k) = generateRandomExtent(model,runTime);
end

% Posterior CRLB
posCRLB = posteriorCRLB(randomExtension,model,measProcess);

%%%

function randomExtension = generateRandomExtent(model,runTime)

n0 = model.n0;
n = model.n;

eigVec(1:2,1:2,1) = [cosd(model.angle) -sind(model.angle);...
                     sind(model.angle) cosd(model.angle)];
semiMin = model.semiMinor;
semiMaj = model.semiMajor;
initialExtension = eigVec(:,:,1)*[semiMaj^2 0; 0 semiMin^2]*eigVec(:,:,1)';     % (55)
randomExtension = zeros(2,2,runTime);
randomExtension(:,:,1) = wishrnd(initialExtension/n0,n0);      % (53b)
for t = 2:runTime    
    sigma = randomExtension(:,:,t-1)/n;
    randomExtension(:,:,t) = wishrnd(sigma,n);     % (57b)    
end

%%%

function [posCRLB] = posteriorCRLB(extension,model,measProcess)

[~, ~, runTime, ~] = size(extension);

n0 = model.n0;
n = model.n;

F = model.F;
Q = model.Q;
B = model.B;
H = model.H;
s = model.s;
R = model.R;
mbar = measProcess.mbar;

eigVec(1:2,1:2,1) = [cosd(model.angle) -sind(model.angle);...
    sind(model.angle) cosd(model.angle)];
semiMin = model.semiMinor;
semiMaj = model.semiMajor;
initialExtension = eigVec(:,:,1)*[semiMaj^2 0; 0 semiMin^2]*eigVec(:,:,1)';

FIMKin = zeros(4,4,runTime);
boundKin = zeros(4,4,runTime);

FIMExt = zeros(3,3,runTime);
boundExt = zeros(3,3,runTime);
boundEig = zeros(2,runTime);

initBound = diag([75 75 15 15]);    % (38)
FIMKin(:,:,1) = initBound^(-1);     % (37)
boundKin(:,:,1) = initBound;

FIMExt(:,:,1) = calculateInitialFIM(initialExtension,n0);
initBound = FIMExt(:,:,1)^(-1);
boundExt(:,:,1) = initBound;
[jacEig1,jacEig2] = jacobianEigPos(extension(:,:,1,:));
boundEig(1,1) = jacEig1'*boundExt(:,:,1)*jacEig1;
boundEig(2,1) = jacEig2'*boundExt(:,:,1)*jacEig2;

for t = 2:runTime-1
    
    % Kinematic
    FIMx = calculateFIMx(squeeze(extension(:,:,t+1,:)),mbar,H,s,R);
    FIMKin(:,:,t) = (B*Q*B'+F*FIMKin(:,:,t-1)^(-1)*F')^(-1)+FIMx;      % Equivalent of (36) using Matrix Inversion Lemma
    boundKin(:,:,t) = FIMKin(:,:,t)^(-1);
    
    % Extended
    [D11,D21,D221,FIMX] = calculateSubmatrices(squeeze(extension(:,:,t,:)),n,model,measProcess);
    D12 = D21';     % (43c)
    D22 = D221+FIMX;
    dummy = FIMExt(:,:,t-1)+D11;
    
    FIMExt(:,:,t) = D22-D21*dummy^(-1)*D12;     % (40b)
    boundExt(:,:,t) = (FIMExt(:,:,t))^(-1);     % (39)
    [jacEig1,jacEig2] = jacobianEigPos(extension(:,:,t,:));
    boundEig(1,t) = jacEig1'*boundExt(:,:,t)*jacEig1;
    boundEig(2,t) = jacEig2'*boundExt(:,:,t)*jacEig2;
end

posCRLB.kinematic = boundKin;
posCRLB.extension = boundExt;
posCRLB.eig = boundEig;

%%%

function [D11,D21,D221,FIMX] = calculateSubmatrices(X,n,model,measProcess)

sumMCD11_1 = zeros(2,2,2,2);
sumMCD11_2 = zeros(2,2,2,2);
sumMCD21_1 = zeros(2,2,2,2);
sumMCD21_2 = zeros(2,2,2,2);
sumMCD22_1 = zeros(2,2,2,2);
sumMCD22_2 = zeros(2,2,2,2);
sumMCD22_3 = zeros(2,2,2,2);
sumMCFIM_1 = zeros(2,2,2,2);
sumMCFIM_2 = zeros(2,2,2,2);

[~, d , MCRun] = size(X);

coefficientD11 = n/4;
coefficientD21 = -(n/4);

c2 = ((n-d)*(n-d-1)*(n-d-3))^(-1);       % (44b)
c1 = (n-d-2)*c2;       % (44a)
coefficientD22_1 = (n^2/4)*(c1*(n-d-1)^2-1);
coefficientD22_2 = c2*n^2*(n-d-1)^2/4;

R = model.R;
s = model.s;
mbar = measProcess.mbar;
coefficientFIM = s^2*mbar/4;

X11 = X(1,1,:);
X12 = X(1,2,:);
X21 = X(2,1,:);
X22 = X(2,2,:);

detX = X11.*X22-X12.*X21;
invX = zeros(d,d,MCRun);
invX(1,1,:) = X22./detX;
invX(1,2,:) = -X12./detX;
invX(2,1,:) = -X21./detX;
invX(2,2,:) = X11./detX;

invsXpR = zeros(d,d,MCRun);
Y11 = s*X11+R(1,1);
Y12 = s*X12+R(1,2);
Y21 = s*X21+R(2,1);
Y22 = s*X22+R(2,2);
detY = Y11.*Y22-Y12.*Y21;
invsXpR(1,1,:) = Y22./detY;
invsXpR(1,2,:) = -Y12./detY;
invsXpR(2,1,:) = -Y21./detY;
invsXpR(2,2,:) = Y11./detY;

for i = 1:2
    for j = 1:2
        for l = 1:2
            for m = 1:2
                sumMCD11_1(i,j,l,m) = sum(invX(i,l,:).*invX(j,m,:));  % Argument of the first expectation in (43a)
                sumMCD11_2(i,j,l,m) = sum(invX(i,m,:).*invX(l,j,:));  % Argument of the second expectation in (43a)
                
                sumMCD21_1(i,j,l,m) = sum(invX(i,l,:).*invX(j,m,:));  % Argument of the first expectation in (43b).  There is a typo in (43b). The first expectation should be E([X^-1]_{il}[X^-1]_{jm}) as written in this line.
                sumMCD21_2(i,j,l,m) = sum(invX(j,l,:).*invX(i,m,:));  % Argument of the second expectation in (43b)
                
                sumMCD22_1(i,j,l,m) = sum(invX(i,j,:).*invX(l,m,:));  % Argument of the first expectation in (43d)
                sumMCD22_2(i,j,l,m) = sum(invX(i,l,:).*invX(j,m,:));  % Argument of the second expectation in (43d)
                sumMCD22_3(i,j,l,m) = sum(invX(l,j,:).*invX(i,m,:));  % Argument of the third expectation in (43d)
                
                sumMCFIM_1(i,j,l,m) = sum(invsXpR(i,l,:).*invsXpR(j,m,:));
                sumMCFIM_2(i,j,l,m) = sum(invsXpR(i,m,:).*invsXpR(l,j,:));
            end
        end
    end
end

expectationD11_1 = sumMCD11_1/MCRun;
expectationD11_2 = sumMCD11_2/MCRun;
d11 = coefficientD11*(expectationD11_1+expectationD11_2);      % (43a)
D11 = [d11(1,1,1,1) d11(1,1,1,2)+d11(1,1,2,1) d11(1,1,2,2);      % (41)
    d11(1,2,1,1)+d11(2,1,1,1) d11(1,2,1,2)+d11(1,2,2,1)+d11(2,1,1,2)+d11(2,1,2,1) d11(1,2,2,2)+d11(2,1,2,2);
    d11(2,2,1,1) d11(2,2,1,2)+d11(2,2,2,1) d11(2,2,2,2)];

expectationD21_1 = sumMCD21_1/MCRun;
expectationD21_2 = sumMCD21_2/MCRun;
d21 = coefficientD21*(expectationD21_1+expectationD21_2);      % (43b)
D21 = [d21(1,1,1,1) d21(1,1,1,2)+d21(1,1,2,1) d21(1,1,2,2);      % (41)
    d21(1,2,1,1)+d21(2,1,1,1) d21(1,2,1,2)+d21(1,2,2,1)+d21(2,1,1,2)+d21(2,1,2,1) d21(1,2,2,2)+d21(2,1,2,2);
    d21(2,2,1,1) d21(2,2,1,2)+d21(2,2,2,1) d21(2,2,2,2)];

expectationD22_1 = sumMCD22_1/MCRun;
expectationD22_2 = sumMCD22_2/MCRun;
expectationD22_3 = sumMCD22_3/MCRun;
d22 = coefficientD22_1*expectationD22_1+coefficientD22_2*(expectationD22_2+expectationD22_3);      % (43d)
D221 = [d22(1,1,1,1) d22(1,1,1,2)+d22(1,1,2,1) d22(1,1,2,2);      % (41)
    d22(1,2,1,1)+d22(2,1,1,1) d22(1,2,1,2)+d22(1,2,2,1)+d22(2,1,1,2)+d22(2,1,2,1) d22(1,2,2,2)+d22(2,1,2,2);
    d22(2,2,1,1) d22(2,2,1,2)+d22(2,2,2,1) d22(2,2,2,2)];

expectationFIM_1 = sumMCFIM_1/MCRun;
expectationFIM_2 = sumMCFIM_2/MCRun;
fim = coefficientFIM*(expectationFIM_1+expectationFIM_2);
FIMX = [fim(1,1,1,1) fim(1,1,1,2)+fim(1,1,2,1) fim(1,1,2,2);      % (41)
    fim(1,2,1,1)+fim(2,1,1,1) fim(1,2,1,2)+fim(1,2,2,1)+fim(2,1,1,2)+fim(2,1,2,1) fim(1,2,2,2)+fim(2,1,2,2);
    fim(2,2,1,1) fim(2,2,1,2)+fim(2,2,2,1) fim(2,2,2,2)];

%%%

function [FIMx] = calculateFIMx(X,mbar,H,s,R)
% This function calculates the last expectation in (36).

[~, d, MCRun] = size(X);

Y11 = s*X(1,1,:)+R(1,1);
Y12 = s*X(1,2,:)+R(1,2);
Y21 = s*X(2,1,:)+R(2,1);
Y22 = s*X(2,2,:)+R(2,2);
detY = Y11.*Y22-Y12.*Y21;
EinvsXpR = zeros(d);
EinvsXpR(1,1) = sum(Y22./detY);
EinvsXpR(1,2) = -sum(Y12./detY);
EinvsXpR(2,1) = -sum(Y21./detY);
EinvsXpR(2,2) = sum(Y11./detY);
EinvsXpR = mbar*EinvsXpR/MCRun;

r1H = H(1,:);
r2H = H(2,:);

FIMx = EinvsXpR(1,1)*r1H'*r1H+EinvsXpR(1,2)*r1H'*r2H+EinvsXpR(2,1)*r2H'*r1H+EinvsXpR(2,2)*r2H'*r2H;

%%%

function result = calculateInitialFIM(Xbar,n)

X = Xbar;
d = length(X);
c2 = ((n-d)*(n-d-1)*(n-d-3))^(-1);   % (48b)
c1 = c2*(n-d-2);       %(48a)
coefficient1 = (n^2/4)*(c1*(n-d-1)^2-1);
coefficient2 = (n^2)*c2*(((n-d-1)^2/4));
invX = X^(-1);

d0 =  zeros(2,2,2,2);
for i = 1:2
    for j = 1:2
        for k = 1:2
            for l = 1:2
                d0(i,j,k,l) = coefficient1*invX(i,j)*invX(k,l)...      % (47)
                    +coefficient2*(invX(i,k)*invX(j,l)+invX(k,j)*invX(i,l));
            end
        end
    end
end
result = [d0(1,1,1,1) d0(1,1,1,2)+d0(1,1,2,1) d0(1,1,2,2);      % (46)
    d0(1,2,1,1)+d0(2,1,1,1) d0(1,2,1,2)+d0(1,2,2,1)+d0(2,1,1,2)+d0(2,1,2,1) d0(1,2,2,2)+d0(2,1,2,2);
    d0(2,2,1,1) d0(2,2,1,2)+d0(2,2,2,1) d0(2,2,2,2)];


%%%

function[jacEig1,jacEig2] = jacobianEigPos(X)

[~, ~, MCRun]=size(X);

X11 = squeeze(X(1,1,:))';
X12 = squeeze(X(1,2,:))';
X22 = squeeze(X(2,2,:))';
detX = X11.*X22-X12.^2;
trX = X11+X22;
sqrtDelta = sqrt(trX.^2-4*detX);

f1 = (trX+sqrtDelta)/2;
f2 = (trX-sqrtDelta)/2;

f1inv = 1./(2*sqrt(f1));
f2inv = 1./(2*sqrt(f2));

temp1 = (X11-X22)./(2*sqrtDelta);
temp2 = 2*X12./sqrtDelta;
temp3 = (X22-X11)./(2*sqrtDelta);

jacEig1 = sum([(0.5+temp1).*f1inv; temp2.*f1inv; (0.5+temp3).*f1inv],2)/MCRun;
jacEig2 = sum([(0.5-temp1).*f2inv; -temp2.*f2inv; (0.5-temp3).*f2inv],2)/MCRun;

