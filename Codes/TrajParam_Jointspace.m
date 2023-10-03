% Calculates parameterized symbolic function of the recorded trajectory 'Q_dim.mat'.
%
% The recorded trajectory is a [DoF+1, N] array, the first row
% is the recorded time (not used), N is the number of datapoints.
%
%
% Written by Giovanni Braglia and Davide Tebaldi, 2023
% University of Modena and Reggio Emilia
% website: https://www.automatica.unimore.it/

clear all
close all
clc


DoF = 7; % Degree of Freedom, 7 for Franka robot
load('Q_dim.mat');


qtr = Q_dim(2:8,:);
qtr = qtr';

l  = length(qtr);
Ts = 0.001;
time_qtr = linspace( 0, l*Ts, l )';



%%

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Application of Spatial Sampling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[tn,sn,xn]= SpatialSampling(time_qtr,qtr,0.006);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 1:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NrFig=1;
Fig_Name='Figure_1';
hh=figure(NrFig); clf;
set(hh,'NumberTitle','off')
set(hh,'Name',[num2str(NrFig) '.' Fig_Name])
% set(hh,'PaperPositionMode','manual','PaperPosition',[5 0 15 17.75]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii=1:size(qtr,2)
    subplot(4,2,ii)
    plot(time_qtr,qtr(:,ii))
    grid on; zoom on; hold on;
    plot(tn,xn(:,ii),'r--')
    xlabel('Time [s]')
    ylabel('rad')
    title(['q' num2str(ii)])
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 2:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NrFig=2;
Fig_Name='Figure_2';
hh=figure(NrFig); clf;
set(hh,'NumberTitle','off')
set(hh,'Name',[num2str(NrFig) '.' Fig_Name])
% set(hh,'PaperPositionMode','manual','PaperPosition',[5 0 15 17.75]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% s_range=sn;
s_range = linspace( 0, sn(end), length(sn) );
s_fin=s_range(end);
time_qtr=tn;
ds=diff(s_range)/(time_qtr(2)-time_qtr(1));
ds2=ds.^2;
dds=diff(ds)/(time_qtr(2)-time_qtr(1));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameterization with Radial Basis Functions (RBF)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 22;  % Number of RBF
L = length( xn );
c = linspace ( 0 - 1/(N-3), s_fin + 1/(N-3), N )' ; % centers of RBF
h = 1/(2*(N));  % variance of RBF


RBF   = cell([N,1]);   % array of RBF functions
dRBF  = cell([N,1]);   % array of first derivatives
ddRBF = cell([N,1]);   % array of second derivatives

for jj=1:N
    RBF{jj}   = @(s) exp(-(s-c(jj)).^2/(2*h));
    dRBF{jj}  = @(s) -(s-c(jj))/h.*RBF{jj}(s);
    ddRBF{jj} = @(s) -1/h.*RBF{jj}(s)+(s-c(jj)).^2/(h.^2).*RBF{jj}(s);
end

sum_rbf   = @(s)RBF{1}(s);
dsum_rbf  = @(s)dRBF{1}(s);
ddsum_rbf = @(s)ddRBF{1}(s);
for jj=2:N
    sum_rbf   = @(s) sum_rbf(s)+RBF{jj}(s);
    dsum_rbf  = @(s) dsum_rbf(s)+dRBF{jj}(s);
    ddsum_rbf = @(s) ddsum_rbf(s)+ddRBF{jj}(s);
end


Phi   = cell([N,1]); % normalized RBF
dPhi  = cell([N,1]); % normalized dRBF
ddPhi = cell([N,1]); % normalized ddRBF

for jj=1:N
    Phi{jj}   = @(s) RBF{jj}(s)./sum_rbf(s);
    dPhi{jj}  = first_der(RBF{jj},dRBF{jj},sum_rbf,dsum_rbf);
    ddPhi{jj} = second_der(RBF{jj},dRBF{jj},ddRBF{jj},sum_rbf,dsum_rbf,ddsum_rbf);
end

for ii=1:size(qtr,2)

    qd=xn(:,ii);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Optimization constraints
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PSI=ones([L,N]);
    dPSI=ones([L,N]);
    ddPSI=ones([L,N]);

    Aeq=zeros([4,N]);  % constraint matrix for initial/final
    beq=zeros([4,1]);  % velocities and accelerations --> NOT USED

    Ayg=ones([2,N]);   % constraint matrix for initial/final positions
    byg=[qd(1),qd(end)]';

    for jj=1:N
        PSI(:,jj)=Phi{jj}(s_range);
        dPSI(:,jj)=dPhi{jj}(s_range);
        ddPSI(:,jj)=ddPhi{jj}(s_range);
        %
        Ayg(1,jj)=Phi{jj}(0);
        Ayg(2,jj)=Phi{jj}(s_fin);
    end


    dqd  = diff(qd)/(Ts);
    ddqd = diff(dqd)/(Ts);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Weights Calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    W = LagrangeApprox( qd, PSI, beq, Aeq, Ayg, byg, N );
    clc;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % create symbolic functions for Simulink ( NORMALIZED VERSION )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  
%   | ---------------------------------------- |
%   | The followings creates symbolic          |
%   | 'file.m' to be used in Simulink or       |
%   | to test  'OptimalPhase_Jointspace.m'.    |
%   | ---------------------------------------- |
%
    q   = @(s) W(1).*Phi{1}(s);
    dq  = @(s) W(1).*dPhi{1}(s);
    ddq = @(s) W(1).*ddPhi{1}(s);

    for jj = 2:N
        q   = @(s) q(s)   + W(jj).*Phi{jj}(s);
        dq  = @(s) dq(s)  + W(jj).*dPhi{jj}(s);
        ddq = @(s) ddq(s) + W(jj).*ddPhi{jj}(s);
    end

    ysym   = q(sym('s'));


    matlabFunction(ysym,'File',['q' num2str(ii) '_of_s']);

    dysym  = dq(sym('s'));
    matlabFunction(dysym,'File',['dq' num2str(ii) '_of_s']);

    ddysym = ddq(sym('s'));
    matlabFunction(ddysym,'File',['ddq' num2str(ii) '_of_s']);

    clear q dq ddq ysym dysym ddysym


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    qtr_par(:,ii)=PSI*W;
    dqtr_par(:,ii)=dPSI*W;
    ddqtr_par(:,ii)=ddPSI*W;

    qtrii_par=qd*0;
    dqtrii_par=qd*0;
    ddqtrii_par=qd*0;
    for jj=1:length(qtrii_par)
        eval(['qtrii_par(jj)=q' num2str(ii) '_of_s(s_range(jj));']);
        eval(['dqtrii_par(jj)=dq' num2str(ii) '_of_s(s_range(jj));']);
        eval(['ddqtrii_par(jj)=ddq' num2str(ii) '_of_s(s_range(jj));']);
    end


    figure(2)
    %
    subplot(3,7,ii);
    plot(s_range,qd,'r--');
    hold on; grid on; zoom on;
    plot(s_range,qtrii_par,'b');
    ylabel('[rad]')
    title('actual (r) and par (b)')
    %
    subplot(3,7,ii+7);
    hold on; grid on; zoom on;
    plot(s_range(2:end),dqtrii_par(2:end).*ds,'b');
    ylabel('[rad/s]')
    title('calc (g) and par (b)')
    %
    subplot(3,7,ii+14);
    hold on; grid on; zoom on;
    plot(s_range(3:end),ddqtrii_par(3:end).*ds2(2:end)+dqtrii_par(3:end).*dds,'b');
    xlabel('Time [s]')
    ylabel('[rad/s^2]')
    title('calc (g) and par (b)')


end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save variables for Optimization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save('Filtered_Trajectory_Q.mat','tn','sn','xn')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Auxiliary Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%---------------------------------------------------------------------
function [tn,sn,xn]= SpatialSampling(t,x,delta)
%---------------------------------------------------------------------
% Returns equally spaced ('delta') points on trajectory 'x'

    [nrow,ncol]=size(x);
    if ncol>nrow
        x=x';
    end
    %
    for jj=1:size(x,2)
        eval(['xn' num2str(jj) '=[x(1,' num2str(jj) ')];']);
    end
    %
    sn=[t(1)];
    tn=[t(1)];
    %
    xcurr=x(1,:);
    for ii=1:size(x,1)-1
        %
        if norm(x(ii+1,:)-xcurr)>delta
            nr_steps=floor(norm(x(ii+1,:)-xcurr)/delta);
            for pp=1:nr_steps
                sn=[sn; sn(end)+delta];
                for jj=1:size(x,2)
                    xincr=xcurr(jj)+delta*(x(ii+1,jj)-xcurr(jj))/norm(x(ii+1,:)-xcurr);
                    eval(['xn' num2str(jj) '=[xn' num2str(jj) '; ' num2str(xincr) '];']);
                end
                xcurr=x(1,:)*0;
                for jj=1:size(x,2)
                    eval(['xcurr(' num2str(jj) ')=xn' num2str(jj) '(end);']);
                end
                nDX=norm(x(ii+1,:)-x(ii,:));
                DT=t(ii+1)-t(ii);
                K=nDX/DT;
                tn=[tn; (norm(xcurr-x(ii,:))+K*t(ii))/K];
            end
        end
    end
    xn=zeros(size(xn1,1),size(x,2));
    for jj=1:size(x,2)
        eval(['xn(:,' num2str(jj) ')=xn' num2str(jj) ';']);
    end

end

%---------------------------------------------------------------------
function dphi = first_der( rbf, drbf, sum_rbf, dsum_rbf )
%---------------------------------------------------------------------
% Calculates analitically the first derivative of normalized RBF

    num = @(s) drbf(s) .* sum_rbf(s) - rbf(s) .* dsum_rbf(s);
    den = @(s) sum_rbf(s).^2;

    dphi = @(s) num(s) ./ den(s);

end


%---------------------------------------------------------------------
function ddphi = second_der( rbf, drbf, ddrbf, sum_rbf, dsum_rbf, ddsum_rbf )
%---------------------------------------------------------------------
% Calculates analitically the second derivative of normalized RBF


    a = @(s) ddrbf(s) .* sum_rbf(s) - rbf(s) .* ddsum_rbf(s);
    b = @(s) drbf(s) .* sum_rbf(s) - rbf(s) .* dsum_rbf(s);


    num = @(s) a(s) .* ( sum_rbf(s).^2 ) -  b(s) .* ( 2 * sum_rbf(s) ) .* dsum_rbf(s);
    den = @(s) ( sum_rbf(s).^4 );

    ddphi = @(s) num(s) ./ den(s);

end


%---------------------------------------------------------------------
function P = LagrangeApprox( yd, PSI, beq, Aeq, Ayg, byg, N )
%---------------------------------------------------------------------
% Computes analitically the weighting coefficients for the
% approximation of 'yd' with RBF basis functions

    L  = length( yd );
    W  = eye( L ) .* 1e0;
    Bi = vertcat( Aeq, Ayg );
    Ri = vertcat( beq, byg );
    Ba = PSI;
    Ra = yd;

    Ba2 = Ba' * W * Ba + eye(N).*realmin;
    BRa = Ba' * W * Ra;

    [ k, ~ ] = size( Bi );
    Q = ( Bi / Ba2 ) * Bi'  + eye(k).*realmin ;

    l = Q \ ( ( Bi / Ba2 ) * BRa - Ri );
    P = ( Ba2 \ BRa ) - ( Ba2 \ Bi' ) * l ;

end
