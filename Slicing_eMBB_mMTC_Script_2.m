%% Description
% This script generates the curves of eMBB data rate r_B versus the maximum
% number of connected MTC devices for a given mMTC data rate r_M.

%% Reset
clearvars
close all
clc

%% Parameters
N=1e4;                                      % Number of Monte Carlo runs
Gamma_B_dB=20;                              % Average received SNR of the eMBB device [dB]
Gamma_B=10^(Gamma_B_dB/10);                 % Average received SNR of the eMBB device
Gamma_m_dB=5;                               % Average received SNR of the mMTC device [dB]
Gamma_m=10^(Gamma_m_dB/10);                 % Average received SNR of the mMTC device
Em=1e-1;                                    % Reliability requirement for mMTC
Eb=1e-3;                                    % Reliability requirement for eMBB
alpha=0:0.1:1;                              % Time sharing factor for the orthogonal slicing
Pm=1;                                       % Normalized transmit power of the mMTC device
M_max=100;                                    % Maximum number of mMTC devices connected to the BS
L=8;                                        % Number of receiver antennas
rm=0.25;                                     % mMTC data rate [bits/s/Hz]
rng('default')                              % Set the random number generator

%% eMBB parameters
gamma_B_min=Gamma_B*gammaincinv(Eb*factorial(L-1)/gamma(L),L,'lower');      % Threshold SNR
gamma_B_tar_max=Gamma_B*factorial(L-1)/igamma(L-1,gamma_B_min/Gamma_B);     % Maximum Target SNR
rb_max=log2(1+gamma_B_tar_max);                                             % Outage rate of the eMBB user
rb_ort=alpha*rb_max;                                                        % eMBB data rate [bits/s/Hz] - Orthogonal Slicing
rb_non=[0:1:rb_max,round(rb_max+0.05,1)];                                   % eMBB data rate [bits/s/Hz] - Non-Orthogonal Slicing

%% Channel realizations for the mMTC users
Hm_max=(1/sqrt(2))*(randn(L,M_max,N)+1i*randn(L,M_max,N));      % Channel coefficients of all the mMTC Devices
Gm_max=sqrt(Gamma_m)*Hm_max;                                    % Channel gains of all the mMTC Devices

%% Channel realizations for the eMBB user
Hb=(1/sqrt(2))*(randn(L,1,N)+1i*randn(L,1,N));              % Channel coefficients of the eMBB device                                          
Gb=sqrt(Gamma_B)*Hb;                                        % Channel gains of the eMBB device
gamma_B=squeeze(Gamma_B*vecnorm(Hb,2,1).^2);                % Received SNR of the eMBB Device

%% Monte Carlo Simulation - Orthogonal Slicing between eMBB and mMTC
M_ort=zeros(1,length(alpha));
parfor k=1:length(alpha)
    for M=1:M_max                                        % Number of mMTC devices connected to the BS      
        Dm=zeros(1,N);                                      % Number of correctly decoded mMTC packets in each timeslot
        for j=1:N
            Gm=Gm_max(:,1:M,j);                               % Select the channel gains of the mMTC devices
            [~,ord]=sort(vecnorm(Gm,2,1).^2,'descend');
            Gm_ord=Gm(:,ord);                               % Sort channel gains of the mMTC devices in descending order
            for m0=1:size(Gm,2)
                if ~isempty(Gm)                             % If there is at least one mMTC user in the timeslot:
                    I=0;
                    for m1=m0+1:size(Gm,2)
                        I=I+abs(Gm_ord(:,m0)'*Gm_ord(:,m1))^2;       % Interference from other active mMTC devices
                    end
                    Sigma_m=(Pm*vecnorm(Gm_ord(:,m0),2,1)^4)/(Pm*I+vecnorm(Gm_ord(:,m0),2,1)^2);    % SINR of the mMTC device being decoded                                        
                   if log2(1+Sigma_m)>=rm/(1-alpha(k))
                      Dm(j)=m0;                             % Update the number of correctly decoded mMTC packets
                   else
                       break;
                   end
               end
            end
        end          
        Pr_Em=1-mean(Dm)/M;                      % Error probability of mMTC
        if Pr_Em<=Em                             % Test if the error probability satisfy the reliability requirement
           M_ort(k)=M;
        else
           break
        end
    end
end

%% Monte Carlo Simulation - Non-Orthogonal Slicing between eMBB and mMTC
M_non=zeros(1,length(rb_non));
parfor x=1:length(rb_non)
    Gm_max_temp=Gm_max;
    gamma_B_temp=gamma_B;
    Gb_temp=Gb;
    gamma_B_min=(2^rb_non(x))-1;
    for M=1:M_max                             % Number of MTC devices connected to the BS        
        flag=0;
        for gamma_B_tar=gamma_B_min:0.5:gamma_B_tar_max 
            Db=zeros(1,N);                              % Variable that indicates if the eMBB device has already been decoded
            Dm=zeros(1,N);                              % Number of correctly decoded mMTC packets in each timeslot
            for j=1:N
                Gm=Gm_max_temp(:,1:M,j);                  % Select the channel gains of the mMTC devices
                [~,ord]=sort(vecnorm(Gm,2,1).^2,'descend');
                Gm_ord=Gm(:,ord);                               % Sort channel gains of the mMTC devices in descending order
                Pb=gamma_B_tar/gamma_B_temp(j);             % Transmission power of the eMBB device
                if ~isempty(Gm)                             % If there is at least one mMTC user in the timeslot:
                    for m0=1:size(Gm,2)
                        if Db(j)==0                         % If the eMBB device has not been decoded yet:                                
                            I=0;
                            for m1=m0+1:size(Gm,2)
                                I=I+abs(Gm_ord(:,m0)'*Gm_ord(:,m1))^2;       % Interference from other active mMTC devices
                            end
                            Sigma_m=(Pm*vecnorm(Gm_ord(:,m0),2,1)^4)/(Pm*I+Pb*(abs(Gm_ord(:,m0)'*Gb_temp(:,1,j))^2)+vecnorm(Gm_ord(:,m0),2,1)^2);    % SINR of the mMTC device being decoded                                                                                           
                            if log2(1+Sigma_m)>=rm
                                Dm(j)=m0;                   % Update the number of correctly decoded mMTC packets
                            else
                                I=0;
                                for m1=m0:size(Gm,2)
                                    I=I+abs(Gb_temp(:,1,j)'*Gm_ord(:,m1))^2;       % Interference from other active mMTC devices
                                end                                                                        
                                Sigma_b=(Pb*vecnorm(Gb_temp(:,1,j),2)^4)/(Pm*I+vecnorm(Gb_temp(:,1,j),2,1)^2);      % SINR of the eMBB device that is being decoded                     
                                if log2(1+Sigma_b)>=rb_non(x)
                                    Db(j)=1;
                                    Sigma_m=(Pm*vecnorm(Gm_ord(:,m0),2,1)^4)/(Pm*I+vecnorm(Gm_ord(:,m0),2,1)^2);    % SINR of the mMTC device being decoded (second attempt, without the interference from eMBB)                                                                                          
                                    if log2(1+Sigma_m)>=rm
                                        Dm(j)=m0;                   % Update the number of correctly decoded mMTC packets
                                    else
                                       break; 
                                    end                                        
                                else
                                    break;
                                end
                            end                            
                        else                                % Else, if the eMBB device has already been decoded:                                
                            I=0;
                            for m1=m0+1:size(Gm,2)                                
                                I=I+abs(Gm_ord(:,m0)'*Gm_ord(:,m1))^2;       % Interference from other active mMTC devices                                
                            end
                            Sigma_m=(Pm*vecnorm(Gm_ord(:,m0),2,1)^4)/(Pm*I+vecnorm(Gm_ord(:,m0),2,1)^2);    % SINR of the mMTC device being decoded                                                                                                                                                                                 
                            if log2(1+Sigma_m)>=rm
                                Dm(j)=m0;                   % Update the number of corrected decoded mMTC packets
                            else
                                break;
                            end
                        end
                    end
                    % After the correct decoding of all active mMTC users, if the eMBB device has not been decoded yet:
                    if m0==size(Gm,2) && Db(j)==0
                        Sigma_b=gamma_B_tar;                 % SNR of the eMBB device being decoded
                        if log2(1+Sigma_b)>=rb_non(x)
                           Db(j)=1;
                        end
                    end
                else                                % In the case of zero active mMTC users in the minislot:
                    Sigma_b=gamma_B_tar;            % SNR of the eMBB device being decoded
                    if log2(1+Sigma_b)>=rb_non(x)
                        Db(j)=1;
                    end
                end                
            end
            Pr_Em=1-mean(Dm)/M;         % Error Probability of mMTC
            Pr_Eb=1-mean(Db);           % Error Probability of eMBB
            if Pr_Em<=Em                % Test if the error probability satisfy the reliability requirement
                if Pr_Eb<=Eb            % Test if the error probability satisfy the reliability requirement
                    M_non(x)=M;
                    flag=1;
                    break; 
                end                
            else
                break;
            end
        end
        if flag==0
            break
        end
    end
end

%% Saving the results
rb_ort_8=rb_ort;
M_ort_8=M_ort;
rb_non_8=rb_non;
M_non_8=M_non;
save('Results_8.mat','rb_ort_8','M_ort_8','rb_non_8','M_non_8')

%% Plotting the curves
fig1=figure(1);
    set(fig1,'Position',[200 200 650 500])
    hold on
    plot(rb_ort_8,M_ort_8,'k--','LineWidth',1.5)
    plot(rb_non_8,M_non_8,'k-','LineWidth',1.5)    
    grid on
    set(gca,'TickLabelInterpreter','latex','FontSize',12)
    xlabel('$r_B$ [bits/s/Hz]','Interpreter','latex','fontsize',16)
    ylabel('$r_M$ [bits/s/Hz]','Interpreter','latex','fontsize',16)
    leg=legend('Orthogonal','Non-Orthogonal');
    set(leg,'Interpreter','latex','fontsize',12)
    