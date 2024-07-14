%% README
%
% INTRO
% This script verifies that the postprocessing procedure, used for
% postprocessing the measurements using the impdance tube, is correct. 
% by R. Liupekevicius 20-12-2021. Eindhoven University of Technology
%
% OBJECTIVE
% Apply two-load method for one frequency/ freq band. A time signal is
% generated via the analytical solution of a forced 1D impedance tube. If
% you'd like to see how coefficients are computed, ask for 
% analytical_T_and_R.nb and it_analytical_sol.nb Mathematica files.
% 
%
%
% INSTRUCTIONS
% Main settings are commented by a sandwitch of '%%%%' together with a
% short description. Please see the example below.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% EXAMPLE of settings to change %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parameter_to_change = 1; % description about 'parameter_to_chage'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% CLEAN WORKSPACE/COMMAND WINDOW AND CLOSE FIGURES
clear all   % clear all variables in the workspace
close all   % close figures. Comment if multiple samples need to be plotted
clc         % clears the command window

%% ADD FUNCTIONS TO PATH
% % add the shared functions to the search path
% % path(path,...
% %     ['C:\GitLab\Tools\impedance_tube_post_processing\shared_functions']);    

%% MEASUREMENT SPECS
% Measurement specifications

% sample thickness is not used directly on the calculation of the
% coefficients however it highly influences them. If the sample's surface
% impedance would like to be computed, the sample thickness is a must.
d         = 0.01;       % sample thickness [m]
temp      = 293.15;     % temperature during measurement [K]



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency range: planewave up to 1.84*c/(2*pi*r) = 5000 [Hz].
% Define a vector of frequencies. These are the frequencies for which 
% the acoustic indicators R, T and A will be computed.
% ( !!! note that the frequency must be equal to the index of f !!!)
f         = 1:1:5000;    %Frequency range [Hz]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% MICROPHONE SPACINGS
% Microphone spacings; values for general transmission set up
s0 = 0.977; s1 = 0.033; s2 = 0.033;
l1 = 0.053; l2 = 0.141;
ze = 1.13+0.73; % measured by Indu, tube length
measures = [s0, s1, l1, l2, s2, ze];

% Distance between mics [0-1, 0-2, 0-3, 0-4]
x = [s0, s0+s1, s0+s1+l1+l2, s0+s1+l1+l2+s2];

%% CONSTANTS
% define parameters in air
constants.P     = 101.325;                      % Atmoshperic pressure [kPa]
constants.T     = temp;                         % Room temperature [K]
constants.c     = 20.047*sqrt(constants.T);     % Speed of sound in air [m/s]
% Density of air [kg/m^3]
constants.rho   = 1.290*(constants.P/101.325)*(273.15/constants.T); 

% define homogeneous layer parameters rho_L and c_L
constants.rho_L = 145; % [kg/m^3]
constants.c_L   = 250; %[ m/s]

%% DEFINE WAVENUMBER IN THE AIR AND IN THE LAYER
% air domain: characteristic impedance(create this variable to make shorter
% expression for the pressure coefficients) and wavenumber in the air
Z0 = constants.rho   * constants.c  ;
k  = 2*pi*f/constants.c;   

% define homogeneous layer characteristic impedance and wavenumber
ZL = constants.rho_L * constants.c_L;
a  = 2*pi*f/constants.c_L; %index corresponds to the frequency

% NOTE: index corresponds to the frequency, i.g.,
% a(100) is the wavenumber of homogeneous layer for f=100 Hz, or,
% k(500) is the wavenumber of air for f=500 Hz.

%% Delany-Bazley model for the layer
% sigma=10000; % flow resitivity
% X  = constants.rho*f/sigma;
% ZL = Z0*(1+0.057*X.^(-0.754)-1i*0.087.*X.^(-0.732)) ; 
% a  = k.*(1+0.0978*X.^(-0.700)-1i*0.189.*X.^(-0.595)); 
% %Note that the real part of the surface impedance is negative at low fre-
% % quencies, which is not physical. This defect in the Delany and Bazley
% % model has been noticed by Miki (1990) who proposed modifications to
% % Equations 4.102. Book Reference: "finite element..." [Atalla, 2015]


%% ANALYTIACL SOLUTION OF FORCED IMPEDANCE TUBE, FREQ. DOMAIN
%  NOTE1:  exp(iwt) convention used
%  NOTE2:  boundary condition downstream may be open or closed tube.

% impedance tube geometry
zs = 1.13; %distance from source(speaker) to sample surface(left)
zb = 0.73; %distance from sample surface(left) to the tube backing
% impedance tube radius (needed for open boundary condition)
r=(40e-3)/2 ; %[m]

%define excitation amplitude
p0 = 2e-5 ; % [Pa]


% analytical reflection and transmission coefficients for homogeneous layer
% check analytical_T_and_R.nb file for more information.
RA = ((-Z0^2 + ZL.^2).*sin(a.*d))./(-2*1i*Z0*ZL.*cos(a.*d) + ...
      (Z0^2 + ZL.^2).*sin(a.*d));
TA = (4*exp(1i*d*(a + k)).*Z0.*ZL)./(-(Z0 - ZL).^2 + ...
      exp(2*1i*a.*d).*(Z0 + ZL).^2);    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% choose if computes analytical solution of forced impedance tube %%%%%%
compute_synthetic_data_generation_tag = 1; %(=1) compute, (=0) don't
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% compute & extract coefficients from struct 'coeffs'
if compute_synthetic_data_generation_tag==1
    % analytical coefficients of forced impedance tube reponse
    coeffs=calc_coef_a(p0,d,zs,zb,r,Z0,ZL,k,a);
    %                  --|---------|---------
    %                  ^      ^         ^
    %                amp.  geometry   material      
    
    % hard termination v=0
    Ah=coeffs.Ah;
    Bh=coeffs.Bh;
    Ch=coeffs.Ch;
    Dh=coeffs.Dh;
    Eh=coeffs.Eh;
    Fh=coeffs.Fh;
    % open termination p=0
    Ao=coeffs.Ao;
    Bo=coeffs.Bo;
    Co=coeffs.Co;
    Do=coeffs.Do;
    Eo=coeffs.Eo;
    Fo=coeffs.Fo;

    %%%%%%% visualization purposes until the end of this codeblock%%%%%%%
    % plot the pressure field for a certain instant of time (say t=0) and
    % frequency of index defined below. 
    
    % choose frequency
    ind= 501; %index of frequency vector is the frequency itself
    
    % mesh-----------
    x1 =-zs:0.001: 0 ;
    x2 = 0 :0.001: d ;
    x3 = d :0.001: zb;
    %-----------------
    %1,2,3 here refers to regions 1 (upstream),2(sample) or 3(downstream)
    
    % complex pressure amplitudes for frequency = ind
    p1h = Ah(ind)*exp(-1i*k(ind).*x1) + Bh(ind)*exp(1i*k(ind).*x1);
    p2h = Eh(ind)*exp(-1i*a(ind).*x2) + Fh(ind)*exp(1i*a(ind).*x2);
    p3h = Ch(ind)*exp(-1i*k(ind).*x3) + Dh(ind)*exp(1i*k(ind).*x3);
    p1o = Ao(ind)*exp(-1i*k(ind).*x1) + Bo(ind)*exp(1i*k(ind).*x1);
    p2o = Eo(ind)*exp(-1i*a(ind).*x2) + Fo(ind)*exp(1i*a(ind).*x2);
    p3o = Co(ind)*exp(-1i*k(ind).*x3) + Do(ind)*exp(1i*k(ind).*x3);
    
    
%     % plot movie hard termination tube
%     frame = 0:1/f(ind)/80:1/f(ind);
%     figure(700)
%     phmax = max([real(p1h) real(p2h) real(p3h)]);
%     for kk=1:length(frame);
%         plot(x1,real(p1h * exp(1i*2*pi*f(ind)*frame(kk))),...
%              x2,real(p2h * exp(1i*2*pi*f(ind)*frame(kk))), ...
%              x3,real(p3h * exp(1i*2*pi*f(ind)*frame(kk)))      );
%         ylim([-phmax phmax]);
%         pause(.05)
%         xlabel('x [m]');
%         ylabel('pressure [Pa]');
%         title('hard(closed)');
%     end     
    
%     % plot movie open termination tube
%     frame = 0:1/f(ind)/80:1/f(ind);
%     figure(700)
%     phmax = max([real(p1o) real(p2o) real(p3o)]);
%     for kk=1:length(frame);
%         plot(x1,real(p1o * exp(1i*2*pi*f(ind)*frame(kk))),...
%              x2,real(p2o * exp(1i*2*pi*f(ind)*frame(kk))), ...
%              x3,real(p3o * exp(1i*2*pi*f(ind)*frame(kk)))      );
%         ylim([-phmax phmax]);
%         pause(.05)
%         xlabel('x [m]');
%         ylabel('pressure [Pa]');
%         title('hard(closed)');
%     end     
    

    %plot impedance tube pressure for t=0 (steady state)
    figure(500+ind) % 500 here is just a large number to saparate the plots
                      % between code blocks
    subplot(2,1,1)
    hold on;
    plot(x1,real(p1h) , "Color","b");
    plot(x2,real(p2h), "Color","b");
    plot(x3,real(p3h), "Color","b");
    xlabel('[m]');
    ylabel('pressure [Pa]');
    title('Closed tube');
    
%     plot(ch(:,1),ch(:,2))
%     legend('analytical','','','comsol');

    hold off;
    grid on;
    box on;
    subplot(2,1,2)
    hold on;
    plot(x1,real(p1o), "Color","b");
    plot(x2,real(p2o), "Color","b");
    plot(x3,real(p3o), "Color","b");
    xlabel('[m]');
    ylabel('pressure [Pa]');
    title('Open tube')
    
%     plot(co(:,1),co(:,2))
%     legend('analytical','','','comsol');

    hold off;
    grid on;
    box on;
end


%% TIME VECTOR WITH SAMPLING FREQUENCY fs
fs   = 48000; % [Hz]
time = (0:fs-1)*1/fs; %timestep of 1/fs

%% GENERATE SYNTHETIC TIME SIGNAL FROM ANALYTICAL SOLUTION
% 'microphone recording' for a period of 1 second

if compute_synthetic_data_generation_tag==1
    
    %define coordinate of each microphone, left sample surface x=0
    synmic0.x =-s0-s1-l1;
    synmic1.x =-s1-l1;
    synmic2.x =-l1;
    synmic3.x = l2;
    synmic4.x = l2+s2;
%     x = [s0, s0+s1, s0+s1+l1+l2, s0+s1+l1+l2+s2];

    
    %initialize to zero pressure mics
    synmic0.ph = zeros(size(time));
    synmic0.po = zeros(size(time));
    synmic1.ph = zeros(size(time));
    synmic1.po = zeros(size(time));
    synmic2.ph = zeros(size(time));
    synmic2.po = zeros(size(time));
    synmic3.ph = zeros(size(time));
    synmic3.po = zeros(size(time));
    synmic4.ph = zeros(size(time));
    synmic4.po = zeros(size(time));
    
   
    % assumes compute_synthetic_data_generation_tag==1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%choose to impose same amplitude response for every frequency(=1); or
    % to impose the same excitation amplitude p0 cte, but the response has
    %  different amplitudes (=0). %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    norm_tag =0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % define frequency content and wavenumbe of the excitation at x=-zs
    f_ex  = 100:1:2000;
%     f_ex  = [100 101 102 103 104 105 500 501 502 1000 1001 1002];
    k_ex  = 2*pi*f_ex/constants.c; 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    if length(f_ex)>length(f)
        error("Excitation has a frequency content outside of the compu"+...
              "ted analytical solution");
    end

    % NOTE FOR THE FOLLOWING CODEBLOCK WORK THE COEFFICIENT INDEX, SAY Ah,
    % NEEDS TO COINCIDE WITH THE FREQUENCY ITSELF. THIS IS POSSIBLE THANKS
    % TO THE DEFINITION ABOVE OF VECTOR f = 1:1:5000 WHERE THE FREQ. AND 
    % THE INDICES ARE THE SAME.

    switch norm_tag
       case 1 % same response amplitude for every frequency
        for ind=1:length(f_ex)

        %microphone 0 ( x<0 ) use solution A and B coeffs
        %hard
        aux       = (Ah(f_ex(ind))*exp(-1i     *k(ind).*synmic0.x) + ...
                     Bh(f_ex(ind))*exp( 1i     *k(ind).*synmic0.x)).*...
                              exp( 1i*2*pi*f_ex(ind).*time     )  ;
        synmic0.ph= synmic0.ph + real(aux)./(abs(aux))*p0;
        %open
        aux       = (Ao(f_ex(ind))*exp(-1i     *k(ind).*synmic0.x) + ...
                     Bo(f_ex(ind))*exp( 1i     *k(ind).*synmic0.x)).*...
                              exp( 1i*2*pi*f_ex(ind).*time     )  ;          
        synmic0.po= synmic0.po + real(aux)./(abs(aux))*p0;

                   
        
        %microphone 1 ( x<0 ) use solution A and B coeffs
        %hard
        aux       = (Ah(f_ex(ind))*exp(-1i     *k(ind).*synmic1.x) + ...
                     Bh(f_ex(ind))*exp( 1i      *k(ind).*synmic1.x)).*...
                              exp( 1i*2*pi*f_ex(ind).*time     )  ;
        synmic1.ph= synmic1.ph + real(aux)./(abs(aux))*p0;
        %open
        aux       = (Ao(f_ex(ind))*exp(-1i     *k(ind).*synmic1.x) + ...
                     Bo(f_ex(ind))*exp( 1i     *k(ind).*synmic1.x)).*...
                              exp( 1i*2*pi*f_ex(ind).*time     )  ;          
        synmic1.po= synmic1.po + real(aux)./(abs(aux))*p0;

                   
        
        %microphone 2 ( x<0 ) use solution A and B coeffs
        %hard
        aux       = (Ah(f_ex(ind))*exp(-1i     *k(ind).*synmic2.x) + ...
                     Bh(f_ex(ind))*exp( 1i     *k(ind).*synmic2.x)).*...
                              exp( 1i*2*pi*f_ex(ind).*time     )  ;
        synmic2.ph= synmic2.ph + real(aux)./(abs(aux))*p0;
        %open
        aux       = (Ao(f_ex(ind))*exp(-1i     *k(ind).*synmic2.x) + ...
                     Bo(f_ex(ind))*exp( 1i     *k(ind).*synmic2.x)).*...
                              exp( 1i*2*pi*f_ex(ind).*time     )  ;           
        synmic2.po= synmic2.po + real(aux)./(abs(aux))*p0;
              
%                 microphone 3 ( x>d ) use solution C and D coeffs
        %hard
        aux       = (Ch(f_ex(ind))*exp(-1i     *k(ind).*synmic3.x) + ...
                     Dh(f_ex(ind))*exp( 1i     *k(ind).*synmic3.x)).*...
                              exp( 1i*2*pi*f_ex(ind).*time     )  ;
        synmic3.ph= synmic3.ph + real(aux)./(abs(aux))*p0;
        %open     
        aux       = (Co(f_ex(ind))*exp(-1i     *k(ind).*synmic3.x) + ...
                     Do(f_ex(ind))*exp( 1i     *k(ind).*synmic3.x)).*...
                              exp( 1i*2*pi*f_ex(ind).*time     )  ;          
        synmic3.po= synmic3.po + real(aux)./(abs(aux))*p0;

                   
        
        %microphone 4 ( x>d ) use solution C and D coeffs
        %hard
        aux       = (Ch(f_ex(ind))*exp(-1i     *k(ind).*synmic4.x) + ...
                     Dh(f_ex(ind))*exp( 1i     *k(ind).*synmic4.x)).*...
                              exp( 1i*2*pi*f_ex(ind).*time     )  ;
        synmic4.ph= synmic4.ph + real(aux)./(abs(aux))*p0;
        %open           
        aux       = (Co(f_ex(ind))*exp(-1i     *k(ind).*synmic4.x) + ...
                     Do(f_ex(ind))*exp( 1i     *k(ind).*synmic4.x)).*...
                              exp( 1i*2*pi*f_ex(ind).*time     )  ;
        synmic4.po= synmic4.po + real(aux)./(abs(aux))*p0;

                   
        end
       case 0 % same input(speaker) amplitude for every frequency
        for ind=1:length(f_ex)

                    %microphone 0 ( x<0 ) use solution A and B coeffs
                    aux       = (Ah(f_ex(ind))*exp(-1i     *k_ex(ind).*synmic0.x) + ...
                                 Bh(f_ex(ind))*exp( 1i     *k_ex(ind).*synmic0.x)).*...
                                          exp( 1i*2*pi*f_ex(ind).*time     )  ;
                    synmic0.ph= synmic0.ph + real(aux);
                    
                    aux       = (Ao(f_ex(ind))*exp(-1i     *k_ex(ind).*synmic0.x) + ...
                                 Bo(f_ex(ind))*exp( 1i     *k_ex(ind).*synmic0.x)).*...
                                          exp( 1i*2*pi*f_ex(ind).*time     )  ;          
                    synmic0.po= synmic0.po + real(aux);
                               
                    
                    %microphone 1 ( x<0 ) use solution A and B coeffs
                    aux       = (Ah(f_ex(ind))*exp(-1i      *k_ex(ind).*synmic1.x) + ...
                                 Bh(f_ex(ind))*exp( 1i      *k_ex(ind).*synmic1.x)).*...
                                          exp( 1i*2*pi*f_ex(ind).*time     )  ;
                    synmic1.ph= synmic1.ph + real(aux);
    
                    aux       = (Ao(f_ex(ind))*exp(-1i     *k_ex(ind).*synmic1.x) + ...
                                 Bo(f_ex(ind))*exp( 1i     *k_ex(ind).*synmic1.x)).*...
                                          exp( 1i*2*pi*f_ex(ind).*time     )  ;          
                    synmic1.po= synmic1.po + real(aux);
                               
                    
                    %microphone 2 ( x<0 ) use solution A and B coeffs
                    aux       = (Ah(f_ex(ind))*exp(-1i     *k_ex(ind).*synmic2.x) + ...
                                 Bh(f_ex(ind))*exp( 1i     *k_ex(ind).*synmic2.x)).*...
                                          exp( 1i*2*pi*f_ex(ind).*time     )  ;
                    synmic2.ph= synmic2.ph + real(aux);
    
                    aux       = (Ao(f_ex(ind))*exp(-1i     *k_ex(ind).*synmic2.x) + ...
                                 Bo(f_ex(ind))*exp( 1i     *k_ex(ind).*synmic2.x)).*...
                                          exp( 1i*2*pi*f_ex(ind).*time     )  ;           
                    synmic2.po= synmic2.po + real(aux);
                               
                    %microphone 3 ( x>d ) use solution C and D coeffs
                    aux       = (Ch(f_ex(ind))*exp(-1i     *k_ex(ind).*synmic3.x) + ...
                                 Dh(f_ex(ind))*exp( 1i     *k_ex(ind).*synmic3.x)).*...
                                          exp( 1i*2*pi*f_ex(ind).*time     )  ;
                    synmic3.ph= synmic3.ph + real(aux);
                               
                    aux       = (Co(f_ex(ind))*exp(-1i     *k_ex(ind).*synmic3.x) + ...
                                 Do(f_ex(ind))*exp( 1i     *k_ex(ind).*synmic3.x)).*...
                                          exp( 1i*2*pi*f_ex(ind).*time     )  ;          
                    synmic3.po= synmic3.po + real(aux);
                               
                    
                    %microphone 4 ( x>d ) use solution C and D coeffs
                    aux       = (Ch(f_ex(ind))*exp(-1i     *k_ex(ind).*synmic4.x) + ...
                                 Dh(f_ex(ind))*exp( 1i     *k_ex(ind).*synmic4.x)).*...
                                          exp( 1i*2*pi*f_ex(ind).*time     )  ;
                    synmic4.ph= synmic4.ph + real(aux);
                               
                    aux       = (Co(f_ex(ind))*exp(-1i     *k_ex(ind).*synmic4.x) + ...
                                 Do(f_ex(ind))*exp( 1i     *k_ex(ind).*synmic4.x)).*...
                                          exp( 1i*2*pi*f_ex(ind).*time     )  ;
                    synmic4.po= synmic4.po + real(aux);
                               
                    end
    end

end

%% COMPUTE TRANSFER FUNCTIONS H_0i
 % assign mic data to struct IRhxc and IRoxc
 IRhxc.mic0 = synmic0.ph;
 IRoxc.mic0 = synmic0.po;

 IRhxc.mic1 = synmic1.ph;
 IRoxc.mic1 = synmic1.po;

 IRhxc.mic2 = synmic2.ph;
 IRoxc.mic2 = synmic2.po;

 IRhxc.mic3 = synmic3.ph;
 IRoxc.mic3 = synmic3.po;

 IRhxc.mic4 = synmic4.ph;
 IRoxc.mic4 = synmic4.po;

 % FAST FOURRIER TRANSFORM AND COMPUTE TRANSFER FUNCTIONS 
 Hh = Fast_Fourier(IRhxc,fs);
 Ho = Fast_Fourier(IRoxc,fs);

%% plot abs H functions of microphones hard/open termination
    L=fs; %length of time signal
    fp=(0:L-1)/L*fs;
    figure(300);
    subplot(2,1,1);
    hold on;
    grid on;
    plot(fp,abs(Hh(:,1)));
    plot(fp,abs(Hh(:,2)));
    plot(fp,abs(Hh(:,3)));
    plot(fp,abs(Hh(:,4)));
    xlim([0 5000]);
    xlabel('frequency[Hz]');
    ylabel('hard : abs Hh')
    legend('0-1','0-2','0-3','0-4');
    title('hard : H');
    hold off;
    subplot(2,1,2);
    hold on;
    grid on;
    plot(fp,abs(Ho(:,1)));
    plot(fp,abs(Ho(:,2)));
    plot(fp,abs(Ho(:,3)));
    plot(fp,abs(Ho(:,4)));
    
    xlim([0 5000]);
    xlabel('frequency[Hz]');
    ylabel('open : abs Ho');
    legend('0-1','0-2','0-3','0-4');
    title('open : H');
    hold off;

%% COMPUTE AMPLITUDE COEFFICIENTS
% the following codeblock is copied from calc_coef_two_load.m
% TODO: use the function calc_coef_two_load.m instead of the following.

 % geometry
 m  = measures; 
 % declare vector with coefficients for each frequency
 TC = zeros(1,length(f)); 
 RC = zeros(1,length(f)); 
 
 countnoise=0; 
 for ff = 1:length(f) %loop on freqs you wish to compute the ac indicators
    
    % RECAP
    % Coordinates of each microphone (note: left sample surface is at x=0)
    %     synmic0.x =-s0-s1-l1;
    %     synmic1.x =-s1-l1;
    %     synmic2.x =-l1;
    %     synmic3.x = l2;
    %     synmic4.x = l2+s2;
    %     x = [s0, s0+s1, s0+s1+l1+l2, s0+s1+l1+l2+s2];
    %     measures = [s0, s1, l1, l2, s2, ze];
    %
    %     equivalent relations to compare
    %     x0 = -m(1)-m(2)-m(3);
    %     x1 = -m(2)-m(3);
    %     x2 = -m(3);
    %  x1-x2 = -m(2);
    %     x3 =  m(4);
    %     x4 =  m(4)+m(5);
    %  x3-x4 = -m(5);

    %Calculation A, B, C and D matrices hard backing
    Ahe = 1j * (Hh(f(ff)+1, 1)*exp( 1j*k(ff)*(-m(3))      )...
               -Hh(f(ff)+1, 2)*exp( 1j*k(ff)*(-m(3)-m(2)) )   )...
                /(    2*sin( k(ff) * (-m(2)) )   );
    
    Bhe = 1j * (Hh(f(ff)+1, 2)*exp(-1j*k(ff)*(-m(3)-m(2)) )...
               -Hh(f(ff)+1, 1)*exp(-1j*k(ff)*(-m(3))      )   )...
                /(    2*sin( k(ff) * (-m(2)) )   );
    
    Che = 1j * (Hh(f(ff)+1,3)*exp( 1j*k(ff)*(m(4)+m(5))  )...
               -Hh(f(ff)+1,4)*exp( 1j*k(ff)* m(4)        )   )...
               /(    2*sin( k(ff) * (-m(5))  )    ); 
    
    Dhe = 1j * (Hh(f(ff)+1,4)*exp(-1j*k(ff)* m(4)        )...
               -Hh(f(ff)+1,3)*exp(-1j*k(ff)* (m(4)+m(5)) )   )...
               /(    2*sin( k(ff) * (-m(5))  )    ); 
  
    
    %Calculation A, B, C and D matrices open backing
    Aoe = 1j * (Ho(f(ff)+1, 1)*exp( 1j*k(ff)*(-m(3))      )...
               -Ho(f(ff)+1, 2)*exp( 1j*k(ff)*(-m(3)-m(2)) )   )...
                /(    2*sin( k(ff) * (-m(2)) )   );
    
    Boe = 1j * (Ho(f(ff)+1, 2)*exp(-1j*k(ff)*(-m(3)-m(2)) )...
               -Ho(f(ff)+1, 1)*exp(-1j*k(ff)*(-m(3))      )   )...
                /(    2*sin( k(ff) * (-m(2)) )   );
    
    Coe = 1j * (Ho(f(ff)+1,3)*exp( 1j*k(ff)*(m(4)+m(5))  )...
               -Ho(f(ff)+1,4)*exp( 1j*k(ff)* m(4)        )   )...
               /(    2*sin( k(ff) * (-m(5))  )    ); 
    
    Doe = 1j * (Ho(f(ff)+1,4)*exp(-1j*k(ff)* m(4)        )...
               -Ho(f(ff)+1,3)*exp(-1j*k(ff)* (m(4)+m(5)) )   )...
               /(    2*sin( k(ff) * (-m(5))  )    ); 
    

    % Compute coefficients R and T
    if(abs(Ahe*Doe-Aoe*Dhe)>1e-10)% if denominator isnt noise
    TC(ff) = (Che*Doe-Coe*Dhe)/(Ahe*Doe-Aoe*Dhe); % Transmission Coefficient [-]
    RC(ff) = (Bhe*Doe-Boe*Dhe)/(Ahe*Doe-Aoe*Dhe); % Reflection coefficient
    else
     countnoise=countnoise+1; 
     %fprintf("calc_coef: noise detected for frequency %f\n", ff);
     TC(ff) = (Che*Doe-Coe*Dhe); % Transmission Coefficient [-]
     RC(ff) = (Bhe*Doe-Boe*Dhe);
    end
    
    % compute absorption coefficient
    AC(ff) = 1-(abs(RC(ff))).^2-(abs(TC(ff))).^2; %Absorbtion coefficient
    
    % Calculation Transmission Loss
    TL(ff) = 20*log10(abs(1/TC(ff)));   % Transmission loss [dB]
 end

 fprintf(['calc_coef_two_load called\n      number of noise detected: '...
          '%d out of %d NFFT data points\n\n'], countnoise, length(f));



%% COMPUTE T AND R FROM AMPLITUDE COEFFICIENTS
% TC, RC matches with TA and RA

 
 if (length(f)>1)
    figure(1) 
    hold on;
    grid on;
    plot(f, abs(RC), 'r--', "LineWidth", 1.5);
    plot(f, abs(RA), 'r'  );
    plot(f, abs(TC), 'b--', "LineWidth", 1.5);
    plot(f, abs(TA), 'b'  );
    legend('RC Coefficients', 'RA coefficient [-]',...
           'TC Coefficient',  'TA coefficient');
           
    hold off;
 else
     TA
     TC
     RA
     RC
 end

 
