%% README
% Post processing script of experimental impedance tube. Adapted from
% MSc L. Schijff, Ashoka karunarathne, N.A. van de Straat, 
% by R. Liupekevicius 20-12-2021. Eindhoven University of Technology
%
% OBJECTIVES
% this script is a tentative to verify the two load postprocessing
% procedure by pluggin
% In case of questions, contact r.liupekevicius.carnielli@tue.nl
%
% INSTRUCTIONS
% Main settings are commented by a sandwitch of '%%%%' with a short
% description on it. Please see the example below.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% EXAMPLE of settings to change %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parameter_to_change = 1; % change the variable 'parameter_to_chage'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%
%
%
%% Impedence tube transmission analysis two load method
% ------------------------------------------------------------------------
clear all   % clear all variables in the workspace
close all   % close figures. Comment if multiple samples need to be plotted
clc         % clears the command window

%% ADD FUNCTIONS TO PATH
% % add the shared functions to the search path
% path(['C:\GitLab\Tools\impedance_tube_post_processing\shared_functions'],...
%      path);       

%% PARAMETERS
% Check and change these parameters for every run!!!
% ------------------------------------------------------------------------
% Data directory specification
% 'Data_path' + 'Folder' + 'Sample' should form the complete directory to the
% file with the mic data

%Path to general data folder
Data_path = erase(pwd,'\test') 

%Specific file of sample type or date
Folder = 'data'                                           


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% choose experimental data to be post processed %%%%%%%%%%%%%%%%%
Sample    = 'sample_wag_PU_ref';   %Sample folder that contains audio files    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Complete path to desired data set
Sample_path = fullfile(Data_path, Folder,Sample)    


%% PLOT SETTINGS
% Plot settings; line color, smoot/raw data, axis ranges
% Change to desired values!!!

                 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Plotting of smooth (=1) or raw(=0) data %%%%%%%%%%%%%%%%%%%%%%
smooth      = 1;                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% next 6 lines are other settings from Lotte's script
line        = '--b';                %Line color coding --> See 'plot' in Matlab documentation
width       = 1.5;                  %Line width
axis_TL     = [100 5000 -1 15];     %Axis range transmission loss - [x-axis lower limit, x-axis upper limit, y-axis lower limit, y-axis upper limit]
axis_coef   = [100 5000 0 1];    	%Axis range coefficients - [x-axis lower limit, x-axis upper limit, y-axis lower limit, y-axis upper limit]

save_results = 1;       %Indicater to save data to .txt file yes=1 and no=0
                        %The saved data is structured as:
                        %column 1: Frequency [Hz]
                        %column 2: Transmission loss [dB]
                        %column 3: Transmission coefficient [-]
                        %column 4: Reflection coefficient [-]
                        %column 5: Absorption coefficient [-]


%% MEASUREMENT SPECS
% Measurement specifications
d         = 0.01;       % sample thickness [m] Not important for two-load method ????? I don't agree in general
temp      = 293.15;     % temperature during measurement [K]

% frequency range, planewave up to 1.84*c/(2*pi*r)
f         = 1:1:5e3;    %Frequency range [Hz]


% impedance tube radius
r=(40e-3)/2 ; %[m]

                    
%% MICROPHONE SPACINGS

% Microphone spacings; values for general transmission set up
s0 = 0.977; s1 = 0.033; s2 = 0.033;
l1 = 0.053; l2 = 0.141;
% ze = 3.0; Lotte's. It looks like this variable is not used
ze = 1.13+0.73; % measured by Indu

measures = [s0, s1, l1, l2, s2, ze];

% Distance between mics [0-1, 0-2, 0-3, 0-4]
x = [s0, s0+s1, s0+s1+l1+l2, s0+s1+l1+l2+s2];


%% CONSTANTS
constants.P     = 101.325;                      % Atmoshperic pressure [kPa]
constants.T     = temp;                         % Room temperature [K]
constants.c     = 20.047*sqrt(constants.T);     % Speed of sound in air [m/s]
constants.rho   = 1.290*(constants.P/101.325)*(273.15/constants.T); % Density of air [kg/m^3]

%% Analytical solution forced impedance tube (closed and open)
% from analytical solution of homogeneous layer ( exp(iwt) convention)

% impedance tube geometry
zs = 1.13; %distance from source(speaker) to sample surface(left)
zb = 0.73; %distance from sample surface(left) to the tube backing
%define excitation amplitude
p0 = 2e-5 ; % [Pa]

% air domain characteristic impedance (define this to make shorter
% expression for the coefficients)
Z0 = constants.rho   * constants.c  ;
k = 2*pi*f/constants.c; % define wave vector in the air


% define homogeneous layer parameters rho_L and c_L
constants.rho_L = 145; % [kg/m^3]
constants.c_L   = 250; %[ m/s]
ZL              = constants.rho_L * constants.c_L;
a               = 2*pi*f/constants.c_L;

% % Delany-Bazley model
% sigma=10000; % flow resitivity
% X  = constants.rho*f/sigma;
% ZL = Z0*(1+0.057*X.^(-0.754)-1i*0.087.*X.^(-0.732)) ; 
% a  = k.*(1+0.0978*X.^(-0.700)-1i*0.189.*X.^(-0.595)); 
% %Note that the real part of the surface impedance is negative at low fre-
% % quencies, which is not physical. This defect in the Delany and Bazley
% % model has been noticed by Miki (1990) who proposed modifications to
% % Equations 4.102.

% analytical reflection and transmission coefficients for homogeneous layer
% check analytical_T_and_R.nb file for more information.
RA = ((-Z0^2 + ZL.^2).*sin(a.*d))./(-2*1i*Z0*ZL.*cos(a.*d) + ...
      (Z0^2 + ZL.^2).*sin(a.*d));
TA = (4*exp(1i*d*(a + k)).*Z0.*ZL)./(-(Z0 - ZL).^2 + ...
      exp(2*1i*a.*d).*(Z0 + ZL).^2);    




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% choose if computes analytical solution of forced impedance tube %%%%%%
compute_synthetic_data_generation_tag = 0; %(=1) compute, (=0) don't
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% extract coefficients from struct
if compute_synthetic_data_generation_tag==1
    % analytical coefficients of forced impedance tube reponse
    coeffs=calc_coef_a(p0,d,zs,zb,r,Z0,ZL,k,a);
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
    % frequency of index defined below. p1h to p3o variables are meant 
    % only for the visualization purposes.
    
    ind= 500;
    x1 =-zs:0.01: 0 ;
    x2 = 0 :0.01: d ;
    x3 = d :0.01: zb;
    %1,2,3 here refers to regions 1 (upstream),2(sample) or 3(downstream)
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

%% TIME VECTOR, SAMPLING FREQUENCY fs
fs   = 48000;
time = (0:fs-1)*1/fs; %timestep of 1/fs

%% Generate synthetic 1-second time signal from analytical solution
% microphone recording for a period of 1 second
% with sampling frequency fs

if compute_synthetic_data_generation_tag==1
    
    %define coordinate of each microphone, left sample surface x=0
    synmic0.x =-s0-s1-l1;
    synmic1.x =-s1-l1;
    synmic2.x =-l1;
    synmic3.x = l2;
    synmic4.x = l2+s2;
    
    
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
    %%%%%%% choose to impose same amplitude for every frequency (=1) or to
    % impose the same excitation amplitude p0 cte, but the response has
    %  different amplitudes (=0). %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    norm_tag =0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    fmin  = 10 ;
    fstep = 10   ;
    fmax  = 2000;
    switch norm_tag
        case 1
                    for ind= fmin:fstep:fmax 
    
                    %microphone 0 ( x<0 ) use solution A and B coeffs
                    %hard
                    aux       = (Ah(ind)*exp(-1i     *k(ind).*synmic0.x) + ...
                                 Bh(ind)*exp( 1i     *k(ind).*synmic0.x)).*...
                                          exp( 1i*2*pi*f(ind).*time     )  ;
                    synmic0.ph= synmic0.ph + real(aux)./(abs(aux))*p0;
                    %open
                    aux       = (Ao(ind)*exp(-1i     *k(ind).*synmic0.x) + ...
                                 Bo(ind)*exp( 1i     *k(ind).*synmic0.x)).*...
                                          exp( 1i*2*pi*f(ind).*time     )  ;          
                    synmic0.po= synmic0.po + real(aux)./(abs(aux))*p0;
    
                               
                    
                    %microphone 1 ( x<0 ) use solution A and B coeffs
                    %hard
                    aux       = (Ah(ind)*exp(-1i     *k(ind).*synmic1.x) + ...
                                 Bh(ind)*exp( 1i      *k(ind).*synmic1.x)).*...
                                          exp( 1i*2*pi*f(ind).*time     )  ;
                    synmic1.ph= synmic1.ph + real(aux)./(abs(aux))*p0;
                    %open
                    aux       = (Ao(ind)*exp(-1i     *k(ind).*synmic1.x) + ...
                                 Bo(ind)*exp( 1i     *k(ind).*synmic1.x)).*...
                                          exp( 1i*2*pi*f(ind).*time     )  ;          
                    synmic1.po= synmic1.po + real(aux)./(abs(aux))*p0;
    
                               
                    
                    %microphone 2 ( x<0 ) use solution A and B coeffs
                    %hard
                    aux       = (Ah(ind)*exp(-1i     *k(ind).*synmic2.x) + ...
                                 Bh(ind)*exp( 1i     *k(ind).*synmic2.x)).*...
                                          exp( 1i*2*pi*f(ind).*time     )  ;
                    synmic2.ph= synmic2.ph + real(aux)./(abs(aux))*p0;
                    %open
                    aux       = (Ao(ind)*exp(-1i     *k(ind).*synmic2.x) + ...
                                 Bo(ind)*exp( 1i     *k(ind).*synmic2.x)).*...
                                          exp( 1i*2*pi*f(ind).*time     )  ;           
                    synmic2.po= synmic2.po + real(aux)./(abs(aux))*p0;
                          
    %                 microphone 3 ( x>d ) use solution C and D coeffs
                    %hard
                    aux       = (Ch(ind)*exp(-1i     *k(ind).*synmic3.x) + ...
                                 Dh(ind)*exp( 1i     *k(ind).*synmic3.x)).*...
                                          exp( 1i*2*pi*f(ind).*time     )  ;
                    synmic3.ph= synmic3.ph + real(aux)./(abs(aux))*p0;
                    %open     
                    aux       = (Co(ind)*exp(-1i     *k(ind).*synmic3.x) + ...
                                 Do(ind)*exp( 1i     *k(ind).*synmic3.x)).*...
                                          exp( 1i*2*pi*f(ind).*time     )  ;          
                    synmic3.po= synmic3.po + real(aux)./(abs(aux))*p0;
    
                               
                    
                    %microphone 4 ( x>d ) use solution C and D coeffs
                    %hard
                    aux       = (Ch(ind)*exp(-1i     *k(ind).*synmic4.x) + ...
                                 Dh(ind)*exp( 1i     *k(ind).*synmic4.x)).*...
                                          exp( 1i*2*pi*f(ind).*time     )  ;
                    synmic4.ph= synmic4.ph + real(aux)./(abs(aux))*p0;
                    %open           
                    aux       = (Co(ind)*exp(-1i     *k(ind).*synmic4.x) + ...
                                 Do(ind)*exp( 1i     *k(ind).*synmic4.x)).*...
                                          exp( 1i*2*pi*f(ind).*time     )  ;
                    synmic4.po= synmic4.po + real(aux)./(abs(aux))*p0;
    
                               
                    end
        case 0
                    for ind= fmin:fstep:fmax
                    %microphone 0 ( x<0 ) use solution A and B coeffs
                    aux       = (Ah(ind)*exp(-1i     *k(ind).*synmic0.x) + ...
                                 Bh(ind)*exp( 1i     *k(ind).*synmic0.x)).*...
                                          exp( 1i*2*pi*f(ind).*time     )  ;
                    synmic0.ph= synmic0.ph + real(aux);
                    
                    aux       = (Ao(ind)*exp(-1i     *k(ind).*synmic0.x) + ...
                                 Bo(ind)*exp( 1i     *k(ind).*synmic0.x)).*...
                                          exp( 1i*2*pi*f(ind).*time     )  ;          
                    synmic0.po= synmic0.po + real(aux);
                               
                    
                    %microphone 1 ( x<0 ) use solution A and B coeffs
                    aux       = (Ah(ind)*exp(-1i      *k(ind).*synmic1.x) + ...
                                 Bh(ind)*exp( 1i      *k(ind).*synmic1.x)).*...
                                          exp( 1i*2*pi*f(ind).*time     )  ;
                    synmic1.ph= synmic1.ph + real(aux);
    
                    aux       = (Ao(ind)*exp(-1i     *k(ind).*synmic1.x) + ...
                                 Bo(ind)*exp( 1i     *k(ind).*synmic1.x)).*...
                                          exp( 1i*2*pi*f(ind).*time     )  ;          
                    synmic1.po= synmic1.po + real(aux);
                               
                    
                    %microphone 2 ( x<0 ) use solution A and B coeffs
                    aux       = (Ah(ind)*exp(-1i     *k(ind).*synmic2.x) + ...
                                 Bh(ind)*exp( 1i     *k(ind).*synmic2.x)).*...
                                          exp( 1i*2*pi*f(ind).*time     )  ;
                    synmic2.ph= synmic2.ph + real(aux);
    
                    aux       = (Ao(ind)*exp(-1i     *k(ind).*synmic2.x) + ...
                                 Bo(ind)*exp( 1i     *k(ind).*synmic2.x)).*...
                                          exp( 1i*2*pi*f(ind).*time     )  ;           
                    synmic2.po= synmic2.po + real(aux);
                               
                    %microphone 3 ( x>d ) use solution C and D coeffs
                    aux       = (Ch(ind)*exp(-1i     *k(ind).*synmic3.x) + ...
                                 Dh(ind)*exp( 1i     *k(ind).*synmic3.x)).*...
                                          exp( 1i*2*pi*f(ind).*time     )  ;
                    synmic3.ph= synmic3.ph + real(aux);
                               
                    aux       = (Co(ind)*exp(-1i     *k(ind).*synmic3.x) + ...
                                 Do(ind)*exp( 1i     *k(ind).*synmic3.x)).*...
                                          exp( 1i*2*pi*f(ind).*time     )  ;          
                    synmic3.po= synmic3.po + real(aux);
                               
                    
                    %microphone 4 ( x>d ) use solution C and D coeffs
                    aux       = (Ch(ind)*exp(-1i     *k(ind).*synmic4.x) + ...
                                 Dh(ind)*exp( 1i     *k(ind).*synmic4.x)).*...
                                          exp( 1i*2*pi*f(ind).*time     )  ;
                    synmic4.ph= synmic4.ph + real(aux);
                               
                    aux       = (Co(ind)*exp(-1i     *k(ind).*synmic4.x) + ...
                                 Do(ind)*exp( 1i     *k(ind).*synmic4.x)).*...
                                          exp( 1i*2*pi*f(ind).*time     )  ;
                    synmic4.po= synmic4.po + real(aux);
                               
                    end
    end

end



 %% LOAD SYNTHETIC DATA (COUPLE TO EXPERIMENTAL POST PROCESSING)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% skip (=1), don't skip (=0) %%%%%%%%%%%%%%%%%%%%%%%%%
skip_cross_correlation = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if compute_synthetic_data_generation_tag==1
    
     switch skip_cross_correlation
        case 1
             %use the following for skiping cross correlation procedure.
             % assign directly to the final struct before fft
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
        case 0
             %  use the following when correlation procedure is intended.
             % assign to variables to be cross corralated
             IRh(:,1) =  synmic0.ph;
             IRo(:,1) =  synmic0.po;
            
             IRh(:,2) =  synmic1.ph;
             IRo(:,2) =  synmic1.po;
            
             IRh(:,3) =  synmic2.ph;
             IRo(:,3) =  synmic2.po;
            
             IRh(:,4) =  synmic3.ph;
             IRo(:,4) =  synmic3.po;
            
             IRh(:,5) =  synmic4.ph;
             IRo(:,5) =  synmic4.po;

            
             
    end
end

%% EXPERIMENTAL POST PROCESSING STARTS FROM THIS LINE ON

%% LOADING + CUT MIC DATA

%if compute_syn... tag is off, then read the measure signal is read
if compute_synthetic_data_generation_tag==0
    % Load full data files hard & open backing
    [fs, IRmich] = Load_mic_data(Sample_path, 'Hard');
    [~, IRmico]  = Load_mic_data(Sample_path, 'Open');
    
    %Cut data to first second
    for ii = 1:length(IRmich(1,:))
        IRh(:,ii) = IRmich(1:fs, ii);
        IRo(:,ii) = IRmico(1:fs, ii);
    end

end

%% PLOT TIME SIGNAL OF MICS, termination hard and open
if compute_synthetic_data_generation_tag==1
    figure(600)
    subplot(2,1,1)
    hold on;
    plot(1000*time,real(synmic0.ph));
    plot(1000*time,real(synmic1.ph));
    plot(1000*time,real(synmic2.ph));
    plot(1000*time,real(synmic3.ph));
    plot(1000*time,real(synmic4.ph));
    % xlim([0 1000*5*1/f(ind)]); %5 periods
    xlabel("time(miliseconds)");
    ylabel("pressure[Pa]");
    title("hard tube ending (closed)");
    hold off;
    grid on;
    legend('0','1','2','3','4');
    subplot(2,1,2)
    hold on;
    plot(1000*time,real(synmic0.po));
    plot(1000*time,real(synmic1.po));
    plot(1000*time,real(synmic2.po));
    plot(1000*time,real(synmic3.po));
    plot(1000*time,real(synmic4.po));
    % xlim([0 1000*5*1/f(ind)]); %5 periods
    xlabel("time(miliseconds)");
    ylabel("pressure[Pa]");
    title("open tube ending");
    hold off;
    grid on;
    legend('0','1','2','3','4');
else %if compute_synthetic_data_generation_tag==0
    figure(600)
    subplot(2,1,1)
    hold on;
    plot(1000*time,IRh(:,1));
    plot(1000*time,IRh(:,2));
    plot(1000*time,IRh(:,3));
    plot(1000*time,IRh(:,4));
    plot(1000*time,IRh(:,5));
    % xlim([0 1000*5*1/f(ind)]); %5 periods
    xlabel("time(miliseconds)");
    ylabel("pressure[Pa]");
    title("hard tube ending (closed)");
    hold off;
    grid on;
    legend('0','1','2','3','4');
    subplot(2,1,2)
    hold on;
    plot(1000*time,IRo(:,1));
    plot(1000*time,IRo(:,2));
    plot(1000*time,IRo(:,3));
    plot(1000*time,IRo(:,4));
    plot(1000*time,IRo(:,5));
    % xlim([0 1000*5*1/f(ind)]); %5 periods
    xlabel("time(miliseconds)");
    ylabel("pressure[Pa]");
    title("open tube ending");
    hold off;
    grid on;
    legend('0','1','2','3','4');
end

%% Plot FFT of raw time signals

if compute_synthetic_data_generation_tag==0
    % plot abs fft of microphones
    L=fs; %length of time signal
    fp=(0:L-1)/L*fs-fs/2;
    figure(300);
    subplot(2,1,1);
    hold on;
    plot(fp,abs(fftshift(2/L*fft(IRh(:,1)))));
    plot(fp,abs(fftshift(2/L*fft(IRh(:,2)))));
    plot(fp,abs(fftshift(2/L*fft(IRh(:,3)))));
    plot(fp,abs(fftshift(2/L*fft(IRh(:,4)))));
    plot(fp,abs(fftshift(2/L*fft(IRh(:,5)))));
    xlim([0 5000]);
    xlabel('frequency[Hz]');
    ylabel('hard : abs 2 fft/L')
    legend('0','1','2','3','4');
    title('hard : fft before mic correlation');
    hold off;
    subplot(2,1,2);
    hold on;
    plot(fp,abs(fftshift(2/L*fft(IRo(:,1)))));
    plot(fp,abs(fftshift(2/L*fft(IRo(:,2)))));
    plot(fp,abs(fftshift(2/L*fft(IRo(:,3)))));
    plot(fp,abs(fftshift(2/L*fft(IRo(:,4)))));
    plot(fp,abs(fftshift(2/L*fft(IRo(:,5)))));
    xlim([0 5000]);
    xlabel('frequency[Hz]');
    ylabel('open : abs 2 fft/L');
    legend('0','1','2','3','4');
    title('open : fft before mic correlation');
    hold off;
else
   % plot abs fft of microphones
    L=fs; %length of time signal
    fp=(0:L-1)/L*fs-fs/2;
    figure(300);
    subplot(2,1,1);
    hold on;
    plot(fp,abs(fftshift(2/L*fft(synmic0.ph))));
    plot(fp,abs(fftshift(2/L*fft(synmic1.ph))));
    plot(fp,abs(fftshift(2/L*fft(synmic2.ph))));
    plot(fp,abs(fftshift(2/L*fft(synmic3.ph))));
    plot(fp,abs(fftshift(2/L*fft(synmic4.ph))));
    xlim([0 5000]);
    xlabel('frequency[Hz]');
    ylabel('hard : abs 2 fft/L')
    legend('0','1','2','3','4');
    title('hard : fft before mic correlation');
    hold off;
    subplot(2,1,2);
    hold on;
    plot(fp,abs(fftshift(2/L*fft(synmic0.po))));
    plot(fp,abs(fftshift(2/L*fft(synmic1.po))));
    plot(fp,abs(fftshift(2/L*fft(synmic2.po))));
    plot(fp,abs(fftshift(2/L*fft(synmic3.po))));
    plot(fp,abs(fftshift(2/L*fft(synmic4.po))));
    xlim([0 5000]);
    xlabel('frequency[Hz]');
    ylabel('open : abs 2 fft/L');
    legend('0','1','2','3','4');
    title('open : fft before mic correlation');
    hold off;
end

%% plot fft of ith microphone REAL and IMAG parts (synthetic signal)
% L=fs; %length of time signal
% fp=(0:L-1)/L*fs-fs/2;
% mici = synmic1.ph;  % pressure time signal of mici, BC hard or open
% figure(300);
% subplot(2,1,1);
% plot(fp(1:length(mici)),real(fftshift(2/L*fft(mici))));
% xlim([-2000 2000]);
% xlabel('frequency[Hz]');
% ylabel('real 2 fft/L')
% subplot(2,1,2);
% plot(fp(1:length(mici)),imag(fftshift(2/L*fft(mici))));
% xlim([-2000 2000]);
% xlabel('frequency[Hz]');
% ylabel('imag 2 fft/L');


%% plot fft of ith microphone ABS and PHASE parts (synthetic signal)
% fp=(0:L-1)/L*fs-fs/2;
% mici = 1;  % e.g. Hi0 is index i
% figure(400+mici);
% subplot(2,1,1);
% plot(fp(1:length(Hh(:,mici))),abs(fftshift(2/L*fft(Hh(:,mici)))));
% xlim([-2000 2000]);
% xlabel('frequency[Hz]');
% ylabel('abs 2 fft/L')
% subplot(2,1,2);
% plot(fp(1:length(Hh(:,mici))),unwrap(180/pi*angle(fftshift(2/L*fft(Hh(:,mici))))));
% xlim([-2000 2000]);
% xlabel('frequency[Hz]');
% ylabel('phase 2 fft/L');




%% CROSS CORRELATION
%Amount of samples for the sound to travel the length of the tube
N = round(ze/constants.c*fs); 

if skip_cross_correlation == 0
    %Cross correlation
    IRhxc = Cross_correlation(IRh, N, x, constants, fs);
    IRoxc = Cross_correlation(IRo, N, x, constants, fs);
    
    %plot fft abs
    fp=(0:L-1)/L*fs-fs/2;
    figure(400);
    subplot(2,1,1)
    hold on;
    plot(fp(1:length(IRhxc.mic0)),abs(fftshift(2/L*fft(IRhxc.mic0))));
    plot(fp(1:length(IRhxc.mic1)),abs(fftshift(2/L*fft(IRhxc.mic1))));
    plot(fp(1:length(IRhxc.mic2)),abs(fftshift(2/L*fft(IRhxc.mic2))));
    plot(fp(1:length(IRhxc.mic3)),abs(fftshift(2/L*fft(IRhxc.mic3))));
    plot(fp(1:length(IRhxc.mic4)),abs(fftshift(2/L*fft(IRhxc.mic4))));
    xlim([0 5000]);
    xlabel('frequency[Hz]');
    ylabel('hard : abs 2 fft/L');
    title('hard : fft after mic correlation')
    legend('0','1','2','3','4');
    hold off;
    subplot(2,1,2)
    hold on;
    plot(fp(1:length(IRoxc.mic0)),abs(fftshift(2/L*fft(IRoxc.mic0))));
    plot(fp(1:length(IRoxc.mic1)),abs(fftshift(2/L*fft(IRoxc.mic1))));
    plot(fp(1:length(IRoxc.mic2)),abs(fftshift(2/L*fft(IRoxc.mic2))));
    plot(fp(1:length(IRoxc.mic3)),abs(fftshift(2/L*fft(IRoxc.mic3))));
    plot(fp(1:length(IRoxc.mic4)),abs(fftshift(2/L*fft(IRoxc.mic4))));
    xlim([0 5000]);
    xlabel('frequency[Hz]');
    ylabel('abs 2 fft/L');
    title('open : fft after mic correlation')
    legend('0','1','2','3','4');
    hold off;
    
    % % %plot fft abs of time signal up to 100 miliseconds
    % i_tmax=5000; % around 100 miliseconds
    % L = i_tmax;
    % fp=(0:L-1)/L*fs-fs/2;
    % figure(400);
    % subplot(2,1,1)
    % hold on;
    % plot(fp(1:length(IRhxc.mic0(1:i_tmax))),abs(fftshift(2/L*fft(IRhxc.mic0(1:i_tmax)))));
    % plot(fp(1:length(IRhxc.mic0(1:i_tmax))),abs(fftshift(2/L*fft(IRhxc.mic1(1:i_tmax)))));
    % plot(fp(1:length(IRhxc.mic0(1:i_tmax))),abs(fftshift(2/L*fft(IRhxc.mic2(1:i_tmax)))));
    % plot(fp(1:length(IRhxc.mic0(1:i_tmax))),abs(fftshift(2/L*fft(IRhxc.mic3(1:i_tmax)))));
    % plot(fp(1:length(IRhxc.mic0(1:i_tmax))),abs(fftshift(2/L*fft(IRhxc.mic4(1:i_tmax)))));
    % xlim([0 5000]);
    % xlabel('frequency[Hz]');
    % ylabel('hard : abs 2 fft/L');
    % title('fft after mic correlation')
    % legend('0','1','2','3','4');
    % hold off;
    % subplot(2,1,2)
    % hold on;
    % plot(fp(1:length(IRhxc.mic0(1:i_tmax))),abs(fftshift(2/L*fft(IRoxc.mic0(1:i_tmax)))));
    % plot(fp(1:length(IRhxc.mic0(1:i_tmax))),abs(fftshift(2/L*fft(IRoxc.mic1(1:i_tmax)))));
    % plot(fp(1:length(IRhxc.mic0(1:i_tmax))),abs(fftshift(2/L*fft(IRoxc.mic2(1:i_tmax)))));
    % plot(fp(1:length(IRhxc.mic0(1:i_tmax))),abs(fftshift(2/L*fft(IRoxc.mic3(1:i_tmax)))));
    % plot(fp(1:length(IRhxc.mic0(1:i_tmax))),abs(fftshift(2/L*fft(IRoxc.mic4(1:i_tmax)))));
    % xlim([0 5000]);
    % xlabel('frequency[Hz]');
    % ylabel('abs 2 fft/L');
    % title('open : fft after mic correlation')
    % legend('0','1','2','3','4');
    % hold off;
end




%% FAST FOURRIER TRANSFORM AND COMPUTE TRANSFER FUNCTIONS H_0i
Hh = Fast_Fourier(IRhxc,fs);
Ho = Fast_Fourier(IRoxc,fs);

%% CALCULATIONS COEFFICIENTS, TRANSMISSION LOSS, CHARACTERISTIC IMPEDENCE
[TC, RC, AC, TL] = calc_coef_two_load(f, Hh, Ho, constants, measures);


%% PLOT FIGURES

if smooth == 1
    % %% SMOOTH DATA
    % TL_smooth = lowpass(TL, 1e-3);
    TC_smooth = lowpass(TC, 1e-3);
    RC_smooth = lowpass(RC, 1e-3);
    AC_smooth = lowpass(AC, 1e-3);

%     plot_coef(f, abs(TC_smooth), 'b', width, axis_coef, 1, 'Transmission coefficient [-]');
%     plot_coef(f, abs(TA), 'k--', width, axis_coef, 1, 'Transmission coefficient [-]');
%     plot_coef(f, abs(RC_smooth), 'r', width, axis_coef, 2, 'Reflection coefficient [-]');
%     plot_coef(f, abs(RA), 'k--', width, axis_coef, 2, 'Reflection coefficient [-]');
% %     plot_coef(f, abs(AC_smooth), line, width, axis_coef, 3, 'Absorption coefficient [-]');
% %     plot_coef(f, TL_smooth, line, width, axis_TL, 4, 'Transmission loss [dB]');
    
    %Transmission, Reflection  in 1 plot
    plot_coef(f, abs(TC_smooth), 'b--', width, axis_coef, 5, 'Coefficients');
    plot_coef(f, abs(TA), 'b', width, axis_coef, 5, 'Transmission coefficient [-]');
    plot_coef(f, abs(RC_smooth), 'r--', width, axis_coef, 5, 'Coefficients');
    plot_coef(f, abs(RA), 'r', width, axis_coef, 5, 'Reflection coefficient [-]');
    legend('Transmission', 'Transmission Analytic', 'Reflection', 'Reflection Analytic')
else
%     plot_coef(f, abs(TC), line, width, axis_coef, 1, 'Transmission coefficient [-]');
%     plot_coef(f, abs(TA), 'k--', width, axis_coef, 1, 'Transmission coefficient [-]');
%     plot_coef(f, abs(RC), line, width, axis_coef, 2, 'Reflection coefficient [-]');
%     plot_coef(f, abs(RA), 'k--', width, axis_coef, 2, 'Reflection coefficient [-]');
%     plot_coef(f, abs(AC), line, width, axis_coef, 3, 'Absorption coefficient [-]');
%     plot_coef(f, TL, line, width, axis_TL, 4, 'Transmission loss [dB]');
    
    %Transmission, Reflection  in 1 plot
    plot_coef(f, abs(TC), 'b--', width, axis_coef, 5, 'Coefficients');
    plot_coef(f, abs(TA), 'b', width, axis_coef, 5, 'Transmission coefficient [-]');
    plot_coef(f, abs(RC), 'r--', width, axis_coef, 5, 'Coefficients');
    plot_coef(f, abs(RA), 'r', width, axis_coef, 5, 'Reflection coefficient [-]');
    legend('Transmission', 'Transmission Analytic', 'Reflection', 'Reflection Analytic')
end



%% TODO list plotting
disp("TODO: tranform codeblock to generate synthetic data into a function");



% %% WRITE RESULTS TO .txt FILE
% 
% data = [f', TL_smooth', TC_smooth', RC_smooth', AC_smooth'];
% 
% 
% writematrix(data, strcat(Sample_path,'/',Sample,'_TwoLoad','_smooth','.txt'))
% if save_results==1
%     if smooth==1
%         data = [f', TL_smooth', TC_smooth', RC_smooth', AC_smooth'];
%         writematrix(data, strcat(Sample_path,'/',Sample,'_TwoLoad','_smooth','.txt'), 'Delimiter', 'space')
%     else 
%         data = [f', TL', TC', RC', AC'];
%         writematrix(data, strcat(Sample_path,'/',Sample,'_TwoLoad','_raw','.txt'), 'Delimiter', 'space')
%     end
% end
%     
%     
% %     %% Comparing COMSOL and experimental data
% %     
% %     figure(1)
% %     hold on;
% %     plot(f,TC);
% %     plot(cf,cT);
% %     hold off;
% %     
    