%% README
% Post processing script of experimental impedance tube. Adapted from
% MSc L. Schijff, Ashoka karunarathne, N.A. van de Straat, 
% by R. Liupekevicius 20-12-2021. Eindhoven University of Technology
%
% INSTRUCTIONS
% Main settings are commented by a sandwitch of '%%%%' with a short
% description. Please see the example below.
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
%% TWO LOAD METHOD POSTPROCESSING
% ------------------------------------------------------------------------
% clear all   % clear all variables in the workspace
% close all   % close figures. Comment if multiple samples need to be plotted
clc         % clears the command window

%% ADD FUNCTIONS TO PATH
% % add the shared functions to the search path
% % path(path, ['shared_functions']);

%% ADD DATA FOLDER TO PATH
% ------------------------------------------------------------------------
% Check and change these parameters for every run
% ------------------------------------------------------------------------
% Data directory specification
% 'Data_path' + 'Folder' + 'Sample' should form the complete directory to the
% file with the mic data

% Path to general data folder 
% (pwd = built in command for current folder path)
Data_path = fullfile(pwd,'data');  
                      


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% choose experimental data to be post processed %%%%%%%%%%%%%%%%%
Sample    = 'nosample1';   %Sample folder that contains audio files    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Complete path to desired data set
Sample_path = fullfile(Data_path,Sample);     
   


%% PLOT SETTINGS
% Plot settings; line color, smoot/raw data, axis ranges       

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Plotting of smooth (=1) or raw(=0) data %%%%%%%%%%%%%%%%%%%%%%
smooth      = 1;             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% line_type        = '--'; 
line_type        = ':'; 
% line_type        = '-.'; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Line width
width       = 1.5;        
%Axis range transmission loss - 
% [x-axis lower limit, x-axis upper limit,...
% y-axis lower limit, y-axis upper limit]
axis_TL     = [100 5000 -1 15];     
%Axis range coefficients - 
% [x-axis lower limit, x-axis upper limit,...
% y-axis lower limit, y-axis upper limit]
axis_coef   = [100 5000 0 1];    

save_results = 1;   %Indicater to save data to .txt file yes=1 and no=0
                    %The saved data is structured as:
                    %column 1: Frequency [Hz]
                    %column 2: Transmission loss [dB]
                    %column 3: Transmission coefficient [-]
                    %column 4: Reflection coefficient [-]
                    %column 5: Absorption coefficient [-]

%% MEASUREMENT SPECS
% Measurement specifications

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% sample thickness [m] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d     = 0.1; % [m] should be smaller than l2      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% temperature during measurement [K]
temp  = 293.15;     

% frequency range, planewave up to 1.84*c/(2*pi*r)
f     = 1:1:5e3;    %Frequency range [Hz]

% impedance tube radius (used for synthetic data generation, open bc)
r     = (40e-3)/2 ; %[m]
                    
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

%% ANALYTICAL R AND T FOR HOMOGENEOUS/ DELANY-BAZLEY LAYER
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% choose the model to compare your sample %%%%%%%%%%%%%%%%%%%%%
% ref_model='model-off';     % no actions required      
% ref_model='homogeneous';     % must define rho_L and c_L   
ref_model='Delany-Bazley'; % must define sigma
% ref_model='Infinite-Plate'; % must define sigma
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch ref_model
    
    case 'model-off'
        
        TA = nan*zeros(1,length(f));
        RA = nan*zeros(1,length(f));

    case 'homogeneous'
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%% define homogeneous layer parameters rho_L and c_L %%%%%%%%
        constants.rho_L = 145; % [kg/m^3] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        constants.c_L   = 250; %[ m/s] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % characteristic impedance and wavenumber
        ZL              = constants.rho_L * constants.c_L;
        a               = 2*pi*f/constants.c_L;
        
        % analytical reflection and transmission coefficients for
        % homogeneous layer check analytical_T_and_R.nb file for more.
        RA = ((-Z0^2 + ZL.^2).*sin(a.*d))./(-2*1i*Z0*ZL.*cos(a.*d) + ...
              (Z0^2 + ZL.^2).*sin(a.*d));
        TA = (4*exp(1i*d*(a + k)).*Z0.*ZL)./(-(Z0 - ZL).^2 + ...
              exp(2*1i*a.*d).*(Z0 + ZL).^2);    

    case 'Delany-Bazley'
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        constants.sigma = 30000; % flow resitivity %%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        sigma=constants.sigma; 
        X  = constants.rho*f/sigma;
        ZL = Z0*(1+0.057*X.^(-0.754)-1i*0.087.*X.^(-0.732)) ; 
        a  = k.*(1+0.0978*X.^(-0.700)-1i*0.189.*X.^(-0.595)); 
        % Note that the real part of the surface impedance is negative
        % at low frequencies, which is not physical. This defect in the
        % Delany and Bazley model has been noticed by Miki (1990) who 
        % proposed modifications to Equations 4.102 from "finite element
        % and boundary..." [N.Atalla, 2015].

        % same as analytical reflection and transmission coefficients for
        % homogeneous layer, however, plugging in themodified
        % charactertistic impedance and wavenumber.
        RA = ((-Z0^2 + ZL.^2).*sin(a.*d))./(-2*1i*Z0*ZL.*cos(a.*d) + ...
              (Z0^2 + ZL.^2).*sin(a.*d));
        TA = (4*exp(1i*d*(a + k)).*Z0.*ZL)./(-(Z0 - ZL).^2 + ...
              exp(2*1i*a.*d).*(Z0 + ZL).^2);    

end

        
%% TIME VECTOR, SAMPLING FREQUENCY fs
fs   = 48000;
time = (0:fs-1)*1/fs; %timestep of 1/fs

%% LOADING + CUT MIC DATA
% Load full data files hard & open backing
[fs, IRmich] = Load_mic_data(Sample_path, 'Hard');
[~, IRmico]  = Load_mic_data(Sample_path, 'Open');

%Cut data to first second
for ii = 1:length(IRmich(1,:))
    IRh(:,ii) = IRmich(1:fs, ii);
    IRo(:,ii) = IRmico(1:fs, ii);
end

%% PLOT TIME SIGNAL OF MICS, termination hard and open
%     figure(600)
%     subplot(2,1,1)
%     hold on;
%     plot(1000*time,IRh(:,1));
%     plot(1000*time,IRh(:,2));
%     plot(1000*time,IRh(:,3));
%     plot(1000*time,IRh(:,4));
%     plot(1000*time,IRh(:,5));
%     % xlim([0 1000*5*1/f(ind)]); %5 periods
%     xlabel("time(miliseconds)");
%     ylabel("pressure[Pa]");
%     title("hard tube ending (closed)");
%     hold off;
%     grid on;
%     legend('0','1','2','3','4');
%     subplot(2,1,2)
%     hold on;
%     plot(1000*time,IRo(:,1));
%     plot(1000*time,IRo(:,2));
%     plot(1000*time,IRo(:,3));
%     plot(1000*time,IRo(:,4));
%     plot(1000*time,IRo(:,5));
%     % xlim([0 1000*5*1/f(ind)]); %5 periods
%     xlabel("time(miliseconds)");
%     ylabel("pressure[Pa]");
%     title("open tube ending");
%     hold off;
%     grid on;
%     legend('0','1','2','3','4');

%% Plot FFT of raw time signals

%     % plot abs fft of microphones
%     L=fs; %length of time signal
%     fp=(0:L-1)/L*fs-fs/2;
%     figure(300);
%     subplot(2,1,1);
%     hold on;
%     plot(fp,abs(fftshift(2/L*fft(IRh(:,1)))));
%     plot(fp,abs(fftshift(2/L*fft(IRh(:,2)))));
%     plot(fp,abs(fftshift(2/L*fft(IRh(:,3)))));
%     plot(fp,abs(fftshift(2/L*fft(IRh(:,4)))));
%     plot(fp,abs(fftshift(2/L*fft(IRh(:,5)))));
%     xlim([0 5000]);
%     xlabel('frequency[Hz]');
%     ylabel('hard : abs 2 fft/L')
%     legend('0','1','2','3','4');
%     title('hard : fft before mic correlation');
%     hold off;
%     subplot(2,1,2);
%     hold on;
%     plot(fp,abs(fftshift(2/L*fft(IRo(:,1)))));
%     plot(fp,abs(fftshift(2/L*fft(IRo(:,2)))));
%     plot(fp,abs(fftshift(2/L*fft(IRo(:,3)))));
%     plot(fp,abs(fftshift(2/L*fft(IRo(:,4)))));
%     plot(fp,abs(fftshift(2/L*fft(IRo(:,5)))));
%     xlim([0 5000]);
%     xlabel('frequency[Hz]');
%     ylabel('open : abs 2 fft/L');
%     legend('0','1','2','3','4');
%     title('open : fft before mic correlation');
%     hold off;

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


%Cross correlation
IRhxc = Cross_correlation(IRh, N, x, constants, fs);
IRoxc = Cross_correlation(IRo, N, x, constants, fs);


%% FAST FOURRIER TRANSFORM AND COMPUTE TRANSFER FUNCTIONS H_0i
Hh = Fast_Fourier(IRhxc,fs);
Ho = Fast_Fourier(IRoxc,fs);

%% CALCULATIONS COEFFICIENTS, TRANSMISSION LOSS, CHARACTERISTIC IMPEDENCE
[TC, RC, AC, TL] = calc_coef_two_load(f, Hh, Ho, constants, measures);


%% PLOT FIGURES
switch ref_model
    case 'model-off'
         if smooth == 1
            %SMOOTH DATA
            TL_smooth = lowpass(TL, 1e-3);
            TC_smooth = lowpass(TC, 1e-3);
            RC_smooth = lowpass(RC, 1e-3);
            AC_smooth = lowpass(AC, 1e-3);
            
        
            %TRANSMISSION AND REFLECTION COEFFICIENTS
            plot_coef(f, abs(TC_smooth), 'b', line_type, width, axis_coef, 5);
            plot_coef(f, abs(RC_smooth), 'r', line_type, width, axis_coef, 5);     
            ylabel('Coefficients[-]');
            if exist('legendc', 'var')
                legendc=[legendc {['T ' Sample],['R' Sample]}];
            else
            legendc={['T' Sample],['R ' Sample]};
            end
            legend(legendc)
            title('Smoothened')
            savefig(gcf,[Sample_path '\_smooth_coeffs']);
        
        
        
            % ABSORPTION COEFFICIENT
            plot_coef(f, AC_smooth, 'k', line_type, width, axis_coef, 6);
            title('Smoothened')
            if exist('legenda', 'var') 
                legenda=[legenda {['A ' Sample]}];
            else
            legenda={['A ' Sample]};
            end
            legend(legenda)
            savefig(gcf,[Sample_path '\_smooth_absorp']);
        
        
          
           % TRANSMISSION LOSS
            plot_coef(f, TL_smooth, 'g', line_type, width, axis_TL, 7);
            title('Smoothened')
            if exist('legendtl', 'var') 
                legendtl=[legendtl {['TL ' Sample]}];
            else
            legendtl={['TL ' Sample]};
            end
            legend(legendtl);
            savefig(gcf,[Sample_path '\_smooth_TL']);
        
         else
            error('No plot available for nonsmooth option. Ask R.Liupekevicius for update.')
         end

    otherwise
        if smooth == 1
            %SMOOTH DATA
            TL_smooth = lowpass(TL, 1e-3);
            TC_smooth = lowpass(TC, 1e-3);
            RC_smooth = lowpass(RC, 1e-3);
            AC_smooth = lowpass(AC, 1e-3);
            
        
            %TRANSMISSION AND REFLECTION COEFFICIENTS
            plot_coef(f, abs(TC_smooth), 'b', line_type, width, axis_coef, 5);
            plot_coef(f, abs(TA),        'b', '-',    width, axis_coef, 5);
            plot_coef(f, abs(RC_smooth), 'r', line_type, width, axis_coef, 5);
            plot_coef(f, abs(RA),        'r', '-',    width, axis_coef, 5);
            ylabel('Coefficients[-]');
            if exist('legendc', 'var')
                legendc=[legendc {['T ' Sample], ['T ' ref_model],...
                                  ['R ' Sample], ['R ' ref_model]}];
            else
            legendc={ ['T' Sample], ['T ' ref_model],...
                      ['R' Sample], ['R ' ref_model]};
            end
            legend(legendc)
            title('Smoothened')
            savefig(gcf,[Sample_path '\_smooth_coeffs']);
        
        
        
            % ABSORPTION COEFFICIENT
            plot_coef(f, AC_smooth, 'k', line_type, width, axis_coef, 6);
            title('Smoothened')
            if exist('legenda', 'var') 
                legenda=[legenda {['A ' Sample]}];
            else
            legenda={['A ' Sample]};
            end
            legend(legenda)
            savefig(gcf,[Sample_path '\_smooth_absorp']);
        
        
          
           % TRANSMISSION LOSS
            plot_coef(f, TL_smooth, 'g', line_type, width, axis_TL, 7);
            plot_coef(f, 20*log10(1./TA), 'g', '-', width, axis_TL, 7);
            title('Smoothened')
            if exist('legendtl', 'var') 
                legendtl=[legendtl {['TL ' Sample],['TL ' ref_model]}];
            else
            legendtl={['TL ' Sample],['TL ' ref_model]};
            end
            legend(legendtl);
            savefig(gcf,[Sample_path '\_smooth_TL']);
        
        else
            
           error('No plot available for nonsmooth option. Ask R.Liupekevicius for update.')
        end
end


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
    