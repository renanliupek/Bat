function [TC, RC, AC, TL, TM1, TM2] = calc_TM(f, Hh, Ho, constants, m, d)

%%------------------------------------------------------------------------
% INPUT
% f = frequency vector
% H = Transfer functions between mic 0 and mic 1-4
% constants = struct with constants used
% d = thickness sample
% m = vector with measures of the tube [s0, s1, l1, l2, s2, ze];
% 
% OUTPUT
% TC = Transmission coefficient
% RC = Reflection coefficient
% AC = Absorption coefficient
% TL = Transmission Loss
% Z  = Characteristic Impedence
%% ---------------------------------------------------------------------

for ff = 1:length(f)
    %Wavenumber
    k = 2*pi*f(ff)/constants.c;
    
    %Calculation A, B, C and D matrices hard backing
    Ah = 1j * (Hh(f(ff)+1, 1)*exp(-1j*k*m(3))...
             -Hh(f(ff)+1, 2)*exp(-1j*k*(m(3)+m(2))))...
        /(2*sin(k*m(2)))   ;                                %Initial wave
    Bh = 1j * (Hh(f(ff)+1,2)*exp(1j*k*(m(3)+m(2)))...
             -Hh(f(ff)+1,1)*exp(1j*k*m(3)))...
        /(2*sin(k*m(2)))  ;                              %Reflected wave
    Ch = 1j * (Hh(f(ff)+1,3)*exp(1j*k*(m(4)+m(5)))...
             -Hh(f(ff)+1,4)*exp(1j*k*m(4)))...
        /(2*sin(k*m(5)));                                %Transmitted wave
    Dh = 1j * (Hh(f(ff)+1,4)*exp(-1j*k*m(4))...
             -Hh(f(ff)+1,3)*exp(-1j*k*(m(4)+m(5))))...
        /(2*sin(k*m(5)));                                  %Reflected wave

    
    %Calculation A, B, C and D matrices open backing
    Ao = 1j * (Ho(f(ff)+1, 1)*exp(-1j*k*m(3))...
             -Ho(f(ff)+1, 2)*exp(-1j*k*(m(3)+m(2))))...
        /(2*sin(k*m(2)));                                    %Initial wave
    Bo = 1j * (Ho(f(ff)+1,2)*exp(1j*k*(m(3)+m(2)))...
             -Ho(f(ff)+1,1)*exp(1j*k*m(3)))...
        /(2*sin(k*m(2)));                                  %Reflected wave
    Co = 1j * (Ho(f(ff)+1,3)*exp(1j*k*(m(4)+m(5)))...
             -Ho(f(ff)+1,4)*exp(1j*k*m(4)))...
        /(2*sin(k*m(5)));                                %Transmitted wave
    Do = 1j * (Ho(f(ff)+1,4)*exp(-1j*k*m(4))...
             -Ho(f(ff)+1,3)*exp(-1j*k*(m(4)+m(5))))...
        /(2*sin(k*m(5)));                                  %Reflected wave
    

    % Calculation coefficients
    TC(ff) = 1/((Ah*Do-Ao*Dh)/(Ch*Do-Co*Dh)); % Transmission Coefficient [-]
    
    RC(ff) = (Bh*Do-Bo*Dh)/(Ah*Do-Ao*Dh); %Reflection coefficient
    
    AC(ff) = 1-(abs(RC(ff))).^2-(abs(TC(ff))).^2; %Absorbtion coefficient
    
    % Calculation Transmission Loss
    TL(ff) = 20*log10(abs(1/TC(ff)));   % Transmission loss [dB]
    

    % Calculation od transfer matrix METHOD 1
    
        % matrix connecting pressure and velocity to wave amplitudes
        % [p v]^T = [E] [A B]^T or [p v]^T = [E] [C D]^T
        rho =constants.rho;
        c   =constants.c;
        E = @(x) [exp(-1i*k*x)           exp(1i*k*x)     ;
                  exp(-1i*k*x)/rho/c    -exp(1i*k*x)/rho/c];

        % matrix L connects amplitudes of each sample size
        % [A B]^T = [L] [C D]^T 
        L = [Ah Ao; Bh Bo] * inv([Ch Co; Dh Do]);

        % Trnsfer matrix calculation
        TM1(:,:,ff) = E(0)*L*inv(E(d));
        TM2(:,:,ff) = E(-d)*L*inv(E(0));

    % Calculation od transfer matrix METHOD 2 (B&K report)




    
end

end