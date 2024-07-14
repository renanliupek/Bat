function [TC, RC, AC, TL] = calc_coef_two_load(f, Hh, Ho, constants, m)

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


%% FORMER FUNCTION
% Aparently there is an extra minus sign on each coeff denominador
%--------------------------------------------------------------------------
% for ff = 1:length(f)
% %Wavenumber
% k = 2*pi*f(ff)/constants.c;
% 
% %Calculation A, B, C and D matrices hard backing
% Ah = 1j * (Hh(f(ff)+1, 1)*exp(-1j*k*m(3))...
%          -Hh(f(ff)+1, 2)*exp(-1j*k*(m(3)+m(2))))...
%     /(2*sin(k*m(2)));                                    %Initial wave
% Bh = 1j * (Hh(f(ff)+1,2)*exp(1j*k*(m(3)+m(2)))...
%          -Hh(f(ff)+1,1)*exp(1j*k*m(3)))...
%     /(2*sin(k*m(2)));                                  %Reflected wave
% Ch = 1j * (Hh(f(ff)+1,3)*exp(1j*k*(m(4)+m(5)))...
%          -Hh(f(ff)+1,4)*exp(1j*k*m(4)))...
%     /(2*sin(k*m(5)));                                %Transmitted wave
% Dh = 1j * (Hh(f(ff)+1,4)*exp(-1j*k*m(4))...
%          -Hh(f(ff)+1,3)*exp(-1j*k*(m(4)+m(5))))...
%     /(2*sin(k*m(5)));                                  %Reflected wave
% 
% 
% %Calculation A, B, C and D matrices open backing
% Ao = 1j * (Ho(f(ff)+1, 1)*exp(-1j*k*m(3))...
%          -Ho(f(ff)+1, 2)*exp(-1j*k*(m(3)+m(2))))...
%     /(2*sin(k*m(2)));                                    %Initial wave
% Bo = 1j * (Ho(f(ff)+1,2)*exp(1j*k*(m(3)+m(2)))...
%          -Ho(f(ff)+1,1)*exp(1j*k*m(3)))...
%     /(2*sin(k*m(2)));                                  %Reflected wave
% Co = 1j * (Ho(f(ff)+1,3)*exp(1j*k*(m(4)+m(5)))...
%          -Ho(f(ff)+1,4)*exp(1j*k*m(4)))...
%     /(2*sin(k*m(5)));                                %Transmitted wave
% Do = 1j * (Ho(f(ff)+1,4)*exp(-1j*k*m(4))...
%          -Ho(f(ff)+1,3)*exp(-1j*k*(m(4)+m(5))))...
%     /(2*sin(k*m(5)));                                  %Reflected wave
% 
% 
% % Calculation coefficients
% TC(ff) = (Ch*Do-Co*Dh)/(Ah*Do-Ao*Dh); % Transmission Coefficient [-]
% RC(ff) = (Bh*Do-Bo*Dh)/(Ah*Do-Ao*Dh); % Reflection coefficient  
% AC(ff) = 1-(abs(RC(ff))).^2-(abs(TC(ff))).^2; %Absorbtion coefficient
% 
% % Calculation Transmission Loss
% TL(ff) = 20*log10(abs(1/TC(ff)));   % Transmission loss [dB]

% end
%--------------------------------------------------------------------------

%% CURRENT FUNCTION
%--------------------------------------------------------------------------
% RECAP
    %     Coordinates of each microphone 0, 1, 2, 3, and 4.
    %     (note: left sample surface is at x=0, right, x=d)
    %     mic0.x =-s0-s1-l1;
    %     mic1.x =-s1-l1;
    %     mic2.x =-l1;
    %     mic3.x = l2;
    %     mic4.x = l2+s2;
    %
    %     m = [s0, s1, l1, l2, s2, ze];
    %
    %     Equivalent relations to compare with "Technical review..."
    %     [Bolton,2007]
    %
    %     x0 = -m(1)-m(2)-m(3);
    %     x1 = -m(2)-m(3);
    %     x2 = -m(3);
    %  x1-x2 = -m(2);
    %     x3 =  m(4);
    %     x4 =  m(4)+m(5);
    %  x3-x4 = -m(5);
%--------------------------------------------------------------------------

countnoise=0;
%--------------------------------------------------------------------------
for ff = 1:length(f) %loop on freqs you wish to compute the ac indicators
    
    % define wavenumber
    k = 2*pi*f(ff)/constants.c;

    %Calculation A, B, C and D matrices hard backing
    Ah = 1j * (Hh(f(ff)+1, 1)*exp( 1j*k*(-m(3))      )...
               -Hh(f(ff)+1, 2)*exp( 1j*k*(-m(3)-m(2)) )   )...
                /(    2*sin( k * (-m(2)) )   );
    
    Bh = 1j * (Hh(f(ff)+1, 2)*exp(-1j*k*(-m(3)-m(2)) )...
               -Hh(f(ff)+1, 1)*exp(-1j*k*(-m(3))      )   )...
                /(    2*sin( k * (-m(2)) )   );
    
    Ch = 1j * (Hh(f(ff)+1,3)*exp( 1j*k*(m(4)+m(5))  )...
               -Hh(f(ff)+1,4)*exp( 1j*k* m(4)        )   )...
               /(    2*sin( k * (-m(5))  )    ); 
    
    Dh = 1j * (Hh(f(ff)+1,4)*exp(-1j*k* m(4)        )...
               -Hh(f(ff)+1,3)*exp(-1j*k* (m(4)+m(5)) )   )...
               /(    2*sin( k * (-m(5))  )    ); 
  
    
    %Calculation A, B, C and D matrices open backing
    Ao = 1j * (Ho(f(ff)+1, 1)*exp( 1j*k*(-m(3))      )...
               -Ho(f(ff)+1, 2)*exp( 1j*k*(-m(3)-m(2)) )   )...
                /(    2*sin( k * (-m(2)) )   );
    
    Bo = 1j * (Ho(f(ff)+1, 2)*exp(-1j*k*(-m(3)-m(2)) )...
               -Ho(f(ff)+1, 1)*exp(-1j*k*(-m(3))      )   )...
                /(    2*sin( k * (-m(2)) )   );
    
    Co = 1j * (Ho(f(ff)+1,3)*exp( 1j*k*(m(4)+m(5))  )...
               -Ho(f(ff)+1,4)*exp( 1j*k* m(4)        )   )...
               /(    2*sin( k * (-m(5))  )    ); 
    
    Do = 1j * (Ho(f(ff)+1,4)*exp(-1j*k* m(4)        )...
               -Ho(f(ff)+1,3)*exp(-1j*k* (m(4)+m(5)) )   )...
               /(    2*sin( k * (-m(5))  )    ); 
    

    % Compute coefficients R and T
    if(abs(Ah*Do-Ao*Dh)>1e-10) % if denominator isnt noise
    TC(ff) = (Ch*Do-Co*Dh)/(Ah*Do-Ao*Dh); % Transmission Coefficient [-]
    RC(ff) = (Bh*Do-Bo*Dh)/(Ah*Do-Ao*Dh); % Reflection coefficient
    else
     countnoise=countnoise+1;    
     %fprintf("noise detected for frequency %f\n", ff);
     TC(ff) = (Ch*Do-Co*Dh); 
     RC(ff) = (Bh*Do-Bo*Dh);
    end
    
    % Compute absorption coefficient
    AC(ff) = 1-(abs(RC(ff))).^2-(abs(TC(ff))).^2; 
    
    % Compute transmission loss [dB]
    TL(ff) = 20*log10(abs(1/TC(ff)));   

end
%--------------------------------------------------------------------------

fprintf(['calc_coef_two_load called\n      number of noise detected: %d out' ...
         ' of %d NFFT data points\n\n'], countnoise, length(f));


end