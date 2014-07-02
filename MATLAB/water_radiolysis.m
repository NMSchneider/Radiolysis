function products = water_radiolysis(t,reactants,doseRate,gValues)
% Water Radiolysis Model
%{
    The MIT License (MIT)
    
    Copyright (c) 2014 Brian J. Mendel, Nicholas M. Schneider
    
    Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions: The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
        
        THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
        
%}



% This function, when used in conjunction with an ODE solver (ode23s and
% ode15s work best).
%
% Concentrations are in MICROMOLAR, and time is in SECONDS.
% 
% 
%"t" is the time span [tinitial tfinal] in SECONDS, "reactants" is a row  
% vector of 16 elements specifying the initial concentrations of radiolysis
% products in MICROMOLAR.




%equilibria%
K = zeros(5,1);

%H2O <=> H+ + OH-&
K(1) = 10^-13.999; %Molar^2

%H2O2 <=> H+ + HO2-%
K(2) = 10^-11.65;%Molar

%OH <=> H+ + O-%
K(3) = 10^-11.9;%Molar

%HO2 <=> H+ + O2-%
K(4) = 10^-4.57;%Molar

%H <=> H+ + e-%
K(5) = 10^-9.77;%Molar

%Rate constants%
%%1/(Ms) unless specified%%

%Molar to micromolar conversion factor%
molarToMicromolar = 10^6;

% Note, some are out of numerical order so they are defined when they are
% used in other rate constants. To convert numbering scheme from what is
% presented here to what it presented in the paper [ref], add 6 to the
% index. (i.e., Reaction 25's rate constant has an index of 19)
k = zeros(73,1);
k(1) = 1.4*10^11;
k(2) = k(1)*K(1)*molarToMicromolar^2;%M/s
k(4) = 5*10^10;
k(3) = k(4)*K(2)*molarToMicromolar;%/s
k(5) = 1.3*10^10;
k(6) = k(5)*(K(1)/K(2))*molarToMicromolar;%%/s
k(7) = 1.9*10;
k(8) = 2.2*10^7;
k(10) = 2.3*10^10;
k(9) = k(10)*K(5)*molarToMicromolar;%/s
k(11) = 1.3*10^10;
k(12) = k(11)*(K(1)/K(3))*molarToMicromolar;%/s
k(14) = 10^11;
k(13) = k(14)*K(3)*molarToMicromolar;%/s
k(16) = 5*10^10;
k(15) = k(16)*K(4)*molarToMicromolar;%/s
k(17) = 5*10^10;
k(18) = k(17)*(K(1)/K(4))*molarToMicromolar;%/s%
k(19) = 3*10^10;
k(20) = 1.1*10^10;
k(21) = 1.3*10^10 ;
k(22) = 2*10^10;
k(23) = 1.9*10^10;
k(24) = 5.5*10^9; 
k(25) = 2.5*10^10 ;
k(26) = 3.5*10^9;
k(27) = 2.2*10^10 ;
k(28) = 1.6*10^10; 
k(29) = 3.6*10^10;
k(30) = 1.1*10;
k(31) = 10^10;
k(32) = 9*10^7;
k(33) = 10^10;
k(34) = 7.8*10^9;
k(35) = 7.0*10^9;
k(36) = 9*10^7;
k(37) = 2.1*10^10;
k(38) = 1.8*10^10;
k(39) = 1.8*10^10;
k(40) = 3.8*10^10;
k(41) = 3.6*10^9;
k(42) = 6*10^9; 
k(43) = 8.2*10^9;
k(44) = 4.3*10^7;
k(45) = 2.7*10^7;
k(46) = 2.5*10^10;
k(47) = 7.5*10^9;
k(48) = 2.6*10^9;
k(49) = 6*10^9;
k(50) = 1.1*10^8;
k(51) = 8*10^7;
k(52) = 7*10^5;
k(53) = 6*10^9;
k(54) = 5*10^-1;
k(55) = 5*10^-1;
k(56) = 6*10^9;
k(57) = 5*10^8;
k(58) = 10^2;
k(59) = 6*10^8;
k(60) = 1.3*10^-1;
k(61) = 1.3*10^-1;
k(62) = 10^4 ;
k(63) = 1.5*10^9;
k(64) = 10^9 ;
k(65) = 3.6*10^9;
k(66) = 8*10^7;
k(67) = 5*10^8;
k(68) = 4*10^8;
k(69) = 7*10^8;
k(70) = 5*10^9;
k(71) = 3.3*10^3*molarToMicromolar;%/s
k(72) = 9*10^10;
k(73) = 1.1*10^5*molarToMicromolar;%/s


%%convert from molar to micromolar
k = k./(molarToMicromolar);

%avogadro's number%
avogadro = 6.022*10^23;

%convert from molecules/(100 eV) to molecules/eV
% gValues = gValues./100;

%conversion from cubic meters to liters%
m3_L = 1000;

%%calculated in micromolar/s
doseRate = doseRate / (1.6e-22);
generationDueToGValues = molarToMicromolar * doseRate.*gValues./(avogadro*m3_L);

%Reactant concentrations (micromolar)%
electrons = reactants(1); 
H_ions    = reactants(2);
OH_ions   = reactants(3);
H2O2      = reactants(4);
HO2_ions  = reactants(5);
H_mono    = reactants(6);
OH        = reactants(7);
O_ions    = reactants(8);
HO2       = reactants(9);
O2_ions   = reactants(10);
O2        = reactants(11);
H2        = reactants(12);
O3_ions   = reactants(13);
O3        = reactants(14);
HO3       = reactants(15);
H2O       = reactants(16);


%Reaction Set%
% Note: H2O is divided out as a reactant
r = zeros(73,1);
r(1)  = k(1)*H_ions*OH_ions;
r(2)  = k(2);
r(3)  = k(3)*H2O2;
r(4)  = k(4)*H_ions*HO2_ions;
r(5)  = k(5)*H2O2*OH_ions;
r(6)  = k(6)*HO2_ions;
r(7)  = k(7)*electrons*H2O;
r(8)  = k(8)*H_mono*OH_ions;
r(9)  = k(9)*H_mono;
r(10) = k(10)*electrons*H_ions;
r(11) = k(11)*OH*OH_ions;
r(12) = k(12)*O_ions;
r(13) = k(13)*OH;
r(14) = k(14)*O_ions*H_ions;
r(15) = k(15)*HO2;
r(16) = k(16)*O2_ions*H_ions;
r(17) = k(17)*HO2*OH_ions;
r(18) = k(18)*O2_ions;
r(19) = k(19)*electrons*OH;
r(20) = k(20)*electrons*H2O2;
r(21) = k(21)*electrons*O2_ions;
r(22) = k(22)*electrons*HO2;
r(23) = k(23)*electrons*O2;
r(24) = k(24)*electrons^2;
r(25) = k(25)*electrons*H_mono;
r(26) = k(26)*electrons*HO2_ions;
r(27) = k(27)*electrons*O_ions;
r(28) = k(28)*electrons*O3_ions;
r(29) = k(29)*electrons*O3;
r(30) = k(30)*H_mono*H2O; 
r(31) = k(31)*H_mono*O_ions;
r(32) = k(32)*H_mono*HO2_ions;
r(33) = k(33)*H_mono*O3_ions;
r(34) = k(34)*H_mono^2;
r(35) = k(35)*H_mono*OH;
r(36) = k(36)*H_mono*H2O2;
r(37) = k(37)*H_mono*O2;
r(38) = k(38)*H_mono*HO2;
r(39) = k(39)*H_mono*O2_ions;
r(40) = k(40)*H_mono*O3;
r(41) = k(41)*OH^2;
r(42) = k(42)*OH*HO2; 
r(43) = k(43)*OH*O2_ions;
r(44) = k(44)*OH*H2;
r(45) = k(45)*OH*H2O2;
r(46) = k(46)*OH*O_ions;
r(47) = k(47)*OH*HO2_ions;
r(48) = k(48)*OH*O3_ions;
r(49) = k(49)*OH*O3_ions;
r(50) = k(50)*OH*O3;
r(51) = k(51)*HO2*O2_ions;
r(52) = k(52)*HO2^2;
r(53) = k(53)*HO2*O_ions;
r(54) = k(54)*HO2*H2O2;
r(55) = k(55)*HO2*HO2_ions;
r(56) = k(56)*HO2*O3_ions;
r(57) = k(57)*HO2*O3;
r(58) = k(58)*O2_ions^2;
r(59) = k(59)*O2_ions*O_ions;
r(60) = k(60)*O2_ions*H2O2;
r(61) = k(61)*O2_ions*HO2_ions;
r(62) = k(62)*O2_ions*O3_ions;
r(63) = k(63)*O2_ions*O3;
r(64) = k(64)*O_ions^2;
r(65) = k(65)*O_ions*O2;
r(66) = k(66)*O_ions*H2;
r(67) = k(67)*O_ions*H2O2;
r(68) = k(68)*O_ions*HO2_ions;
r(69) = k(69)*O_ions*O3_ions;
r(70) = k(70)*O_ions*O3;
r(71) = k(71)*O3_ions;
r(72) = k(72)*O3_ions*H_ions;
r(73) = k(73)*HO3;

products = zeros(length(reactants),1);

%%electrons%%
products(1) = -r(7) + r(8) + r(9) - r(10) - r(19) - r(20) - r(21) - r(22)...
    -r(23) - 2*r(24) - r(25) - r(26) - r(27) - r(28) - r(29);

%%H+%%
products(2) = -r(1) + r(2) + r(3) - r(4) + r(9) - r(10) + r(13) - r(14)...
    + r(15) - r(16) + r(49) - r(72);

%%OH-%%
products(3) = -r(1) + r(2) - r(5) + r(6) + r(7) - r(8) - r(11) + r(12) -...
    r(17)+ r(18) + r(19) + r(20) + r(21) + 2*r(24) + r(25) + r(26) + ...
    2*r(27) + 2*r(28) + r(31) + r(32) + r(33) + r(43) + r(47) + r(48) + ...
    r(53) + r(55) + r(56) + 2*r(58) + 2*r(59) + r(60) + r(61) + 2*r(62)...
    + r(64) + r(66) + r(68);

%%H2O2%%
products(4) = - r(3) + r(4) -r(5) + r(6) - r(20) - r(36) + r(38) + r(41)...
    - r(45) + r(52) - r(54) + r(58) - r(60) - r(67)   ;

%%HO2-%%
products(5) = r(3) - r(4) + r(5) - r(6) + r(21) + r(22) - r(26) - r(32)...
    + r(39) + r(46) - r(47) + r(51) - r(55) - r(61) + r(64) - r(68);

%%H%%
products(6) = r(7) - r(8) - r(9) + r(10) - r(25) - r(30) - r(31) - r(32)...
    - r(33) - 2*r(34) - r(35) - r(36) - r(37) - r(38) - r(39) - r(40) + r(44)...
    + r(66);

%%OH%%
products(7) = -r(11) + r(12) - r(13) + r(14) - r(19) + r(20) + r(30) ...
    + r(32) - r(35) + r(36) - 2*r(41) - r(42) - r(43) - r(44) - r(45) -...
    r(46) - r(47) - r(48) - r(49) - r(50) + r(54) + r(55) + r(60) + r(72)...
    + r(73);

%%O-%%
products(8) = r(11) - r(12) + r(13) - r(14) + r(26) - r(27) - r(31) - r(46)...
    - r(53) - r(59) + r(61) - 2*r(64) - r(65) - r(66) - r(67) - r(68) -...
    r(69) - r(70) + r(71) ;

%%HO2%%
products(9) = -r(15) + r(16) - r(17) + r(18) - r(22) + r(37) - r(38) - ...
    r(42) + r(45) + r(47) + r(50) - r(51) - 2*r(52) - r(53) - r(54) - r(55)...
    - r(56) - r(57);

%%O2-%%
products(10) = r(15) - r(16) + r(17) - r(18) - r(21) + r(23) - r(39) - ...
    r(43) + 2*r(49) - r(51) - 2*r(58) - r(59) - r(60) - r(61) - r(62) - r(63)...
    + r(67) + r(68) + 2*r(69) + r(70);

%%O2%%
products(11) = -r(23) + r(28) + r(33) - r(37) + r(42) + r(43) + r(50) ...
    + r(51) + r(52) + r(53) + r(54) + r(55) + 2*r(56) + r(57) + r(58) + r(59)...
    + r(60) + r(61) + 2*r(62) + r(63) - r(65) + r(70) + r(71) + r(72) + r(73);

%%H2%%
products(12) = r(24) + r(25) + r(30) + r(34) - r(44) - r(66);

%%O3-%%
products(13) = -r(28) + r(29) - r(33) - r(48) - r(49) - r(56) - r(62) ...
    + r(63) + r(65) - r(69) - r(71) - r(72);

%%O3%%
products(14) = -r(29) - r(40) + r(48) - r(50) - r(57) - r(63) - r(70);

%%HO3%%
products(15) = r(40) + r(57) - r(73);

%%H2O%%
products(16) = r(1) - r(2) + r(5) - r(6) - r(7) + r(8) + r(11) - r(12) + ...
     r(17) - r(18) - r(21) - 2*r(24) - r(25) - r(27) - r(28) - r(30) + r(35)...
     + r(36) + r(42) + r(44) + r(45) + r(54) - 2*r(58) - r(59) - r(62) - r(64) ...
     + r(67);

products = products + generationDueToGValues;

