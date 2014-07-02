% Water Radiolysis Model 
%{
    The MIT License (MIT)

    Copyright (c) 2014 Nicholas M. Schneider

    Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions: The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
        
%}

% Species Tracking Information. This is the mapping used to convert indexs
% to species names: 
% Index : Species
%  1 : eh- (hydrated electrons)
%  2 : H+
%  3 : OH-
%  4 : H2O2
%  5 : HO2
%  6 : H
%  7 : OH
%  8 : O-
%  9 : HO2
% 10 : O2-
% 11 : 02
% 12 : H2
% 13 : O3-
% 14 : O3
% 15 : HO3
% 16 : H2O

clear all
close all
clc

% Initial concentrations (micromolar)%
% Order: [electrons H+ OH- H2O2 HO2- H OH O- HO2 O2- O2 H2 O3- O3 HO3 H2O];
initialConcentration = [0 10^-7 10^-7 0 0 0 0 0 0 0 0 0 0 0 0 55.56]*10^6;

% number of species%
numberOfSpecies = length(initialConcentration) ;

%Set of Dose Rate in Gy/s (A)%
doseRate = 0.25;%10^(10.);

% G-Values (molecules/100 eV)%
% Deafault G-Values taken from Hill and Smith at 300KeV [ref] %
gValues = zeros(numberOfSpecies,1);
% gValues(1)  = 3.47 / 100;  % Hydrated Electrons
% gValues(2)  = 4.42 / 100;  % H+
% gValues(3)  = 0.95 / 100;  % OH-
% gValues(4)  = 0.47 / 100;  % H2O2
% gValues(6)  = 1.00 / 100;  % H
% gValues(7)  = 3.63 / 100;  % OH
% gValues(9)  = 0.08 / 100;  % HO2
% gValues(12) = 0.17 / 100;  % H2
% gValues(16) = - ( gValues(3) + 2*gValues(4) + gValues(7)+ 2*gValues(9) );
% Note the G-Value for Water can be taken from either an atomic balance on
% hydrogen or oxygen. Above (gValues(16)) is for oxygen.


%If you want to test the Pastina model again, use these Gvalues and set
%current to .25/(7.5*10^16) to account for the conversion to dosage
gValues(1) = 2.6 / 100;%
gValues(2) = 3.1 / 100;%
gValues(3) = .5 / 100;%
gValues(4) = .7 / 100;%
gValues(6) = .66 / 100;%
gValues(7) = 2.7 / 100;%
gValues(9) = .02 / 100;%
gValues(12) = .45 / 100;%
gValues(16) = -gValues(7)-2*gValues(4)-2*gValues(9) - gValues(3);


%time span (s)%
timeSpan = [0 10^3];


     
[times,concentrations] = ode23s(@water_radiolysis, timeSpan, initialConcentration,...
    [], doseRate , gValues);




