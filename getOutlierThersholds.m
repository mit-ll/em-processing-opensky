% Copyright 2018 - 2020, MIT Lincoln Laboratory
% SPDX-License-Identifier: BSD-2-Clause
function [maxAlt_ft_msl,outlierSpeed_kt,outlierTurnRate_deg_s,outlierVertRate_ft_s] = getOutlierThersholds(acType)

switch acType
    case 'FixedWingMultiEngine'
        maxAlt_ft_msl = 85000;
        outlierVertRate_ft_s = 6000 / 60;
        outlierSpeed_kt = 600;  % assume
        outlierTurnRate_deg_s = 10; % assume
    case 'FixedWingSingleEngine'
        % Cessna 206 Series (27k ft): https://www.easa.europa.eu/sites/default/files/dfu/EASA.IM_.A.053_Issue%207_20180621.pdf
        % Cessna 208: https://rgl.faa.gov/Regulatory_and_Guidance_Library/rgMakeModel.nsf/0/084331c4c5578a27862581ca00733289/$FILE/A37CE_Rev_22.pdf
        % Cessna 210R: https://rgl.faa.gov/Regulatory_and_Guidance_Library/rgMakeModel.nsf/0/47f698e4f039604186257ed2006ca246/$FILE/3A21_Rev_49.pdf
        % Extra 300L: https://rgl.faa.gov/Regulatory_and_Guidance_Library/rgMakeModel.nsf/0/0540a38cd12850bc8625805d0058ff51/$FILE/A67EU_Rev_12.pdf
        % Cessna Denali and Pilatus PC-12 advertise 30K
        maxAlt_ft_msl = 30000 + 1000; % 25K was normally listed (Extra was 16K), so 30K is conservative
        outlierVertRate_ft_s = 6000 / 60;
        outlierSpeed_kt = 400; % P-51 mustang top speed was 380 knots
        outlierTurnRate_deg_s = 10; % assume
    case 'Glider'
        % DG-800: https://rgl.faa.gov/Regulatory_and_Guidance_Library/rgMakeModel.nsf/0/19e41de68d031e6e862580f00053dc74/$FILE/G01CE_Rev_6.pdf
        % DG Flugzeugbau DG-800: https://www.easa.europa.eu/sites/default/files/dfu/EASA-TCDS-A.067_DG_Flugzeugbau_DG--800-04-02112010.pdf
        maxAlt_ft_msl = 18000; % FAA TCDS don't state, we assume
        outlierVertRate_ft_s = 2000 / 60; % assume
        outlierSpeed_kt = 200; % assume
        outlierTurnRate_deg_s = 10; % assume
    case 'Gyroplane'
        maxAlt_ft_msl = 18000; % assume
        outlierVertRate_ft_s = 2000 / 60; % assume
        outlierSpeed_kt = 130; % assume
        outlierTurnRate_deg_s = 12; % assume
    case 'PoweredParachute'
        maxAlt_ft_msl = 18000; % assume
        outlierVertRate_ft_s = 2000 / 60; % assume
        outlierSpeed_kt = 100; % assume
        outlierTurnRate_deg_s = 12; % assume
    case 'Rotorcraft'
        % Bell 412 EP: https://rgl.faa.gov/Regulatory_and_Guidance_Library/rgMakeModel.nsf/0/df904d3232a709818625851f00686d16/$FILE/H4SW_Rev36.pdf
        % Airbus AS350 (16k ft): https://www.easa.europa.eu/sites/default/files/dfu/tcds_easa_r008_ah_as350_ec130_issue_14.pdf
        % Bell 212/412 Series (20k ft): https://www.easa.europa.eu/sites/default/files/dfu/TCDS_EASA_IM_R106_Bell_212_412_Issue_03.pdf
        % Bell 222/230/430 (20k ft): https://www.easa.europa.eu/sites/default/files/dfu/TCDS_EASA_IM_R114_Bell222_230_430_Issue_01.pdf
        % Robinson R44 (14k ft): https://www.easa.europa.eu/sites/default/files/dfu/TCDS_EASA_IM_R121_RHC_R44_Issue_06.pdf
        maxAlt_ft_msl = 20000 + 1000; % FAA TCDS don't state, we assume
        outlierVertRate_ft_s = 2000 / 60; % assume
        outlierSpeed_kt = 250; % Eurocopter X3 unofficial speed record is 255
        outlierTurnRate_deg_s = 12; % assume
    case 'WeightShiftControl'
        maxAlt_ft_msl = 18000; % assume
        outlierVertRate_ft_s = 2000 / 60; % assume
        outlierSpeed_kt = 100; % assume
        outlierTurnRate_deg_s = 12; % assume
    otherwise
        maxAlt_ft_msl = 85000;
        outlierVertRate_ft_s = 6000 / 60;
        outlierSpeed_kt = 600;  % assume
        outlierTurnRate_deg_s = 12; % assume
end