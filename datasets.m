%List of datasets where the saturating stimulus lasts only 30s, not 5 min
function dataparms = datasets()
    %% No Background, MeAsp Dose Response
    dataparms(1).Removal = 0;
    dataparms(1).backConc = 0; %Background concentration for fitting hill function
    dataparms(1).xlabels = '[MeAsp]'; %xlabel on plots
    dataparms(1).Lplot = 10.^[-3:0.01:2]; %Range of plots
%     dataparms(1).Files = ["./Data/230417_FOV1.mat"];
    dataparms(1).Files = ["./Data/220106_FOV1.mat", "./Data/230417_FOV1.mat" ];
%     dataparms(1).Files = ["./Data/220106_FOV1.mat", "./Data/230417_FOV1.mat", "./Data/210802_FOV1.mat", "./Data/210805_FOV1.mat"];
    dataparms(1).OutFolder = 'MeAspDR_OneConcentrationSets/'; %Where to output data to   
    dataparms(1).concLevels = [0.2, 0.5, 1, 2, 4;...
                                0.2, 0.5, 1, 2, 4; ...23
                                0.2, 0.5, 1, 2, 4; ...
                                0.2, 0.5, 1, 2, 4];
    dataparms(1).Hillp0 = [log(1.2), log(1.3), 1]; %Log of estimated parameters to fit hill function [log(n), log(K)]
    dataparms(1).Normp0 = [log(1.3), 1, 1]; %[log(mean), sigma]
    
    %% 1uM Ser Background, MeAsp Dose Response
    dataparms(2).Removal = 0;
    dataparms(2).backConc = 0; %Background concentration for fitting hill function
    dataparms(2).xlabels = '[MeAsp]'; %xlabel on plots
    dataparms(2).Lplot = 10.^[-3:0.01:2]; %Range of plots
    dataparms(2).Files = ["./Data/210810_FOV1.mat", "./Data/210811_FOV1.mat"];
%     dataparms(2).Files = ["./Data/210811_FOV1.mat"];
    dataparms(2).OutFolder = '1SerBackMeAspDR/'; %Where to output data to   
    dataparms(2).concLevels = [0.2, 0.5, 1, 2, 4; 0.2, 0.5, 1, 2, 4];
    dataparms(2).Hillp0 = [log(1.2), log(3), 1]; %Log of estimated parameters to fit hill function
    dataparms(2).Normp0 = [log(1.2), 0.6, 1];
    
    %% 100uM NaGlu Background, MeAsp Dose Response
    dataparms(3).Removal = 0;
    dataparms(3).backConc = 0; %Background concentration for fitting hill function
    dataparms(3).xlabels = '[MeAsp]'; %xlabel on plots
    dataparms(3).Lplot = 10.^[-3:0.01:2]; %Range of plots
    dataparms(3).Files = ["./Data/210812_FOV1.mat", "./Data/210813_FOV1.mat", "./Data/230705_FOV1.mat"];
%     dataparms(3).Files = ["./Data/230705_FOV1.mat"];
    dataparms(3).OutFolder = '100GluBackMeAspDR/'; %Where to output data to   
    dataparms(3).concLevels = [0.2, 0.5, 1, 2, 4; 0.2, 0.5, 1, 2, 4; 0.2, 0.5, 1, 2, 4];
    dataparms(3).Hillp0 = [log(1.2), log(3), 1]; %Log of estimated parameters to fit hill function
    dataparms(3).Normp0 = [-.1, .8, 1];
    
    %% 100uM MeAsp Background, MeAsp Dose Response
    dataparms(4).Removal = 0;
    dataparms(4).backConc = 100; %Background concentration for fitting hill function
    dataparms(4).xlabels = '[MeAsp]'; %xlabel on plots
    dataparms(4).Lplot = 10.^[1:0.01:3]; %Range of plots
    dataparms(4).Files = ["./Data/210816_FOV1.mat", "./Data/230717_FOV1.mat", "./Data/230718_FOV1.mat"];
%     dataparms(4).Files = ["./Data/230717_FOV1.mat", "./Data/230718_FOV1.mat"];
    dataparms(4).OutFolder = '100MeAspBackMeAspDR/'; %Where to output data to   
    dataparms(4).concLevels = [102, 105, 110, 120, 140; 102, 105, 110, 120, 140; 102, 105, 110, 120, 140];
    dataparms(4).Hillp0 = [log(2), log(10^2), 1]; %Log of estimated parameters to fit hill function
    dataparms(4).Normp0 = [log(110), log(2), 1];
    
    %% 10uM L-Asp Background, MeAsp Dose Response
    dataparms(5).Removal = 0;
    dataparms(5).backConc = 0; %Background concentration for fitting hill function
    dataparms(5).xlabels = '[MeAsp]'; %xlabel on plots
    dataparms(5).Lplot = 10.^(-2:0.01:3); %Range of plots
    dataparms(5).Files = ["./Data/210818_FOV1.mat", "./Data/210908_FOV1.mat"]; %Where to load data from
    dataparms(5).OutFolder = '10LAspBackMeAspDose/'; %Where to output data to
    dataparms(5).concLevels = [5, 10, 20, 40, 80; 5, 10, 20, 40, 80];
    dataparms(5).Hillp0 = [log(1.2), log(30), 1]; %Log of estimated parameters to fit hill function
    dataparms(5).Normp0 = [2, 30, 1];
    
    %% 0 Background, Ser Dose Response
    dataparms(6).Removal = 0;
    dataparms(6).backConc = 0; %Background concentration for fitting hill function
    dataparms(6).xlabels = '[Ser]'; %xlabel on plots
    dataparms(6).Lplot = 10.^[-4:0.01:2]; %Range of plots
%     dataparms(6).Files = ["./Data/210910_FOV1.mat", "./Data/211006_FOV1.mat", "./Data/210913_FOV1.mat", "./Data/220111_FOV1.mat"]; %Where to load data from
%     dataparms(6).OutFolder = '0BackSerDose/'; %Where to output data to
%     dataparms(6).concLevels = [.005, 0.01, 0.02, 0.05, 0.1; .005, 0.01, 0.02, 0.05, 0.1; .005, 0.01, 0.02, 0.05, 0.1; .005, 0.01, 0.02, 0.05, 0.1];
    dataparms(6).Files = ["./Data/210910_FOV1.mat",  "./Data/210913_FOV1.mat", "./Data/220111_FOV1.mat"]; %Where to load data from
    dataparms(6).OutFolder = '0BackSerDose/'; %Where to output data to
    dataparms(6).concLevels = [.005, 0.01, 0.02, 0.05, 0.1; .005, 0.01, 0.02, 0.05, 0.1; .005, 0.01, 0.02, 0.05, 0.1];
    dataparms(6).Hillp0 = [log(.05), log(3), 1]; %Log of estimated parameters to fit hill function
    dataparms(6).Normp0 = [log(0.03), 1, 1];

    %% 1uM Ser Background, Ser Dose Response
    dataparms(7).Removal = 0;
    dataparms(7).backConc = 1; %Background concentration for fitting hill function
    dataparms(7).xlabels = '[Ser]'; %xlabel on plots
    dataparms(7).Lplot = 10.^[-4:0.01:2]; %Range of plots
    dataparms(7).Files = ["./Data/210914_FOV1.mat", "./Data/210916_FOV1.mat"]; %Where to load data from
    dataparms(7).OutFolder = '1SerBackSerDose/'; %Where to output data to
    dataparms(7).concLevels = [1.02, 1.05, 1.1, 1.2, 1.4; 1.02, 1.05, 1.1, 1.2, 1.4];
    dataparms(7).Hillp0 = [log(1.07), log(1.12), 1]; %Log of estimated parameters to fit hill function
    dataparms(7).Normp0 = [1.07, 1, 1];
    
    %% 100uM MeAsp Background, Ser Dose Response
    dataparms(8).Removal = 0;
    dataparms(8).backConc = 0; %Background concentration for fitting hill function
    dataparms(8).xlabels = '[Ser]'; %xlabel on plots
    dataparms(8).Lplot = 10.^[-4:0.01:2]; %Range of plots
%     dataparms(8).Files = ["./Data/210920_FOV1.mat", "./Data/210921_FOV1.mat"]; %Where to load data from
%     dataparms(8).OutFolder = '100MeAspBackSerDose/'; %Where to output data to
%     dataparms(8).concLevels = [.005, 0.01, 0.02, 0.05, 0.1; .005, 0.01, 0.02, 0.05, 0.1];
    dataparms(8).Files = ["./Data/210921_FOV1.mat"]; %Where to load data from
    dataparms(8).OutFolder = '100MeAspBackSerDose/'; %Where to output data to
    dataparms(8).concLevels = [.005, 0.01, 0.02, 0.05, 0.1];
    dataparms(8).Hillp0 = [log(1.2), log(.1), 1]; %Log of estimated parameters to fit hill function
    dataparms(8).Normp0 = [log(0.1), 1, 1];
    
    %% 0 background, Glu Dose Response
    dataparms(9).Removal = 0;
    dataparms(9).backConc = 0; %Background concentration for fitting hill function
    dataparms(9).xlabels = '[Glu]'; %xlabel on plots
    dataparms(9).Lplot = 10.^[-3:0.01:4]; %Range of plots
    dataparms(9).Files = ["./Data/211012_FOV1.mat", "./Data/211011_FOV1.mat"]; %Where to load data from
    dataparms(9).OutFolder = '0BackGluDose/'; %Where to output data to
    dataparms(9).concLevels = [5, 10, 20, 50, 100; 5, 10, 20, 50, 100];
    dataparms(9).Hillp0 = [log(10), log(3), 1]; %Log of estimated parameters to fit hill function
    dataparms(9).Normp0 = [1, 0.5, 1];
    
    %% 100uM MeAsp background, Glu Dose Response
    dataparms(10).Removal = 0;
    dataparms(10).backConc = 0; %Background concentration for fitting hill function
    dataparms(10).xlabels = '[Glu]'; %xlabel on plots
    dataparms(10).Lplot = 10.^[-3:0.01:4]; %Range of plots
    dataparms(10).Files = ["./Data/211013_FOV1.mat"]; %Where to load data from
    dataparms(10).OutFolder = '100MeAspBackGluDose/'; %Where to output data to
    dataparms(10).concLevels = [5, 10, 20, 50, 100];
    dataparms(10).Hillp0 = [log(10), log(3), 1]; %Log of estimated parameters to fit hill function
    dataparms(10).Normp0 = [1, 0.5, 1];
    
    %% 1mM Glu background, MeAsp Dose Response
    dataparms(11).Removal = 0;
    dataparms(11).backConc = 0; %Background concentration for fitting hill function
    dataparms(11).xlabels = '[MeAsp]'; %xlabel on plots
    dataparms(11).Lplot = 10.^[-3:0.01:4]; %Range of plots
    dataparms(11).Files = ["./Data/211020_FOV1.mat", "./Data/211021_FOV1.mat"]; %Where to load data from
%     dataparms(11).Files = ["./Data/211020_FOV1.mat"]; %Where to load data from    
    dataparms(11).OutFolder = '1mMGluBackMeAspDose/'; %Where to output data to
    dataparms(11).concLevels = [0.5, 1, 2, 4, 8; 0.5, 1, 2, 4, 8];
%     dataparms(11).concLevels = [0.5, 1, 2, 4, 8];
    dataparms(11).Hillp0 = [log(1.2), log(3), 1]; %Log of estimated parameters to fit hill function
    dataparms(11).Normp0 = [log(1.1), 1];
    
    %% 0 Background, L-Asp Dose Response
    dataparms(12).Removal = 0;
    dataparms(12).backConc = 0; %Background concentration for fitting hill function
    dataparms(12).xlabels = '[L-Asp]'; %xlabel on plots
    dataparms(12).Lplot = 10.^[-4:0.01:3]; %Range of plots
    dataparms(12).Files = ["./Data/211022_FOV1.mat", "./Data/211025_FOV1.mat"]; %Where to load data from
    dataparms(12).OutFolder = '0BackLAspDose/'; %Where to output data to
    dataparms(12).concLevels = [0.02, 0.05, 0.1, 0.2, 0.4; 0.02, 0.05, 0.1, 0.2, 0.4];
    dataparms(12).Hillp0 = [log(0.04), log(3), 1]; %Log of estimated parameters to fit hill function
    dataparms(12).Normp0 = [log(0.02), 1, 1];
    
    %% 100 uM MeAsp background, L-Asp Dose Response
    dataparms(13).Removal = 0;
    dataparms(13).backConc = 0; %Background concentration for fitting hill function
    dataparms(13).xlabels = '[L-Asp]'; %xlabel on plots
    dataparms(13).Lplot = 10.^[-4:0.01:3]; %Range of plots
    dataparms(13).Files = ["./Data/211026_FOV1.mat", "./Data/211030_FOV2.mat"]; %Where to load data from
    dataparms(13).OutFolder = '100MeAspBackLAspDose/'; %Where to output data to
    dataparms(13).concLevels = [0.05, 0.1, 0.2, 0.4, 0.8; 0.05, 0.1, 0.2, 0.4, 0.8];
    dataparms(13).Hillp0 = [log(2), log(.01), 1]; %Log of estimated parameters to fit hill function
    dataparms(13).Normp0 = [log(.15), 2, 1];
    
    %% 100 uM MeAsp background, Glu Dose Response
    dataparms(14).Removal = 0;
    dataparms(14).backConc = 0; %Background concentration for fitting hill function
    dataparms(14).xlabels = '[Glu]'; %xlabel on plots
    dataparms(14).Lplot = 10.^[0:0.01:5]; %Range of plots
%     dataparms(14).Files = ["./Data/211027_FOV1.mat"]; %Where to load data from
%     dataparms(14).Files = ["./Data/211028_FOV1.mat", "./Data/211029_FOV1.mat", "./Data/211102_FOV1.mat"]; %Where to load data from
    dataparms(14).Files = ["./Data/211102_FOV1.mat", "./Data/211103_FOV1.mat"]; %Where to load data from  
    dataparms(14).OutFolder = '100MeAspBackGluDose/'; %Where to output data to
    dataparms(14).concLevels = [30, 75, 150, 300, 600; 30, 75, 150, 300, 600];
%     dataparms(14).concLevels = [50, 100, 200, 400, 800; 50, 100, 200, 400, 800; 30, 75, 150, 300, 600];
    dataparms(14).Hillp0 = [log(2), log(150), 1]; %Log of estimated parameters to fit hill function
    dataparms(14).Normp0 = [log(150), 5, 1];
    
    %% 100 uM MeAsp Background, analyzed as 0 background
    dataparms(15).Removal = 0;
    dataparms(15).backConc = 0; %Background concentration for fitting hill function
    dataparms(15).xlabels = '[MeAsp]'; %xlabel on plots
    dataparms(15).Lplot = 10.^[-3:0.01:3]; %Range of plots
    dataparms(15).Files = ["./Data/210816_FOV1.mat"];
    dataparms(15).OutFolder = '100MeAspBackMeAspDR_AsIf0/'; %Where to output data to   
    dataparms(15).concLevels = [2, 5, 10, 20, 40];
    dataparms(15).Hillp0 = [log(2), log(10), 1]; %Log of estimated parameters to fit hill function
    dataparms(15).Normp0 = [log(10), log(2), 1];
    
    %% 10 uM MeAsp Background, MeAsp dose response
    dataparms(16).Removal = 0;
    dataparms(16).backConc = 10; %Background concentration for fitting hill function
    dataparms(16).xlabels = '[MeAsp]'; %xlabel on plots
    dataparms(16).Lplot = 10.^[-3:0.01:3]; %Range of plots
    dataparms(16).Files = ["./Data/220302_FOV1.mat", "./Data/220303_FOV1.mat"];
    dataparms(16).OutFolder = '10MeAspBackMeAspDR/'; %Where to output data to   
    dataparms(16).concLevels = [10.5, 11, 12, 14, 18; 10.5, 11, 12, 14, 18];
    dataparms(16).Hillp0 = [log(2), log(13), 1]; %Log of estimated parameters to fit hill function
    dataparms(16).Normp0 = [log(12), 1, 1];
    
    %% 10 uM MeAsp + 10uM L-Asp Background, MeAsp dose response
    dataparms(17).Removal = 0;
    dataparms(17).backConc = 10; %Background concentration for fitting hill function
    dataparms(17).xlabels = '[MeAsp]'; %xlabel on plots
    dataparms(17).Lplot = 10.^[-3:0.01:3]; %Range of plots
    dataparms(17).Files = ["./Data/220324_FOV1.mat", "./Data/230511_FOV1.mat"];
    dataparms(17).OutFolder = '10MeAsp10LAspBackMeAspDR/'; %Where to output data to   
    dataparms(17).concLevels = [15, 20, 30, 50, 90; 15, 20, 30, 50, 90];
    dataparms(17).Hillp0 = [log(2), log(25), 1]; %Log of estimated parameters to fit hill function
    dataparms(17).Normp0 = [log(25), log(2), 1];

    %% 100uM MeAsp + 1uM Ser background, MeAsp dose response
    dataparms(18).Removal = 0;
    dataparms(18).backConc = 100; %Background concentration for fitting hill function
    dataparms(18).xlabels = '[MeAsp]'; %xlabel on plots
    dataparms(18).Lplot = 10.^[1:0.01:3]; %Range of plots
    dataparms(18).Files = ["./Data/220426_FOV1.mat"];
    dataparms(18).OutFolder = '100MeAsp1SerBackMeAspDR/'; %Where to output data to   
    dataparms(18).concLevels = [102, 105, 110, 120, 140];
    dataparms(18).Hillp0 = [log(2), log(10^2), 1]; %Log of estimated parameters to fit hill function
    dataparms(18).Normp0 = [log(110), log(2), 1];

    %% 100uM MeAsp + 1uM Ser background, Ser Dose Response
    dataparms(19).Removal = 0;
    dataparms(19).backConc = 1; %Background concentration for fitting hill function
    dataparms(19).xlabels = '[Ser]'; %xlabel on plots
    dataparms(19).Lplot = 10.^[-4:0.01:2]; %Range of plots
    dataparms(19).Files = ["./Data/220428_FOV1.mat"]; %Where to load data from
    dataparms(19).OutFolder = '100MeAsp1SerBackSerDR/'; %Where to output data to
    dataparms(19).concLevels = [1.02, 1.05, 1.1, 1.2, 1.4];
    dataparms(19).Hillp0 = [log(1.07), log(3), 1]; %Log of estimated parameters to fit hill function
    dataparms(19).Normp0 = [log(1.07), log(2), 1];

    %% 0.3uM MeAsp background, meAsp dose response
    dataparms(20).Removal = 0;
    dataparms(20).backConc = 0.3; %Background concentration for fitting hill function
    dataparms(20).xlabels = '[MeAsp]'; %xlabel on plots
    dataparms(20).Lplot = 10.^[-3:0.01:2]; %Range of plots
    dataparms(20).Files = ["./Data/220615_FOV1.mat", "./Data/230410_FOV1.mat"]; %Where to load data from
%     dataparms(20).Files = ["./Data/220615_FOV1.mat"]; %Where to load data from
    dataparms(20).OutFolder = '0_3MeAspBackMeAspDR/'; %Where to output data to
    dataparms(20).concLevels = [0.5, 0.8, 1.3, 2.3, 4.3; 0.5, 0.8, 1.3, 2.3, 4.3];
    dataparms(20).Hillp0 = [log(1.07), log(1), 1]; %Log of estimated parameters to fit hill function
    dataparms(20).Normp0 = [log(1.07), log(2), 1];

    %% 10uM L-asp, 1 uM MeAsp Background, meAsp dose response
    dataparms(21).Removal = 0;
    dataparms(21).backConc = 1;
    dataparms(21).xlabels = '[MeAsp]';
    dataparms(21).Lplot = 10.^[-3:0.01:2]; %Range of plots    
    dataparms(21).Files = ["./Data/230330_FOV1.mat", "./Data/230403_FOV1.mat"]; %Where to load data from
    dataparms(21).OutFolder = '10LAsp_1meAspBackMeAspDR/'; %Where to output data to
    dataparms(21).concLevels = [5, 10, 20, 40, 80; 5, 10, 20, 40, 80];
    dataparms(21).Hillp0 = [log(1.07), log(30), 1]; %Log of estimated parameters to fit hill function
    dataparms(21).Normp0 = [log(15), log(2), 1];

    %% 0.1 uM MeAsp Background, meAsp dose response
    dataparms(22).Removal = 0;
    dataparms(22).backConc = 0.1;
    dataparms(22).xlabels = '[MeAsp]';
    dataparms(22).Lplot = 10.^[-3:0.01:2]; %Range of plots    
%     dataparms(22).Files = ["./Data/230404_FOV1.mat"]; %Where to load data from
%     dataparms(22).Files = ["./Data/230404_FOV1.mat", "./Data/230413_FOV1.mat", "./Data/230414_FOV1.mat"]; %Where to load data from
    dataparms(22).Files = ["./Data/230830_FOV1.mat", "./Data/230830_FOV2.mat", "./Data/230831_FOV1.mat", "./Data/230831_FOV2.mat"]; %Where to load data from
    dataparms(22).OutFolder = '0_1meAspBackMeAspDR/'; %Where to output data to
    dataparms(22).concLevels = [0.3, 0.6, 1.1, 2.1, 4.1; 0.3, 0.6, 1.1, 2.1, 4.1; 0.3, 0.6, 1.1, 2.1, 4.1; 0.3, 0.6, 1.1, 2.1, 4.1];
    dataparms(22).Hillp0 = [log(1.07), log(1), 1]; %Log of estimated parameters to fit hill function
    dataparms(22).Normp0 = [log(1), log(2), 1];

    %% 1uM MeAsp Background, meAsp dose
    dataparms(23).Removal = 0;
    dataparms(23).backConc = 1;
    dataparms(23).xlabels = '[MeAsp]';
    dataparms(23).Lplot = 10.^[-3:0.01:2]; %Range of plots    
%     dataparms(23).Files = ["./Data/230425_FOV1.mat", "./Data/230426_FOV1.mat"]; %Where to load data from
    dataparms(23).Files = ["./Data/230428_FOV1.mat", "./Data/230429_FOV1.mat"]; %Where to load data from
    dataparms(23).OutFolder = '1meAspBackMeAspDR/'; %Where to output data to
    dataparms(23).concLevels = [1.2, 1.5, 2, 3, 5; 1.2, 1.5, 2, 3, 5];
    dataparms(23).Hillp0 = [log(1.07), log(1), 1]; %Log of estimated parameters to fit hill function
    dataparms(23).Normp0 = [log(1), log(2), 1];

    %% 2uM L-Asp + 1uM MeAsp Background, meAsp dose
    dataparms(24).Removal = 0;
    dataparms(24).backConc = 1;
    dataparms(24).xlabels = '[MeAsp]';
    dataparms(24).Lplot = 10.^[-3:0.01:2]; %Range of plots    
%     dataparms(24).Files = ["./Data/230504_FOV1.mat", "./Data/230504_FOV2.mat","./Data/230505_FOV1.mat", "./Data/230505_FOV2.mat" ]; %Where to load data from
    dataparms(24).Files = ["./Data/230504_FOV1.mat", "./Data/230505_FOV1.mat"]; %Where to load data from
    dataparms(24).OutFolder = '2LAsp_1MeAspBack_MeAspDR/'; %Where to output data to
    dataparms(24).concLevels = [2, 4, 8, 16, 32; 2, 4, 8, 16, 32; 2, 4, 8, 16, 32; 2, 4, 8, 16, 32];
    dataparms(24).Hillp0 = [log(2.1), log(1), 1]; %Log of estimated parameters to fit hill function
    dataparms(24).Normp0 = [log(1), log(2.1), 1];

    %% 2uM L-Asp + 10uM MeAsp Background, meAsp dose
    dataparms(25).Removal = 0;
    dataparms(25).backConc = 10;
    dataparms(25).xlabels = '[MeAsp]';
    dataparms(25).Lplot = 10.^[-3:0.01:2]; %Range of plots    
    dataparms(25).Files = ["./Data/230509_FOV1.mat", "./Data/230510_FOV1.mat","./Data/230510_FOV2.mat"]; %Where to load data from
%     dataparms(25).Files = ["./Data/230509_FOV1.mat", "./Data/230510_FOV1.mat"]; %Where to load data from
    dataparms(25).OutFolder = '2LAsp_10MeAspBack_MeAspDR/'; %Where to output data to
    dataparms(25).concLevels = [12,14, 18, 26, 42; 12,14, 18, 26, 42; 12,14, 18, 26, 42];
    dataparms(25).Hillp0 = [log(2.1), log(12), 1]; %Log of estimated parameters to fit hill function
    dataparms(25).Normp0 = [log(15), 2, 1];
    
    %% 1uM L-Asp, L-Asp dose
    dataparms(26).Removal = 0;
    dataparms(26).backConc = 1;
    dataparms(26).xlabels = '[L-Asp]';
    dataparms(26).Lplot = 10.^[-3:0.01:2]; %Range of plots    
%     dataparms(26).Files = ["./Data/230523_FOV1.mat", "./Data/230523_FOV2.mat","./Data/230524_FOV1.mat", "./Data/230524_FOV2.mat"]; %Where to load data from
    dataparms(26).Files = ["./Data/230523_FOV1.mat","./Data/230524_FOV1.mat"]; %Where to load data from
    dataparms(26).OutFolder = '1LAsp_LAspDR/'; %Where to output data to
    dataparms(26).concLevels = [1.05, 1.1, 1.2, 1.4, 1.8; 1.05, 1.1, 1.2, 1.4, 1.8; 1.05, 1.1, 1.2, 1.4, 1.8; 1.05, 1.1, 1.2, 1.4, 1.8];
    dataparms(26).Hillp0 = [log(2.1), log(1.1), 1]; %Log of estimated parameters to fit hill function
    dataparms(26).Normp0 = [log(1), log(1.1), 1];

    %% 0.5 uM L-Asp, L-Asp Dose
    dataparms(27).Removal = 0;
    dataparms(27).backConc = 0.5;
    dataparms(27).xlabels = '[L-Asp]';
    dataparms(27).Lplot = 10.^[-3:0.01:2]; %Range of plots    
    dataparms(27).Files = ["./Data/230526_FOV1.mat", "./Data/230526_FOV2.mat"]; %Where to load data from
    dataparms(27).OutFolder = '0_5LAsp_LAspDR/'; %Where to output data to
    dataparms(27).concLevels = [0.52, 0.55, 0.6, 0.7, 0.9; 0.52, 0.55, 0.6, 0.7, 0.9];
    dataparms(27).Hillp0 = [log(2.1), log(0.6), 1]; %Log of estimated parameters to fit hill function
    dataparms(27).Normp0 = [0.5, 0.5, 1];

    %% 0.1 uM L-Asp, L-Asp Dose
    dataparms(28).Removal = 0;
    dataparms(28).backConc = 0.1;
    dataparms(28).xlabels = '[L-Asp]';
    dataparms(28).Lplot = 10.^[-3:0.01:2]; %Range of plots    
    dataparms(28).Files = ["./Data/230601_FOV1.mat", "./Data/230602_FOV1.mat"]; %Where to load data from
    dataparms(28).OutFolder = '0_1LAsp_LAspDR/'; %Where to output data to
    dataparms(28).concLevels = [0.12, 0.15, 0.2, 0.3, 0.5; 0.12, 0.15, 0.2, 0.3, 0.5;];
    dataparms(28).Hillp0 = [log(.9), log(0.2), 1]; %Log of estimated parameters to fit hill function
    dataparms(28).Normp0 = [-.5, .8, 1];


    %% 0.05 uM L-Asp, L-Asp Dose
    dataparms(29).Removal = 0;
    dataparms(29).backConc = 0.05;
    dataparms(29).xlabels = '[L-Asp]';
    dataparms(29).Lplot = 10.^[-3:0.01:2]; %Range of plots    
    dataparms(29).Files = ["./Data/230605_FOV1.mat", "./Data/230607_FOV1.mat", "./Data/230608_FOV1.mat"]; %Where to load data from
    dataparms(29).OutFolder = '0_05LAsp_LAspDR/'; %Where to output data to
    dataparms(29).concLevels = [.07, .1, .15, .25, .45; .07, .1, .15, .25, .45; .07, .1, .15, .25, .45];
    dataparms(29).Hillp0 = [log(2.1), log(0.6), 1]; %Log of estimated parameters to fit hill function
    dataparms(29).Normp0 = [0.05, 1, 1];

    %% 1 uM L-Asp + 100uM MeAsp, L-Asp Dose
    dataparms(30).Removal = 0;
    dataparms(30).backConc = 1;
    dataparms(30).xlabels = '[L-Asp]';
    dataparms(30).Lplot = 10.^[-3:0.01:2]; %Range of plots    
%     dataparms(30).Files = ["./Data/230615_FOV1.mat", "./Data/230615_FOV2.mat", "./Data/230620_FOV1.mat"]; %Where to load data from
    dataparms(30).Files = ["./Data/230615_FOV1.mat",  "./Data/230620_FOV1.mat"];
    dataparms(30).OutFolder = '1LAsp_100meAsp_LAspDR/'; %Where to output data to
    dataparms(30).concLevels = [1.1, 1.2, 1.4, 1.8, 2.6; 1.1, 1.2, 1.4, 1.8, 2.6; 1.1, 1.2, 1.4, 1.8, 2.6];
    dataparms(30).Hillp0 = [log(2.1), log(0.6), 1]; %Log of estimated parameters to fit hill function
    dataparms(30).Normp0 = [0.05, 1, 1];

    %% 10 uM L-Asp + 100uM MeAsp, L-Asp Dose
    dataparms(31).Removal = 0;
    dataparms(31).backConc = 10;
    dataparms(31).xlabels = '[L-Asp]';
    dataparms(31).Lplot = 10.^[-3:0.01:2]; %Range of plots    
    dataparms(31).Files = ["./Data/230621_FOV1.mat", "./Data/230623_FOV1.mat"]; %Where to load data from
    dataparms(31).OutFolder = '10LAsp_100meAsp_LAspDR/'; %Where to output data to
    dataparms(31).concLevels = [10.5, 11, 12, 14, 18; 10.5, 11, 12, 14, 18];
    dataparms(31).Hillp0 = [log(2.1), log(0.6), 1]; %Log of estimated parameters to fit hill function
    dataparms(31).Normp0 = [0.05, 1, 1];

    %% 10 uM L-Asp, L-Asp Dose
    dataparms(32).Removal = 0;
    dataparms(32).backConc = 10;
    dataparms(32).xlabels = '[L-Asp]';
    dataparms(32).Lplot = 10.^[-3:0.01:2]; %Range of plots    
    dataparms(32).Files = ["./Data/230627_FOV1.mat", "./Data/230627_FOV2.mat", "./Data/230628_FOV1.mat"]; %Where to load data from
    dataparms(32).OutFolder = '10LAsp_LAspDR/'; %Where to output data to
    dataparms(32).concLevels = [10.2, 10.5, 11, 12, 14; 10.2, 10.5, 11, 12, 14; 10.2, 10.5, 11, 12, 14];
    dataparms(32).Hillp0 = [log(2.1), log(0.6), 1]; %Log of estimated parameters to fit hill function
    dataparms(32).Normp0 = [0.05, 1, 1];

    %% 1uM L_Asp, meAsp Dose
    dataparms(33).Removal = 0;
    dataparms(33).backConc = 0;
    dataparms(33).xlabels = '[meAsp]';
    dataparms(33).Lplot = 10.^[-3:0.01:2]; %Range of plots    
    dataparms(33).Files = ["./Data/230629_FOV1.mat", "./Data/230704_FOV1.mat"]; %Where to load data from
    dataparms(33).OutFolder = '1LAsp_meAspDR/'; %Where to output data to
    dataparms(33).concLevels = [0.5, 1, 2, 4, 8; 0.5, 1, 2, 4, 8];
    dataparms(33).Hillp0 = [log(2.1), log(0.6), 1]; %Log of estimated parameters to fit hill function
    dataparms(33).Normp0 = [1, 1, 1];

    %% 0.01uM L_Asp, LAsp Dose
    dataparms(34).Removal = 0;
%     dataparms(34).backConc = 0.01;
    dataparms(34).backConc = 0;
    dataparms(34).xlabels = '[L-Asp]';
    dataparms(34).Lplot = 10.^[-3:0.01:2]; %Range of plots    
    dataparms(34).Files = ["./Data/230706_FOV1.mat", "./Data/230707_FOV1.mat"]; %Where to load data from
    dataparms(34).OutFolder = '0_01LAsp_LAspDR/'; %Where to output data to
%     dataparms(34).concLevels = [0.03, 0.06, 0.11, 0.21, 0.41; 0.03, 0.06, 0.11, 0.21, 0.41];
    dataparms(34).concLevels = [0.02, 0.05, 0.1, 0.2, 0.4; 0.02, 0.05, 0.1, 0.2, 0.4];
    dataparms(34).Hillp0 = [log(.9), log(0.2), 1]; %Log of estimated parameters to fit hill function
    dataparms(34).Normp0 = [-.5, .8, 1];

    %% 100uM meAsp, LAsp Dose
    dataparms(35).Removal = 0;
    dataparms(35).backConc = 0.0;
    dataparms(35).xlabels = '[L-Asp]';
    dataparms(35).Lplot = 10.^[-3:0.01:2]; %Range of plots    
    dataparms(35).Files = ["./Data/230710_FOV1.mat"]; %Where to load data from
    dataparms(35).OutFolder = '100meAsp_LAspDR/'; %Where to output data to
    dataparms(35).concLevels = [.05, 0.1, 0.2, 0.4, 0.8];
    dataparms(35).Hillp0 = [log(.9), log(0.2), 1]; %Log of estimated parameters to fit hill function
    dataparms(35).Normp0 = [-.5, .8, 1];

    %% 0.01uM Lasp, meAsp Dose
    dataparms(36).Removal = 0;
    dataparms(36).backConc = 0;
    dataparms(36).xlabels = '[meAsp]';
    dataparms(36).Lplot = 10.^[-3:0.01:2]; %Range of plots    
    dataparms(36).Files = ["./Data/230711_FOV1.mat"]; %Where to load data from
    dataparms(36).OutFolder = '0_01LAsp_meAspDR/'; %Where to output data to
    dataparms(36).concLevels = [0.2, 0.5, 1, 2, 4];
    dataparms(36).Hillp0 = [log(.9), log(0.2), 1]; %Log of estimated parameters to fit hill function
    dataparms(36).Normp0 = [-.1, .8, 1];

    %% 100uM MeAsp + 10uM LAsp, LAsp Dose
    dataparms(37).Removal = 0;
    dataparms(37).backConc = 10;
    dataparms(37).xlabels = '[LAsp]';
    dataparms(37).Lplot = 10.^[-3:0.01:2]; %Range of plots    
    dataparms(37).Files = ["./Data/230714_FOV1.mat", "./Data/230715_FOV1.mat"]; %Where to load data from
    dataparms(37).OutFolder = '10LAsp_100meAsp_LAspDR_V2/'; %Where to output data to
    dataparms(37).concLevels = [10.1, 10.2, 10.5, 11, 12; 10.1, 10.2, 10.5, 11, 12];
    dataparms(37).Hillp0 = [log(.9), log(0.2), 1]; %Log of estimated parameters to fit hill function
    dataparms(37).Normp0 = [log(10.1), 2, 1];

    %% 0.01uM MeAsp, meAsp Dose
    dataparms(38).Removal = 0;
    dataparms(38).backConc = .01;
%     dataparms(38).backConc = .0;
    dataparms(38).xlabels = '[meAsp]';
    dataparms(38).Lplot = 10.^[-3:0.01:2]; %Range of plots    
%     dataparms(38).Files = ["./Data/230815_FOV1.mat", "./Data/230815_FOV2.mat",...
%         "./Data/230816_FOV1.mat", "./Data/230816_FOV2.mat"]; %Where to load data from
    dataparms(38).Files = ["./Data/230815_FOV1.mat", "./Data/230815_FOV2.mat"]; %Where to load data from
    dataparms(38).OutFolder = '0_01meAsp_meAspDR/'; %Where to output data to
    dataparms(38).concLevels = [0.21, 0.51, 1.01, 2.01, 4.01; ...
        0.21, 0.51, 1.01, 2.01, 4.01;...
        0.21, 0.51, 1.01, 2.01, 4.01; ...
        0.21, 0.51, 1.01, 2.01, 4.01];
    dataparms(38).Hillp0 = [log(.9), log(0.2), 1]; %Log of estimated parameters to fit hill function
    dataparms(38).Normp0 = [log(1.21), 2, 1];

    %% 10uM L-Asp + 100uM meAsp, meAsp dose
    dataparms(39).Removal = 0;
    dataparms(39).backConc = 100;
    dataparms(39).xlabels = '[LAsp]';
    dataparms(39).Lplot = 10.^[1:0.01:3]; %Range of plots    
    dataparms(39).Files = ["./Data/230822_FOV1.mat", "./Data/230823_FOV1.mat", "./Data/230823_FOV2.mat"]; %Where to load data from
    dataparms(39).OutFolder = '10LAsp_100meAsp_MeAspDR/'; %Where to output data to
    dataparms(39).concLevels = [105, 110, 120, 140, 180; 105, 110, 120, 140, 180; 105, 110, 120, 140, 180];
    dataparms(39).Hillp0 = [log(2), log(10^2), 1]; %Log of estimated parameters to fit hill function
    dataparms(39).Normp0 = [log(110), log(2), 1];

%     %% 10uM L-Asp + 100uM meAsp, meAsp dose as 0 background
%     dataparms(39).Removal = 0;
%     dataparms(39).backConc = 0;
%     dataparms(39).xlabels = '[LAsp]';
%     dataparms(39).Lplot = 10.^[0:0.01:2]; %Range of plots    
%     dataparms(39).Files = ["./Data/230822_FOV1.mat", "./Data/230823_FOV1.mat", "./Data/230823_FOV2.mat"]; %Where to load data from
%     dataparms(39).OutFolder = '10LAsp_100meAsp_MeAspDR/'; %Where to output data to
%     dataparms(39).concLevels = [5, 10, 20, 40, 80; 5, 10, 20, 40, 80; 05, 10, 20, 40, 80];
%     dataparms(39).Hillp0 = [log(2), log(10), 1]; %Log of estimated parameters to fit hill function
%     dataparms(39).Normp0 = [log(10), log(2), 1];

    %% 0.1uM L-Asp + 100uM meAsp, Lasp dose
    dataparms(40).Removal = 0;
    dataparms(40).backConc = 0.1;
    dataparms(40).xlabels = '[L-Asp]';
    dataparms(40).Lplot = 10.^[-3:0.01:2]; %Range of plots    
    dataparms(40).Files = ["./Data/231017_FOV1.mat", "./Data/231017_FOV2.mat", "./Data/231018_FOV1.mat", "./Data/231018_FOV2.mat"]; %Where to load data from
%     dataparms(40).Files = ["./Data/231017_FOV1.mat", "./Data/231017_FOV2.mat", "./Data/231018_FOV1.mat"]; %Where to load data from
    dataparms(40).OutFolder = '100meAsp_0_1LAsp_LAspDR/'; %Where to output data to
    dataparms(40).concLevels = [0.2, 0.3, 0.6, 1.1, 2.1; 0.2, 0.3, 0.6, 1.1, 2.1; 0.2, 0.3, 0.6, 1.1, 2.1; 0.2, 0.3, 0.6, 1.1, 2.1];
%     dataparms(35).concLevels = [0.02, 0.05, 0.1, 0.2, 0.4; 0.02, 0.05, 0.1, 0.2, 0.4];
    dataparms(40).Hillp0 = [log(.9), log(0.2), 1]; %Log of estimated parameters to fit hill function
    dataparms(40).Normp0 = [log(0.3), 2, 1];

    %% 10uM Serine background, L-Asp Dose Response
    dataparms(41).Removal = 0;
    dataparms(41).backConc = 0; %Background concentration for fitting hill function
    dataparms(41).xlabels = '[L-Asp]'; %xlabel on plots
    dataparms(41).Lplot = 10.^[-4:0.01:3]; %Range of plots
    dataparms(41).Files = ["./Data/231207_FOV1.mat", "./Data/231207_FOV2.mat", "./Data/231212_FOV1.mat", "./Data/231212_FOV2.mat"]; %Where to load data from
%     dataparms(41).Files = ["./Data/231212_FOV1.mat", "./Data/231212_FOV2.mat"]; %Where to load data from
    dataparms(41).OutFolder = '10serBackLAspDose/'; %Where to output data to
    dataparms(41).concLevels = [0.02, 0.05, 0.1, 0.2, 0.4; 0.02, 0.05, 0.1, 0.2, 0.4; 0.02, 0.05, 0.1, 0.2, 0.4; 0.02, 0.05, 0.1, 0.2, 0.4];
    dataparms(41).Hillp0 = [log(0.04), log(3), 1]; %Log of estimated parameters to fit hill function
    dataparms(41).Normp0 = [log(0.03), 1, 1];
 
    %% 0 background, Maltose Dose Response
    dataparms(42).Removal = 0;
    dataparms(42).backConc = 0; %Background concentration for fitting hill function
    dataparms(42).xlabels = '[Maltose]'; %xlabel on plots
    dataparms(42).Lplot = 10.^[-4:0.01:3]; %Range of plots
%     dataparms(42).Files = ["./Data/210610_FOV1.mat", "./Data/210610_FOV2.mat"]; %Where to load data from
    dataparms(42).Files = ["./Data/240131_FOV1.mat", "./Data/240131_FOV2.mat", "./Data/240201_FOV1.mat", "./Data/240201_FOV2.mat"];
    dataparms(42).OutFolder = '0backMaltoseDose/'; %Where to output data to
%     dataparms(42).concLevels = [0.0316, 0.1, 0.316, 1, 1.316; 0.0316, 0.1, 0.316, 1, 1.316];
    dataparms(42).concLevels = [0.05, 0.1, 0.2, 0.5, 1; 0.05, 0.1, 0.2, 0.5, 1; 0.05, 0.1, 0.2, 0.5, 1; 0.05, 0.1, 0.2, 0.5, 1];
    dataparms(42).Hillp0 = [log(1.2), log(0.3), 1]; %Log of estimated parameters to fit hill function [n, K]
    dataparms(42).Normp0 = [log(0.3), 1, 1]; %[log(mean), sigma]

    %% 1uM Maltose background, Maltose Dose Response
    dataparms(43).Removal = 0;
    dataparms(43).backConc = 1; %Background concentration for fitting hill function
    dataparms(43).xlabels = '[Maltose]'; %xlabel on plots
    dataparms(43).Lplot = 10.^[-4:0.01:3]; %Range of plots
    dataparms(43).Files = ["./Data/210616_FOV1.mat", "./Data/210616_FOV2.mat", "./Data/240202_FOV1.mat", "./Data/240202_FOV2.mat"]; %Where to load data from
    dataparms(43).OutFolder = '1MaltosebackMaltoseDose/'; %Where to output data to
    dataparms(43).concLevels = [1.1, 1.2, 1.4, 1.8, 2.0; 1.1, 1.2, 1.4, 1.8, 2.0; 1.1, 1.2, 1.4, 1.8, 2.0; 1.1, 1.2, 1.4, 1.8, 2.0];
    dataparms(43).Hillp0 = [log(2), log(1.2), 1]; %Log of estimated parameters to fit hill function [n, K]
    dataparms(43).Normp0 = [log(1.2), 1, 1]; %[log(mean), sigma]

    %% 1uM Maltose background, meAsp Dose Response
    dataparms(44).Removal = 0;
    dataparms(44).backConc = 0; %Background concentration for fitting hill function
    dataparms(44).xlabels = '[MeAsp]'; %xlabel on plots
    dataparms(44).Lplot = 10.^[-4:0.01:3]; %Range of plots
%     dataparms(44).Files = ["./Data/240129_FOV1.mat", "./Data/240129_FOV2.mat", "./Data/240130_FOV1.mat", "./Data/240130_FOV2.mat"]; %Where to load data from
    dataparms(44).Files = ["./Data/240129_FOV1.mat", "./Data/240129_FOV2.mat", "./Data/240130_FOV1.mat"]; %Where to load data from
    dataparms(44).OutFolder = '1MaltosebackMeAspDose/'; %Where to output data to
    dataparms(44).concLevels = [0.2, 0.5, 1, 2, 4; 0.2, 0.5, 1, 2, 4; 0.2, 0.5, 1, 2, 4; 0.2, 0.5, 1, 2, 4];
    dataparms(44).Hillp0 = [log(2), log(1.1), 1]; %Log of estimated parameters to fit hill function [n, K]
    dataparms(44).Normp0 = [log(1.3), 2, 1]; %[log(mean), sigma]

    %% 100uM MeAsp + 1uM Maltose Background, MeAsp Dose Response
    dataparms(45).Removal = 0;
    dataparms(45).backConc = 100; %Background concentration for fitting hill function
    dataparms(45).xlabels = '[MeAsp]'; %xlabel on plots
    dataparms(45).Lplot = 10.^[1:0.01:3]; %Range of plots
    dataparms(45).Files = ["./Data/240229_FOV1.mat", "./Data/240229_FOV1.mat"];
    dataparms(45).OutFolder = '100MeAsp_1uMMaltose_BackMeAspDR/'; %Where to output data to   
    dataparms(45).concLevels = [102, 105, 110, 120, 140; 102, 105, 110, 120, 140; 102, 105, 110, 120, 140];
    dataparms(45).Hillp0 = [log(2), log(105), 1]; %Log of estimated parameters to fit hill function
    dataparms(45).Normp0 = [log(103), 2, 1];

end