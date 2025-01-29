% guess = zeros(sizeY, sizeX, 13);
% guess(:,:,1)=pkf(:,:,1)/1e9;
% guess(:,:,2)=pkf(:,:,2)/1e9;
% guess(:,:,3)=0.001; %prev -0.01
% 
% guess(:,:,4)=1/3*(pkf(:,:,5)/1e9-pkf(:,:,6)/1e9);
% guess(:,:,5)=2/3*(pkf(:,:,5)/1e9-pkf(:,:,6)/1e9);
% guess(:,:,6)=1/3*(pkf(:,:,5)/1e9+pkf(:,:,6)/1e9);
% guess(:,:,7)=2/3*(pkf(:,:,5)/1e9+pkf(:,:,6)/1e9);
% 
% guess(:,:,8)=pkf(:,:,3)/1e9; %bc peakfinder already converts FWHM to std
% guess(:,:,9)=pkf(:,:,4)/1e9;
% guess(:,:,10)=pkf(:,:,3)/1e9; 
% guess(:,:,11)=pkf(:,:,4)/1e9;
% 
% guess(:,:,12) = 1;
% guess(:,:,13) = 0;

constraint = zeros(sizeY, sizeX, 26);
constraint(:,:,1) = pkf(:,:,1)/1e9-0.1;
constraint(:,:,2) = pkf(:,:,1)/1e9+0.1;

constraint(:,:,3) = pkf(:,:,2)/1e9-0.1;
constraint(:,:,4) = pkf(:,:,2)/1e9+0.1;

constraint(:,:,5) = -0.02;
constraint(:,:,6) = 0.02; %prev 0.05 span

constraint(:,:,7) = -0.03;
constraint(:,:,8) = -0.0001;

constraint(:,:,9) = -0.03; %for g1 the height gaps can be up to 0.015
constraint(:,:,10) = -0.0001;

constraint(:,:,11) = -0.03;
constraint(:,:,12) = -0.0001;

constraint(:,:,13) = -0.03; %made it rlly small (prev 0003
constraint(:,:,14) = -0.0001;

constraint(:,:,15) = 0.003;
constraint(:,:,16) = 0.025;

constraint(:,:,17) = -0.01;
constraint(:,:,18) = 0.01;

constraint(:,:,19) = 0.003;
constraint(:,:,20) = 0.025;

constraint(:,:,21) = -0.01;
constraint(:,:,22) = 0.01;

constraint(:,:,23) = 0.98;
constraint(:,:,24) = 1.02;

constraint(:,:,25) = -0.002;
constraint(:,:,26) = 0.002;


% constraint = zeros(sizeY, sizeX, 26);
% constraint(:,:,1) = 2.85;
% constraint(:,:,2) = 3;
% constraint(:,:,3) = 0;
% constraint(:,:,4) = 0.4;
% constraint(:,:,5) = 0;
% constraint(:,:,6) = 0.05;
% constraint(:,:,7) = -0.01;
% constraint(:,:,8) = 0.01;
% constraint(:,:,9) = -0.04;
% constraint(:,:,10) = 0;
% constraint(:,:,11) = -0.01;
% constraint(:,:,12) = 0.01;
% constraint(:,:,13) = 0;
% constraint(:,:,14) = 0.01;
% constraint(:,:,15) = -0.01;
% constraint(:,:,16) = 0.01;
% constraint(:,:,17) = 0;
% constraint(:,:,18) = 0.04;
% constraint(:,:,19) = -0.01;
% constraint(:,:,20) = 0.01;
% constraint(:,:,21) = -0.04;
% constraint(:,:,22) = 0;
% constraint(:,:,23) = -0.01;
% constraint(:,:,24) = 0.01;
% constraint(:,:,25) = 0.98;
% constraint(:,:,26) = 1.02;