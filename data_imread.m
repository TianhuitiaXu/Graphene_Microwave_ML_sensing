%% DC data + DC plot
% Graphene 1,3,4,6,9,10 -- 6 devices
% 0,9e-6,...5.4e-3 -- 9 concentrations
% regression: 9 * 122*6
% classification: 9*6 * 122
concen_num  = 9;
S_DC_0 = [Gr1_DI_water(:,1);Gr3_DI_water(:,1);Gr4_DI_water(:,1);Gr6_DI_water(:,1);...
    Gr9_DI_water(:,1);Gr10_DI_water(:,1);];
S_DC_1 = [Gr1_9ug_dL(:,1);Gr3_9ug_dL(:,1);Gr4_9ug_dL(:,1);Gr6_9ug_dL(:,1);...
    Gr9_9ug_dL(:,1);Gr10_9ug_dL(:,1);];
S_DC_4 = [Gr1_54ug_dL(:,1);Gr3_54ug_dL(:,1);Gr4_54ug_dL(:,1);Gr6_54ug_dL(:,1);...
    Gr9_54ug_dL(:,1);Gr10_54ug_dL(:,1);];
S_DC_6 = [Gr1_540ug_dL(:,1);Gr3_540ug_dL(:,1);Gr4_540ug_dL(:,1);Gr6_540ug_dL(:,1);...
    Gr9_540ug_dL(:,1);Gr10_540ug_dL(:,1);];
S_DC_8 = [Gr1_5400ug_dL(:,1);Gr3_5400ug_dL(:,1);Gr4_5400ug_dL(:,1);Gr6_5400ug_dL(:,1);...
    Gr9_5400ug_dL(:,1);Gr10_5400ug_dL(:,1);];
S_DC_3 = [Gr1_36ug_dL(:,1);Gr3_36ug_dL(:,1);Gr4_36ug_dL(:,1);Gr6_36ug_dL(:,1);...
    Gr9_36ug_dL(:,1);Gr10_36ug_dL(:,1);];
S_DC_2 = [Gr1_18ug_dL(:,1);Gr3_18ug_dL(:,1);Gr4_18ug_dL(:,1);Gr6_18ug_dL(:,1);...
    Gr9_18ug_dL(:,1);Gr10_18ug_dL(:,1);];
S_DC_5 = [Gr1_180ug_dL(:,1);Gr3_180ug_dL(:,1);Gr4_180ug_dL(:,1);Gr6_180ug_dL(:,1);...
    Gr9_180ug_dL(:,1);Gr10_180ug_dL(:,1);];
S_DC_7 = [Gr1_1800ug_dL(:,1);Gr3_1800ug_dL(:,1);Gr4_1800ug_dL(:,1);Gr6_1800ug_dL(:,1);...
    Gr9_1800ug_dL(:,1);Gr10_1800ug_dL(:,1);];
S_DC = [S_DC_0';S_DC_1';S_DC_2';S_DC_3';S_DC_4';S_DC_5';S_DC_6';S_DC_7';S_DC_8'];
Vg_water = Gr1_DI_water(62:122,2);
Vg_glucose = Gr1_9ug_dL(62:122,2);

% S_DC_clas = zeros(9*6,122);
% for i = 1:6
%     S_DC_clas(1+(i-1)*9:9*i,:) = S_DC(:,1+122*(i-1):122*i);
% end
% 
x = Gr1_9ug_dL(:,2);
pp = cell(1,9);
Y = zeros(1,9);
for j = 1%:6
%     figure(1)
    for i = 1:9
%         y = 20*log10(S_DC(i,62+(j-1)*122:122*j)*1e6);
        y = 20*log10(S_DC(i,62+(j-1)*122:122*j)*1e6./ Vg_water');
        Y(i) = y(11);
        pp{1,i} = plot(x(62:122)',y,'-o','MarkerIndices',1:10:length(y),'linewidth',2);%'MarkerSize',6 % x(61:122)',
        xlabel('Vg (V)','FontSize',12);%Frequency (GHz) Concentrations (g/dl)
        ylabel('I (dBμ)','FontSize',12);  
        legend([pp{1,1} pp{1,2} pp{1,3} pp{1,4} pp{1,5} pp{1,6} pp{1,7} pp{1,8} pp{1,9}],{'0' '9e-6' '1.8e-5' '3.6e-5' '5.4e-5' '1.8e-4'...
            '5.4e-4' '1.8e-3' '5.4e-3'},  ...
        'location','NE','FontSize',12);
        grid on
        set(gcf,'unit','normalized','position',[0.2,0.2,0.57,0.43]);
%         axis square
        hold on
    end
    hold off
end
figure(2)
S_target_expe = [0;9e-6;1.8e-5;3.6e-5;5.4e-5;1.8e-4;5.4e-4;1.8e-3;5.4e-3];%
plot(S_target_expe,Y,'-o','MarkerSize',6,'linewidth',2);
xlabel('Concentrations (g/dl)','FontSize',12);%Frequency (GHz) Concentrations (g/dl)
ylabel('I (dBμ) at Vg = 0.595V','FontSize',12); 
ax = gca; % current axes
ax.FontSize = 12;
% axis square

% % 
% c = polyfit(x_concen,s21,1);
% % Display evaluated equation y = m*x + b
% disp(['Equation is y = ' num2str(c(1)) '*x + ' num2str(c(2))])
% % Evaluate fit equation using polyval
% y_est = polyval(c,x_concen);
% % Add trend line to plot
% hold on
% plot(x_concen,y_est,'r--','LineWidth',2)
% hold off
%% Simulation data with different chemical potential brought by concentrations
% for all different concentrations using data at 0V
clear all
sample_num = 4004;
S21 = zeros(10,1001);
S11 = zeros(10,1001);
S_all = zeros(10,sample_num);
S_test = zeros(10,sample_num);
S_target = [0;0.025;0.05;0.1;0.15;0.2;0.25;0.3;0.4;0.5];
j=1;
for i =  1:9% devices 56:65(used on the data of 0-0.103268)
%    for j=1:2
    ad1 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\simulation-data\0.014\0.04-0.16914_',int2str(i),'.s2p'];%int2str(j),-0.07-0.04;0.5gdl_chemicalpotential0.0951_
%     ad1 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\simulation-data\0.04\10chemical-potential_0-0.103268_',int2str(i),'.s2p'];
    
%     AD = {ad1,ad2,ad3,ad4,ad5,ad6,ad7,ad8,ad9,ad10};%,ad7,ad8,ad9,ad10,ad11,ad12,ad13,ad14,ad15,ad16,ad17,ad18,ad19,ad20,ad13,ad14,ad15,ad16,ad17,ad18,ad19
%     for j = 1:10 % concentrations
        s = sparameters(ad1);
%         C{j,1} = s; 
        s11 = rfparam(s,1,1);
        s11_re = real(s11);
        s11_imag = imag(s11);
        s11_amp = abs(s11);
        s11_phase = angle(s11);
        s21 = rfparam(s,2,1);
        s21_re = real(s21);
        s21_imag = imag(s21);
        s21_amp = abs(s21);
        s21_phase = angle(s21);
%          if j == 1
            s21_db = 20*log10(abs(s21));
            S21(i,:) = s21_db;%-55
            s11_db = 20*log10(abs(s11));
            S11(i,:) = s11_db;%-55
%          end
        ss=cat(1,s11_imag,s11_re,s21_imag,s21_re);

%         if j<2
            S_all(i,1+sample_num*(j-1):sample_num*j) = ss';%-55
%         else
%             S_test(i,1+sample_num*(j-2):sample_num*(j-1)) = ss';%-55
%         end
%     end
end  
% ans=fft(S21(1,:));L=length(ans);P2 = abs(ans/L);
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% stem(P1);
%% build experimental data for training, validation, and testing
clear all
concen_num = 12;%28
test_num = 0;
S_voltage = zeros(concen_num,18*4004);% 3 concentrations, 18 voltages
S_voltage_test = zeros(concen_num,18*4004);% 3 concentrations, 18 voltages
S21_Voltage = cell(1,concen_num);
S21_voltage = zeros(18,1001);
% S_voltage_Target = [0.025;0.05;0.1;0.2;0.25;0.3;0.4;0.5;0.75;1;1.5;2;4;6;8;10;12;16];

S_all_expe = zeros(concen_num-test_num,4004);
S_test_expe = zeros(concen_num,4004);%*5
S_test_expe1 = zeros(concen_num,4004);

% S_all_expe_devi169 = zeros(concen_num,4004*3);
% 
% S_test_expe_devi169 = zeros(concen_num,4004*3);
% 
% for each concentration, 8 * 3(different devices) * 18(voltages) = 432
S_all_expe_clas = zeros(concen_num*5,4004);
S_test_expe_clas = zeros(concen_num*2,4004);

% S_target_expe = [0;0.025;0.05;0.1;0.15;0.2;0.25;0.3;0.4;0.5];%;1;1.5;2;4;6;8;12;16
% S_target_expe = [0;9e-6;1.8e-5;3.6e-5;5.4e-5;1.8e-4;5.4e-4;1.8e-3;5.4e-3;0.025;0.05;0.1;0.15;0.2;0.25;0.3;0.4;0.5;0.75;1;1.5;2;4;6;8;10;12;16];
S_target_expe = [0;0.5;0.75;1;1.5;2;4;6;8;10;12;16];
S_test_target = [0.05;0.5;4];%3.6e-5;1.8e-3;0.05;0.5;6
S_train_target = [0;0.025;0.05;0.1;0.15;0.2;0.25;0.3;0.4;0.5];%[0;9e-6;1.8e-5;3.6e-5;5.4e-5;1.8e-4];%;5.4e-4;1.8e-3;5.4e-3;0.025;0.05;0.1;0.15;0.2;0.25;...
%      0.3;0.4;0.5;0.75;1;1.5;2;4;6;8;10;12;16];%9e-6;1.8e-5;5.4e-5;1.8e-4];%;...
% S_target_expe = [0;9e-6;1.8e-5;3.6e-5;5.4e-5;1.8e-4];
% Target_clas = linspace(1,15,15);
S_target_expe_clas = zeros(concen_num*5,1);
% Target_clas_whole = zeros(75,1);

% classification Target
for m = 1:5
    S_target_expe_clas((m-1)*concen_num+1:m*concen_num) = S_target_expe;
%     Target_clas_whole((m-1)*15+1:m*15) = Target_clas;
end
C = cell(1,6);
AD = cell(1,15);

% %using data only at 0.28v and 0.32v
% S = zeros(concen_num,4004);
% S_test = zeros(concen_num,4004);
% for i = 1: concen_num
%     for j = 1:2 % voltages
%         for k = 1:2 % train/test
%             if j == 1
%                 ad0 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\0new\Graphene1_no',int2str(k),'_0.28Vdown.s2p'];
% %                 ad1 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\9e-06\Graphene1_no',int2str(k),'_0.28V.s2p'];
% %                 ad2 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\1.8e-05\Graphene1_no',int2str(k),'_0.28V.s2p'];
% %                 ad3 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\3.6e-05\Graphene1_no',int2str(k),'_0.28V.s2p'];
% %                 ad4 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\5.4e-05\Graphene1_no',int2str(k),'_0.28V.s2p'];
% %                 ad5 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\1.8e-04\Graphene1_no',int2str(k),'_0.28V.s2p'];
% %                 ad6 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\0.025\Graphene1_no',int2str(k),'_0.28Vdown.s2p'];
% %                 ad7 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\0.05\Graphene1_no',int2str(k),'_0.28Vdown.s2p'];
% %                 ad8 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\0.1\Graphene1_no',int2str(k),'_0.28Vdown.s2p'];
% %                 ad9 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\0.15\Graphene1_no',int2str(k),'_0.28Vdown.s2p'];
% %                 ad10 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\0.2\Graphene1_no',int2str(k),'_0.28Vdown.s2p'];
% %                 ad11 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\0.25\Graphene1_no',int2str(k),'_0.28Vdown.s2p'];
% %                 ad12 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\0.3\Graphene1_no',int2str(k),'_0.28Vdown.s2p'];
% %                 ad13 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\0.4\Graphene1_no',int2str(k),'_0.28Vdown.s2p'];
%                 ad1 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\0.5\Graphene1_no',int2str(k),'_0.28Vdown.s2p'];
%                 ad9 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\0.75\Graphene1_no',int2str(k),'_0.28Vdown.s2p'];
%                 ad10 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\1\Graphene1_no',int2str(k),'_0.28Vdown.s2p'];
%                 ad11 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\1.5\Graphene1_no',int2str(k),'_0.28Vdown.s2p'];
%                 ad12 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\2\Graphene1_no',int2str(k),'_0.28Vdown.s2p'];
%                 ad13 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\4\Graphene1_no',int2str(k),'_0.28Vdown.s2p'];
%                 ad14 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\6\Graphene1_no',int2str(k),'_0.28Vdown.s2p'];
%                 ad15 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\8\Graphene1_no',int2str(k),'_0.28Vdown.s2p'];
%                 ad16 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\10\Graphene1_no',int2str(k),'_0.28Vdown.s2p'];
%                 ad17 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\12\Graphene1_no',int2str(k),'_0.28Vdown.s2p'];
%                 ad18 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\16\Graphene1_no',int2str(k),'_0.28Vdown.s2p'];
%             else
%                 ad0 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\0new\Graphene1_no',int2str(k),'_0.32Vdown.s2p'];
% %                 ad1 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\9e-06\Graphene1_no',int2str(k),'_0.0V.s2p'];
% %                 ad2 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\1.8e-05\Graphene1_no',int2str(k),'_0.0V.s2p'];
% %                 ad3 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\3.6e-05\Graphene1_no',int2str(k),'_0.0V.s2p'];
% %                 ad4 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\5.4e-05\Graphene1_no',int2str(k),'_0.0V.s2p'];
% %                 ad5 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\1.8e-04\Graphene1_no',int2str(k),'_0.0V.s2p'];
% %                 ad6 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\0.025\Graphene1_no',int2str(k),'_0.0V.s2p'];
% %                 ad7 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\0.05\Graphene1_no',int2str(k),'_0.0V.s2p'];
% %                 ad8 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\0.1\Graphene1_no',int2str(k),'_0.0V.s2p'];
% %                 ad9 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\0.15\Graphene1_no',int2str(k),'_0.0V.s2p'];
% %                 ad10 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\0.2\Graphene1_no',int2str(k),'_0.0V.s2p'];
% %                 ad11 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\0.25\Graphene1_no',int2str(k),'_0.0V.s2p'];
% %                 ad12 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\0.3\Graphene1_no',int2str(k),'_0.0V.s2p'];
% %                 ad13 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\0.4\Graphene1_no',int2str(k),'_0.0V.s2p'];
%                 ad1 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\0.5\Graphene1_no',int2str(k),'_0.32Vdown.s2p'];
%                 ad9 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\0.75\Graphene1_no',int2str(k),'_0.32Vdown.s2p'];
%                 ad10 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\1\Graphene1_no',int2str(k),'_0.32Vdown.s2p'];
%                 ad11 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\1.5\Graphene1_no',int2str(k),'_0.32Vdown.s2p'];
%                 ad12 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\2\Graphene1_no',int2str(k),'_0.32Vdown.s2p'];
%                 ad13 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\4\Graphene1_no',int2str(k),'_0.32Vdown.s2p'];
%                 ad14 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\6\Graphene1_no',int2str(k),'_0.32Vdown.s2p'];
%                 ad15 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\8\Graphene1_no',int2str(k),'_0.32Vdown.s2p'];
%                 ad16 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\10\Graphene1_no',int2str(k),'_0.32Vdown.s2p'];
%                 ad17 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\12\Graphene1_no',int2str(k),'_0.32Vdown.s2p'];
%                 ad18 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\16\Graphene1_no',int2str(k),'_0.32Vdown.s2p'];
%             end
%             AD = {ad0,ad1,ad9,ad10,ad11,ad12,ad13,ad14,ad15,ad16,ad17,ad18};%,ad9,ad10,ad11,ad12,ad13,ad14,ad15,ad16,ad17,ad18
%             s = sparameters(AD{1,i});
%             s11 = rfparam(s,1,1);
%             s11_re = real(s11);
%             s11_imag = imag(s11);
%             s11_amp = abs(s11);
%             s11_phase = angle(s11);
%             s21 = rfparam(s,2,1);
%             s21_re = real(s21);
%             s21_imag = imag(s21);
%             s21_amp = abs(s21);
%             s21_phase = angle(s21);
%             ss=cat(1,s11_imag,s11_re,s21_imag,s21_re);
%             if k==1%k == 4 || k==5
%                 S(i,1+4004*(j-1):4004*j) = ss;
%             else
% %             elseif k==1 
% %                 S_test(i,1:4004) = ss;
% %             elseif k==8
% %                 S_test(i,4005:8008) = ss;
% %             end
% %             else
%                 S_test(i,1+4004*(j-1):4004*j) = ss;
%             end
%         end
%     end
% end
    
% % with 18 voltage change
% voltage = linspace(0,0.68,18);
% voltage = double(vpa(sym(voltage),2));
% for i = 1:18
%     for j = 1:2
%         if i == 1
%             ad0 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\0new\Graphene1_no',int2str(j),'_0.0V.s2p'];
% %             ad1 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\0.025\Graphene1_no',int2str(j),'_0.0V.s2p'];
% %             ad2 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\0.05\Graphene1_no',int2str(j),'_0.0V.s2p'];
% %             ad3 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\0.1\Graphene1_no',int2str(j),'_0.0V.s2p'];
% %             ad34 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\0.15\Graphene1_no',int2str(j),'_0.0V.s2p'];
% %             ad4 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\0.2\Graphene1_no',int2str(j),'_0.0V.s2p'];
% %             ad5 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\0.25\Graphene1_no',int2str(j),'_0.0V.s2p'];
% %             ad6 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\0.3\Graphene1_no',int2str(j),'_0.0V.s2p'];
% %             ad7 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\0.4\Graphene1_no',int2str(j),'_0.0V.s2p'];
% %             ad8 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\0.5\Graphene1_no',int2str(j),'_0.0V.s2p'];
% %             ad9 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\0.75\Graphene1_no',int2str(j),'_0.0V.s2p'];
% %             ad10 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\1\Graphene1_no',int2str(j),'_0.0V.s2p'];
% %             ad11 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\1.5\Graphene1_no',int2str(j),'_0.0V.s2p'];
% %             ad12 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\2\Graphene1_no',int2str(j),'_0.0V.s2p'];
% %             ad13 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\4\Graphene1_no',int2str(j),'_0.0V.s2p'];
% %             ad14 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\6\Graphene1_no',int2str(j),'_0.0V.s2p'];
% %             ad15 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\8\Graphene1_no',int2str(j),'_0.0V.s2p'];
% %             ad16 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\10\Graphene1_no',int2str(j),'_0.0V.s2p'];
% %             ad17 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\12\Graphene1_no',int2str(j),'_0.0V.s2p'];
% %             ad18 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\16\Graphene1_no',int2str(j),'_0.0V.s2p'];
%             ad19 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\9e-06\Graphene1_no',int2str(j),'_0.0V.s2p'];
%             ad20 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\1.8e-05\Graphene1_no',int2str(j),'_0.0V.s2p'];
%             ad21 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\3.6e-05\Graphene1_no',int2str(j),'_0.0V.s2p'];
%             ad22 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\5.4e-05\Graphene1_no',int2str(j),'_0.0V.s2p'];
%             ad23 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\1.8e-04\Graphene1_no',int2str(j),'_0.0V.s2p'];
%             ad24 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\5.4e-04\Graphene1_no',int2str(j),'_0.0V.s2p'];
% %             ad25 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\1.8e-03\Graphene1_no',int2str(j),'_0.0V.s2p'];
% %             ad26 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\5.4e-03\Graphene1_no',int2str(j),'_0.0V.s2p'];
%         else
%             ad0 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\0new\Graphene1_no',int2str(j),'_',num2str(voltage(i)),'V.s2p'];
% %             ad1 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\0.025\Graphene1_no',int2str(j),'_',num2str(voltage(i)),'V.s2p'];
% %             ad2 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\0.05\Graphene1_no',int2str(j),'_',num2str(voltage(i)),'V.s2p'];
% %             ad3 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\0.1\Graphene1_no',int2str(j),'_',num2str(voltage(i)),'V.s2p'];
% %             ad34 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\0.15\Graphene1_no',int2str(j),'_',num2str(voltage(i)),'V.s2p'];
% %             ad4 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\0.2\Graphene1_no',int2str(j),'_',num2str(voltage(i)),'V.s2p'];
% %             ad5 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\0.25\Graphene1_no',int2str(j),'_',num2str(voltage(i)),'V.s2p'];
% %             ad6 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\0.3\Graphene1_no',int2str(j),'_',num2str(voltage(i)),'V.s2p'];
% %             ad7 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\0.4\Graphene1_no',int2str(j),'_',num2str(voltage(i)),'V.s2p'];
% %             ad8 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\0.5\Graphene1_no',int2str(j),'_',num2str(voltage(i)),'V.s2p'];
% %             ad9 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\0.75\Graphene1_no',int2str(j),'_',num2str(voltage(i)),'V.s2p'];
% %             ad10 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\1\Graphene1_no',int2str(j),'_',num2str(voltage(i)),'V.s2p'];
% %             ad11= ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\1.5\Graphene1_no',int2str(j),'_',num2str(voltage(i)),'V.s2p'];
% %             ad12 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\2\Graphene1_no',int2str(j),'_',num2str(voltage(i)),'V.s2p'];
% %             ad13 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\4\Graphene1_no',int2str(j),'_',num2str(voltage(i)),'V.s2p'];
% %             ad14 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\6\Graphene1_no',int2str(j),'_',num2str(voltage(i)),'V.s2p'];
% %             ad15 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\8\Graphene1_no',int2str(j),'_',num2str(voltage(i)),'V.s2p'];
% %             ad16 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\10\Graphene1_no',int2str(j),'_',num2str(voltage(i)),'V.s2p'];
% %             ad17 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\12\Graphene1_no',int2str(j),'_',num2str(voltage(i)),'V.s2p'];
% %             ad18 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\16\Graphene1_no',int2str(j),'_',num2str(voltage(i)),'V.s2p'];
%             ad19 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\9e-06\Graphene1_no',int2str(j),'_',num2str(voltage(i)),'V.s2p'];
%             ad20 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\1.8e-05\Graphene1_no',int2str(j),'_',num2str(voltage(i)),'V.s2p'];
%             ad21 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\3.6e-05\Graphene1_no',int2str(j),'_',num2str(voltage(i)),'V.s2p'];
%             ad22 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\5.4e-05\Graphene1_no',int2str(j),'_',num2str(voltage(i)),'V.s2p'];
%             ad23 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\1.8e-04\Graphene1_no',int2str(j),'_',num2str(voltage(i)),'V.s2p'];
% %             ad24 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\5.4e-04\Graphene1_no',int2str(j),'_',num2str(voltage(i)),'V.s2p'];
% %             ad25 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\1.8e-03\Graphene1_no',int2str(j),'_',num2str(voltage(i)),'V.s2p'];
% %             ad26 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\5.4e-03\Graphene1_no',int2str(j),'_',num2str(voltage(i)),'V.s2p'];
%         end
%         AD = {ad0,ad19,ad20,ad21,ad22,ad23};%,ad5,ad6,ad7,ad8,ad9,ad10,ad11,ad12,ad13,ad14,ad15,ad16,ad17,ad18,ad19,ad20,ad21,ad22,ad23
%         for k = 1:concen_num
%             s = sparameters(AD{1,k});
%             s11 = rfparam(s,1,1);
%             s11_re = real(s11);
%             s11_imag = imag(s11);
%             s11_amp = abs(s11);
%             s11_phase = angle(s11);
%             s21 = rfparam(s,2,1);
%             s21_re = real(s21);
%             s21_imag = imag(s21);
%             s21_amp = abs(s21);
%             s21_phase = angle(s21);
%             
%             ss=cat(1,s11_imag,s11_re,s21_imag,s21_re);
%             if j==1
%                 S_voltage(k,1+4004*(i-1):4004*i) = ss';
%             elseif j==2
%                 S_voltage_test(k,1+4004*(i-1):4004*i) = ss';
%             end                     
%         end
%     end
% end

% % for 3 different devices at 0V
% device_index = [1;6;9];
% % S21 = zeros(10,1001);
% for i = 1:3 % different graphene overlapping devices
%     for j = 1:6
%         ad0 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\0new\Graphene',int2str(device_index(i)),'_no',int2str(j),'_0.0V.s2p'];
% %         ad1 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\0.025\Graphene',int2str(device_index(i)),'_no',int2str(j),'_0.0V.s2p'];
% %         ad2 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\0.05\Graphene',int2str(device_index(i)),'_no',int2str(j),'_0.0V.s2p'];
% %         ad3 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\0.1\Graphene',int2str(device_index(i)),'_no',int2str(j),'_0.0V.s2p'];
% %         ad4 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\0.15\Graphene',int2str(device_index(i)),'_no',int2str(j),'_0.0V.s2p'];
% %         ad5 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\0.2\Graphene',int2str(device_index(i)),'_no',int2str(j),'_0.0V.s2p'];
% %         ad6 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\0.25\Graphene',int2str(device_index(i)),'_no',int2str(j),'_0.0V.s2p'];
% %         ad7 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\0.3\Graphene',int2str(device_index(i)),'_no',int2str(j),'_0.0V.s2p'];
% %         ad8 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\0.4\Graphene',int2str(device_index(i)),'_no',int2str(j),'_0.0V.s2p'];
% %         ad9 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\0.5\Graphene',int2str(device_index(i)),'_no',int2str(j),'_0.0V.s2p'];
% %         ad10 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\0.75\Graphene',int2str(device_index(i)),'_no',int2str(j),'_0.0V.s2p'];
% %         ad11 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\1\Graphene',int2str(device_index(i)),'_no',int2str(j),'_0.0V.s2p'];
% %         ad12 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\1.5\Graphene',int2str(device_index(i)),'_no',int2str(j),'_0.0V.s2p'];
% %             ad12 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\2\Graphene',int2str(device_index(i)),'_no',int2str(j),'_0.0V.s2p'];
% %             ad13 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\4\Graphene',int2str(device_index(i)),'_no',int2str(j),'_0.0V.s2p'];
% %             ad14 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\6\Graphene',int2str(device_index(i)),'_no',int2str(j),'_0.0V.s2p'];
% %             ad15 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\8\Graphene',int2str(device_index(i)),'_no',int2str(j),'_0.0V.s2p'];
% %             ad16 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\10\Graphene',int2str(device_index(i)),'_no',int2str(j),'_0.0V.s2p'];
% %             ad17 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\12\Graphene',int2str(device_index(i)),'_no',int2str(j),'_0.0V.s2p'];
% %             ad18 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\16\Graphene',int2str(device_index(i)),'_no',int2str(j),'_0.0V.s2p'];
%         ad2 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\9e-06\Graphene',int2str(device_index(i)),'_no',int2str(j),'_0.0V.s2p'];
%         ad3 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\1.8e-05\Graphene',int2str(device_index(i)),'_no',int2str(j),'_0.0V.s2p'];
%         ad4 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\3.6e-05\Graphene',int2str(device_index(i)),'_no',int2str(j),'_0.0V.s2p'];
%         ad5= ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\5.4e-05\Graphene',int2str(device_index(i)),'_no',int2str(j),'_0.0V.s2p'];
%         ad6 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\1.8e-04\Graphene',int2str(device_index(i)),'_no',int2str(j),'_0.0V.s2p'];
%     AD = {ad0,ad2,ad3,ad4,ad5,ad6};%,ad9,ad10,ad11,ad12,ad13,ad14,ad15,ad16,ad17,ad18
%         for k = 1:concen_num
%             s = sparameters(AD{1,k});
%             s11 = rfparam(s,1,1);
%             s11_re = real(s11);
%             s11_imag = imag(s11);
%             s11_amp = abs(s11);
%             s11_phase = angle(s11);
%             s21 = rfparam(s,2,1);
%             s21_re = real(s21);
%             s21_imag = imag(s21);
%             s21_amp = abs(s21);
%             s21_phase = angle(s21);
%             ss=cat(1,s11_imag,s11_re,s21_imag,s21_re);
% %             if i==3
% %                 s21_db = 20*log10(abs(s21));
% %                 S21(k,:) = s21_db;
% %             end
%             if j==1
%                 S_all_expe_devi169(k,1+4004*(i-1):4004*i) = ss';
%             elseif j==6
%                 S_test_expe_devi169(k,1+4004*(i-1):4004*i) = ss';
%             end                     
%         end
%     end
% end

% for all different concentrations using data at 0V
S21 = zeros(concen_num,1001);
for i = 1:8 % devices
    ad1 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\0new\Graphene1_no',int2str(i),'_0.0V.s2p'];
    ad2 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\9e-06\Graphene1_no',int2str(i),'_0.0V.s2p'];
    ad3 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\1.8e-05\Graphene1_no',int2str(i),'_0.0V.s2p'];
    ad4 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\3.6e-05\Graphene1_no',int2str(i),'_0.0V.s2p'];
    ad5= ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\5.4e-05\Graphene1_no',int2str(i),'_0.0V.s2p'];
    ad6 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\1.8e-04\Graphene1_no',int2str(i),'_0.0V.s2p'];
    ad7 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\5.4e-04\Graphene1_no',int2str(i),'_0.0V.s2p'];
    ad8 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\1.8e-03\Graphene1_no',int2str(i),'_0.0V.s2p'];
    ad9 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\5.4e-03\Graphene1_no',int2str(i),'_0.0V.s2p'];
    ad10 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\0.025\Graphene1_no',int2str(i),'_0.0V.s2p'];
    ad11 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\0.05\Graphene1_no',int2str(i),'_0.0V.s2p'];
    ad12 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\0.1\Graphene1_no',int2str(i),'_0.0V.s2p'];
    ad13 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\0.15\Graphene1_no',int2str(i),'_0.0V.s2p'];
    ad14= ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\0.2\Graphene1_no',int2str(i),'_0.0V.s2p'];
    ad15 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\0.25\Graphene1_no',int2str(i),'_0.0V.s2p'];
    ad16 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\0.3\Graphene1_no',int2str(i),'_0.0V.s2p'];
    ad17 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\0.4\Graphene1_no',int2str(i),'_0.0V.s2p'];
    ad18 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\0.5\Graphene1_no',int2str(i),'_0.0V.s2p'];
    ad19 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\0.75\Graphene1_no',int2str(i),'_0.0V.s2p'];
    ad20 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\1\Graphene1_no',int2str(i),'_0.0V.s2p'];
    ad21 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\1.5\Graphene1_no',int2str(i),'_0.0V.s2p'];
    ad22 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\2\Graphene1_no',int2str(i),'_0.0V.s2p'];
    ad23 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\4\Graphene1_no',int2str(i),'_0.0V.s2p'];
    ad24 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\6\Graphene1_no',int2str(i),'_0.0V.s2p'];
    ad25 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\8\Graphene1_no',int2str(i),'_0.0V.s2p'];
    ad26 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\10\Graphene1_no',int2str(i),'_0.0V.s2p'];
    ad27 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\12\Graphene1_no',int2str(i),'_0.0V.s2p'];
    ad28 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\project1\experimental-data\16\Graphene1_no',int2str(i),'_0.0V.s2p'];
    AD = {ad1,ad18,ad19,ad20,ad21,ad22,ad23,ad24,ad25,ad26,ad27,ad28};%,ad7,ad8,ad9,ad10,ad11,ad12,ad13,ad14,ad15,ad16,ad17,...
%          ad18,ad19,ad20,ad21,ad22,ad23,ad24,ad25,ad26,ad27,ad28};%,ad7,ad8,ad9,...
%         ad10,ad11,ad12,ad13,ad14,ad15,ad16,ad17,ad18,...
%         ad19,ad20,ad21,ad22,ad23,ad24,ad25,ad26,ad27,ad28};%
%     k = 1;m = 1;
    for j = 1:concen_num % concentrations
        s = sparameters(AD{1,j});
%         C{j,1} = s; 
        s11 = rfparam(s,1,1);
        s11_re = real(s11);
        s11_imag = imag(s11);
        s11_amp = abs(s11);
        s11_phase = angle(s11);
        s21 = rfparam(s,2,1);
        s21_re = real(s21);
        s21_imag = imag(s21);
        s21_amp = abs(s21);
        s21_phase = angle(s21);
        if i == 1
            s21_db = 20*log10(abs(s21));
            S21(j,:) = s21_db;
        end
        ss=cat(1,s11_imag,s11_re,s21_imag,s21_re);%,s21_re,s21_imag
        
%         if i<6 %i==2  %j==3 || j==10 || j==14%4 || j==8 || j==11 || j==18 || j==24 %i < 6
%             S_all_expe_clas(concen_num*(i-1)+j,:) = ss;%S_all_expe_clas(concen_num*(i-1)+j,:) = ss;S_test_expe(j,1:4004) = ss
% %             k = k+1;
% %         elseif i==5
% %             S_test_expe(j,4005:8008) = ss;
%         elseif i==6 
%             S_test_expe_clas(j,:) = ss;%S_test_expe_clas(j,:) = ss;S_all_expe(j,1:4004) = ss
% %             m = m+1;
% %         elseif i==8
% %             S_all_expe(j,4005:8008) = ss;
%         elseif i==8
%             S_test_expe_clas(j+28,:) = ss;
%         end

        if i==1
            S_all_expe(j,1+4004*(i-1):4004*i) = ss';
%         end
        elseif i==6
            S_test_expe(j,:) = ss';
        elseif i==8
            S_test_expe1(j,:) = ss';
        end                     
     end
end

%% feature extraction
% concen_num = 20;
% S_all_expe = S_all_expe(1:20,:);
Mean = zeros(concen_num*5,2);
Variance = zeros(concen_num*5,2);
PeaktoRMS = zeros(concen_num*5,2);
AmplitudeRange = zeros(concen_num*5,2);
STD = zeros(concen_num*5,2);
SS = [];
for i = 1: concen_num
    for j = 1:2
        if j==1
            for k = 1:5
%                 SS = [SS; S_all_expe_clas(i+28*(k-1),:)];%S_all_expe(i,:);
                S_current = S_all_expe_clas(i+28*(k-1),:);
                Mean(i+28*(k-1),j) = mean(S_current,'all');
                Variance(i+28*(k-1),j) = var(S_current,0,'all');
                PeaktoRMS(i+28*(k-1),j) = peak2rms(S_current,'all');
                AmplitudeRange(i+28*(k-1),j) = max(S_current,[],'all') - min(S_current,[],'all');
                STD(i+28*(k-1),j) = std(S_current,0,'all');
            end
        else
            for k = 1:2
            SS = S_test_expe_clas(i+28*(k-1),:);
            Mean(i+28*(k-1),j) = mean(SS,'all');
            Variance(i+28*(k-1),j) = var(SS,0,'all');
            PeaktoRMS(i+28*(k-1),j) = peak2rms(SS,'all');
            AmplitudeRange(i+28*(k-1),j) = max(SS,[],'all') - min(SS,[],'all');
            STD(i+28*(k-1),j) = std(SS,0,'all');
            end
        end
    end
end
for j = 1:2
    if j == 1
        Features_tra = [Mean(:,j) Variance(:,j) PeaktoRMS(:,j) AmplitudeRange(:,j) STD(:,j)];
    else 
        Features_test = [Mean(1:56,j) Variance(1:56,j) PeaktoRMS(1:56,j) AmplitudeRange(1:56,j) STD(1:56,j)];
    end
end
%%
% [S_all_newCole_norm,PS] = mapminmax(S_all_regr_newCole',-1,1);
S_all_newCole_norm = zscore(S_all_regr_newCole,0,1);
% s = svd(S_all_regr_newCole);
% s_norm = svd(S_all_newCole_norm);
% plot(s,'-o');
% hold on 
% plot(s_norm,'-*');
%% data for validation and testing
% S_test_Thru = zeros(7,4004);
% S_test_Thru_1 = zeros(7,4004);
% S_test_gra001_0 = zeros(7,12012);
% S_test_gra001_1 = zeros(7,12012);
S_train_newCole = zeros(32,4004);
S_test_newCole = zeros(32,4004);
S_test_newCole_1 = zeros(32,4004);
concen = linspace(0.5,16,32);
for k=1
    for j=1:32 %0:6 
    %     switch(i)
    %         case 0
    %             s_test = sparameters('C:\Users\93554\Desktop\Project 1\CST\Temp\waveguides\coplanar_tapered\SP_gra\test\SiO2+normal_scam2_7.7_26.s2p');
    %         case 1
    %             s_test = sparameters('C:\Users\93554\Desktop\Project 1\CST\Temp\waveguides\coplanar_tapered\SP_gra\test\strip_scam2_6.3_26.s2p');
    %         case 2
    %             s_test = sparameters('C:\Users\93554\Desktop\Project 1\CST\Temp\waveguides\coplanar_tapered\SP_gra\test\SiO2+normal_scam2_6.3.s2p');
    %         case 3
    %             s_test = sparameters('C:\Users\93554\Desktop\Project 1\CST\Temp\waveguides\coplanar_tapered\SP_gra\test\SiO2_scam2_6.3.s2p');
    %         case 4
    %             s_test = sparameters('C:\Users\93554\Desktop\Project 1\CST\Temp\waveguides\coplanar_tapered\SP_gra\test\water_scam2_6.3.s2p');
    %     end
%         s_test_Thru= sparameters(['C:\Users\93554\Desktop\Project 1\CST\Temp\waveguides\coplanar_tapered_copycopy\SP_comparision\Thru_airORglucose\scam_2_7\glucose_',int2str(i),'.s2p']);

%         switch(j) %addre for graphene, addre1 for Thru
%             case 0
%                 addre = 'C:\Users\93554\Desktop\Project 1\CST\Temp\waveguides\coplanar_tapered_copycopy\SPGlu\chemical-0.4\0.25\graphene0.01\scam_4_3-4_scam_1_4-6_';
%                 addre0 = 'C:\Users\93554\Desktop\Project 1\CST\Temp\waveguides\coplanar_tapered_copycopy\SPGlu\chemical-0.4\0.25\graphene0.01\scam_1_5-6_';
%                 addre1 = 'C:\Users\93554\Desktop\Project 1\CST\Temp\waveguides\coplanar_tapered_copycopy\SP_comparision\Thru_airORglucose\0.25\scam_4_3-4_scam_1_4-6_';
%                 addre11 = 'C:\Users\93554\Desktop\Project 1\CST\Temp\waveguides\coplanar_tapered_copycopy\SP_comparision\Thru_airORglucose\0.25\scam_1_4-5_scam_2_7.6-8_';
%                 addre2 = 'C:\Users\93554\Desktop\Project 1\CST\Temp\waveguides\coplanar_tapered_copycopy\SPGlu\chemical-0.4\100mg\scam_1_5-6_';
%             case 1
%                 addre = 'C:\Users\93554\Desktop\Project 1\CST\Temp\waveguides\coplanar_tapered_copycopy\SPGlu\chemical-0.4\0.5\graphene0.01\scam_4_3-4_scam_1_4-6_';
%                 addre0 = 'C:\Users\93554\Desktop\Project 1\CST\Temp\waveguides\coplanar_tapered_copycopy\SPGlu\chemical-0.4\0.5\graphene0.01\scam_1_5-6_';
%                 addre1 = 'C:\Users\93554\Desktop\Project 1\CST\Temp\waveguides\coplanar_tapered_copycopy\SP_comparision\Thru_airORglucose\0.5\scam_4_3-4_scam_1_4-6_';
%                 addre11 = 'C:\Users\93554\Desktop\Project 1\CST\Temp\waveguides\coplanar_tapered_copycopy\SP_comparision\Thru_airORglucose\0.5\scam_1_4-5_scam_2_7.6-8_';
%                 addre2 = 'C:\Users\93554\Desktop\Project 1\CST\Temp\waveguides\coplanar_tapered_copycopy\SPGlu\chemical-0.4\200mg\scam_1_5-6_';
%             case 2
%                 addre = 'C:\Users\93554\Desktop\Project 1\CST\Temp\waveguides\coplanar_tapered_copycopy\SPGlu\chemical-0.4\1\graphene0.01\scam_4_3-4_scam_1_4-6_';
%                 addre0 = 'C:\Users\93554\Desktop\Project 1\CST\Temp\waveguides\coplanar_tapered_copycopy\SPGlu\chemical-0.4\1\graphene0.01\scam_1_5-6_';
%                 addre1 = 'C:\Users\93554\Desktop\Project 1\CST\Temp\waveguides\coplanar_tapered_copycopy\SP_comparision\Thru_airORglucose\1\scam_4_3-4_scam_1_4-6_';
%                 addre11 = 'C:\Users\93554\Desktop\Project 1\CST\Temp\waveguides\coplanar_tapered_copycopy\SP_comparision\Thru_airORglucose\1\scam_1_4-5_scam_2_7.6-8_';
%                 addre2 = 'C:\Users\93554\Desktop\Project 1\CST\Temp\waveguides\coplanar_tapered_copycopy\SPGlu\chemical-0.4\250mg\scam_1_5-6_';
%             case 3
%                 addre = 'C:\Users\93554\Desktop\Project 1\CST\Temp\waveguides\coplanar_tapered_copycopy\SPGlu\chemical-0.4\2\graphene0.01\scam_4_3-4_scam_1_4-6_';
%                 addre0 = 'C:\Users\93554\Desktop\Project 1\CST\Temp\waveguides\coplanar_tapered_copycopy\SPGlu\chemical-0.4\2\graphene0.01\scam_1_5-6_';
%                 addre1 = 'C:\Users\93554\Desktop\Project 1\CST\Temp\waveguides\coplanar_tapered_copycopy\SP_comparision\Thru_airORglucose\2\scam_4_3-4_scam_1_4-6_';
%                 addre11 = 'C:\Users\93554\Desktop\Project 1\CST\Temp\waveguides\coplanar_tapered_copycopy\SP_comparision\Thru_airORglucose\2\scam_1_4-5_scam_2_7.6-8_';
%                 addre2 = 'C:\Users\93554\Desktop\Project 1\CST\Temp\waveguides\coplanar_tapered_copycopy\SPGlu\chemical-0.4\300mg\scam_1_5-6_';
%             case 4
%                 addre = 'C:\Users\93554\Desktop\Project 1\CST\Temp\waveguides\coplanar_tapered_copycopy\SPGlu\chemical-0.4\4\graphene0.01\scam_4_3-4_scam_1_4-6_';
%                 addre0 = 'C:\Users\93554\Desktop\Project 1\CST\Temp\waveguides\coplanar_tapered_copycopy\SPGlu\chemical-0.4\4\graphene0.01\scam_1_5-6_';
%                 addre1 = 'C:\Users\93554\Desktop\Project 1\CST\Temp\waveguides\coplanar_tapered_copycopy\SP_comparision\Thru_airORglucose\4\scam_4_3-4_scam_1_4-6_';
%                 addre11 = 'C:\Users\93554\Desktop\Project 1\CST\Temp\waveguides\coplanar_tapered_copycopy\SP_comparision\Thru_airORglucose\4\scam_1_4-5_scam_2_7.6-8_';
%                 addre2 = 'C:\Users\93554\Desktop\Project 1\CST\Temp\waveguides\coplanar_tapered_copycopy\SPGlu\chemical-0.4\400mg\scam_1_5-6_';
%             case 5
%                 addre = 'C:\Users\93554\Desktop\Project 1\CST\Temp\waveguides\coplanar_tapered_copycopy\SPGlu\chemical-0.4\8\graphene0.01\scam_4_3-4_scam_1_4-6_';
%                 addre0 = 'C:\Users\93554\Desktop\Project 1\CST\Temp\waveguides\coplanar_tapered_copycopy\SPGlu\chemical-0.4\8\graphene0.01\scam_1_5-6_';
%                 addre1 = 'C:\Users\93554\Desktop\Project 1\CST\Temp\waveguides\coplanar_tapered_copycopy\SP_comparision\Thru_airORglucose\8\scam_4_3-4_scam_1_4-6_';
%                 addre11 = 'C:\Users\93554\Desktop\Project 1\CST\Temp\waveguides\coplanar_tapered_copycopy\SP_comparision\Thru_airORglucose\8\scam_1_4-5_scam_2_7.6-8_';
%                 addre2 = 'C:\Users\93554\Desktop\Project 1\CST\Temp\waveguides\coplanar_tapered_copycopy\SPGlu\chemical-0.4\500mg\scam_1_5-6_';
%             case 6
%                 addre = 'C:\Users\93554\Desktop\Project 1\CST\Temp\waveguides\coplanar_tapered_copycopy\SPGlu\chemical-0.4\16\graphene0.01\scam_4_3-4_scam_1_4-6_';
%                 addre0 = 'C:\Users\93554\Desktop\Project 1\CST\Temp\waveguides\coplanar_tapered_copycopy\SPGlu\chemical-0.4\16\graphene0.01\scam_1_5-6_';
%                 addre1 = 'C:\Users\93554\Desktop\Project 1\CST\Temp\waveguides\coplanar_tapered_copycopy\SP_comparision\Thru_airORglucose\16\scam_4_3-4_scam_1_4-6_';
%                 addre11 = 'C:\Users\93554\Desktop\Project 1\CST\Temp\waveguides\coplanar_tapered_copycopy\SP_comparision\Thru_airORglucose\16\scam_1_4-5_scam_2_7.6-8_';
%         end

        if mod(concen(j),1) == 0
            addre = ['C:\Users\93554\OneDrive - University of Cambridge\projects\Project 1\CST\Temp\waveguides\coplanar_tapered_copycopy\SPGlu\chemical-0.4\',int2str(concen(j)),'gdl\scam_1_4-5_scam_2_7.6-8_',int2str(k),'.s2p'];
            addre2 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\Project 1\CST\Temp\waveguides\coplanar_tapered_copycopy\SPGlu\chemical-0.4\',int2str(concen(j)),'gdl\scam_4_3-4_scam_1_5-6_',int2str(k),'.s2p'];
            addre22 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\Project 1\CST\Temp\waveguides\coplanar_tapered_copycopy\SPGlu\chemical-0.4\',int2str(concen(j)),'gdl\scam_2_6.4-7.6_',int2str(k),'.s2p'];
        else
            addre = ['C:\Users\93554\OneDrive - University of Cambridge\projects\Project 1\CST\Temp\waveguides\coplanar_tapered_copycopy\SPGlu\chemical-0.4\',sprintf('%.1f',(concen(j))),'gdl\scam_1_4-5_scam_2_7.6-8_',int2str(k),'.s2p'];
            addre2 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\Project 1\CST\Temp\waveguides\coplanar_tapered_copycopy\SPGlu\chemical-0.4\',sprintf('%.1f',(concen(j))),'gdl\scam_4_3-4_scam_1_5-6_',int2str(k),'.s2p'];
            addre22 = ['C:\Users\93554\OneDrive - University of Cambridge\projects\Project 1\CST\Temp\waveguides\coplanar_tapered_copycopy\SPGlu\chemical-0.4\',sprintf('%.1f',(concen(j))),'gdl\scam_2_6.4-7.6_',int2str(k),'.s2p'];
        end

%         addre3 = ['C:\Users\93554\Desktop\Project 1\CST\Temp\waveguides\coplanar_tapered_copycopy\SPGlu\chemical-0.4\',int2str(j),'gdl\scam_1_4-5_scam_2_7.6-8_',int2str(k),'.s2p'];
%         s_test_gra001_1 = sparameters([addre,int2str(k),'.s2p']);
%         s11_t_gra = rfparam(s_test_gra001_1,1,1);
%         s11_t_gra_amp = abs(s11_t_gra);
%         s11_t_gra_phase = angle(s11_t_gra);
%         s21_t_gra = rfparam(s_test_gra001_1,2,1);
%         s21_t_gra_amp = abs(s21_t_gra);
%         s21_t_gra_phase = angle(s21_t_gra);
%         ss_t_gra=cat(1,s11_t_gra_amp,s11_t_gra_phase,s21_t_gra_amp,s21_t_gra_phase);
%         S_test_gra001_1(j+1,1+4004*(k-1):4004*k)= ss_t_gra';
%         
%         s_test_gra001 = sparameters([addre0,int2str(k),'.s2p']);
%         s11_t_gra_0 = rfparam(s_test_gra001,1,1);
%         s11_t_gra_amp_0 = abs(s11_t_gra_0);
%         s11_t_gra_phase_0 = angle(s11_t_gra_0);
%         s21_t_gra_0 = rfparam(s_test_gra001,2,1);
%         s21_t_gra_amp_0 = abs(s21_t_gra_0);
%         s21_t_gra_phase_0 = angle(s21_t_gra_0);
%         ss_t_gra_0=cat(1,s11_t_gra_amp_0,s11_t_gra_phase_0,s21_t_gra_amp_0,s21_t_gra_phase_0);
%         S_test_gra001_0(j+1,1+4004*(k-1):4004*k)= ss_t_gra_0';
        
%         s_test_Thru = sparameters([addre1,int2str(k),'.s2p']);
%         s11_t = rfparam(s_test_Thru,1,1);
%         s11_t_amp = imag(s11_t);
%         s11_t_phase = real(s11_t);
%         s21_t = rfparam(s_test_Thru,2,1);
%         s21_t_amp = imag(s21_t);
%         s21_t_phase = real(s21_t);
%         ss_t=cat(1,s11_t_amp,s11_t_phase,s21_t_amp,s21_t_phase);
%         S_test_Thru(j+1,1+4004*(k-1):4004*k)= ss_t';
%         
%         s_test_Thru_1 = sparameters([addre11,int2str(k),'.s2p']);
%         s11_t_1 = rfparam(s_test_Thru_1,1,1);
%         s11_t_amp_1 = imag(s11_t_1);
%         s11_t_phase_1 = real(s11_t_1);
%         s21_t_1 = rfparam(s_test_Thru_1,2,1);
%         s21_t_amp_1 = imag(s21_t_1);
%         s21_t_phase_1 = real(s21_t_1);
%         ss_t_1=cat(1,s11_t_amp_1,s11_t_phase_1,s21_t_amp_1,s21_t_phase_1);
%         S_test_Thru_1(j+1,1+4004*(k-1):4004*k)= ss_t_1';
%         
        s_train_newCole = sparameters(addre);
        s_test_newCole = sparameters(addre2);
        s_test_newCole_1 = sparameters(addre22);
        
        s11_r_newCole = rfparam(s_train_newCole,1,1);
        s11_r_amp_newCole = imag(s11_r_newCole);
        s11_r_phase_newCole = real(s11_r_newCole);
        s21_r_newCole = rfparam(s_train_newCole,2,1);
        s21_r_amp_newCole = imag(s21_r_newCole);
        s21_r_phase_newCole = real(s21_r_newCole);

        s11_t_newCole = rfparam(s_test_newCole,1,1);
        s11_t_amp_newCole = imag(s11_t_newCole);
        s11_t_phase_newCole = real(s11_t_newCole);
        s21_t_newCole = rfparam(s_test_newCole,2,1);
        s21_t_amp_newCole = imag(s21_t_newCole);
        s21_t_phase_newCole = real(s21_t_newCole);
        
        s11_t_newCole_1 = rfparam(s_test_newCole_1,1,1);
        s11_t_amp_newCole_1 = imag(s11_t_newCole_1);
        s11_t_phase_newCole_1 = real(s11_t_newCole_1);
        s21_t_newCole_1 = rfparam(s_test_newCole_1,2,1);
        s21_t_amp_newCole_1 = imag(s21_t_newCole_1);
        s21_t_phase_newCole_1 = real(s21_t_newCole_1);
        
        ss_r_newCole=cat(1,s11_r_amp_newCole,s11_r_phase_newCole,s21_r_amp_newCole,s21_r_phase_newCole);
        ss_t_newCole=cat(1,s11_t_amp_newCole,s11_t_phase_newCole,s21_t_amp_newCole,s21_t_phase_newCole);
        ss_t_newCole_1=cat(1,s11_t_amp_newCole_1,s11_t_phase_newCole_1,s21_t_amp_newCole_1,s21_t_phase_newCole_1);
    %     ss_t = [s21_t_imag(1); s21_t_re(1)];
        S_train_newCole(j,1+4004*(k-1):4004*k)= ss_r_newCole';
        S_test_newCole(j,1+4004*(k-1):4004*k)= ss_t_newCole';
        S_test_newCole_1(j,1+4004*(k-1):4004*k)= ss_t_newCole_1';
    end
end
% %normalisation
% S_test_newCole_norm = zscore(S_test_newCole,0,1);
% S_test_newCole_1_norm = zscore(S_test_newCole_1,0,1);

% S_test_newCole = S_test_newCole';
% S_test_newCole_1 = S_test_newCole_1';
% S_test_newCole_norm = mapminmax('apply',S_test_newCole,PS);
% S_test_newCole_1_norm = mapminmax('apply',S_test_newCole_1,PS);


% pred = ImagRe_gra001_24024_RLR_PCA2.predictFcn(S_test_gra001);
% pred = ImagRe_Thru_24024_RLR_PCA1.predictFcn(S_test_gra001);
% pred = gra001_24024_LR_PCA2.predictFcn(S_test_gra001);
% fit_Thru = Thru_LR_PCA1.predictFcn(S_test_Thru);
% fitt_Thru = LR_Thru_noPCA.predictFcn(S_test_Thru);
% fittt_Thru = LR_nonorm4004_Thru_noPCA_2k.predictFcn(S_test_Thru);

% yfit = Tree_S21imag_2festure.predictFcn(ss_t');
%% z-score_norm
% test = [1,2,3;4,5,6];
% Z = zscore(test,0,2);
% Z_S_all_regr_newCole_Thru = zscore(S_all_regr_newCole_Thru,0,1);
% S_all_expe0 = S_all_expe(1:7,:);
% Z_S_voltage = zscore(S_voltage,0,1);
% Z_S_voltage_test = zscore(S_voltage_test,0,1);

MU = zeros(1,4004);
SIGMA = zeros(1,4004);
% Z_S_test_newCole = zeros(32,4004);
% Z_S_test_newCole1 = zeros(32,4004);
% Z_S_train_newCole = zeros(32,4004);
Z_S_all_expe = zeros(concen_num,4004);%concen_num-3,732
Z_S_test_expe = zeros(concen_num,4004);
Z_S_test_expe1 = zeros(concen_num,4004);
Z_S_test_expe2 = zeros(concen_num,4004);
for i = 1:4004
    [Z_S_all_expe(:,i),MU(i),SIGMA(i)] = zscore(S_all_expe(:,i));%S_DC(1:6,:) ,0,'all'
% % Z_S_test_expe = zscore(S_test_expe,0,1);
%     Z_S_test_newCole(:,i) = (S_test_newCole(:,i) - MU(i))/SIGMA(i);
%     Z_S_test_newCole1(:,i) = (S_test_newCole_1(:,i) - MU(i))/SIGMA(i);
    Z_S_test_expe(:,i) = (S_test_expe(:,i) - MU(i))/SIGMA(i);
%     Z_S_test_expe2(:,i) = (S(:,i) - MU(i))/SIGMA(i);
     Z_S_test_expe1(:,i) = (S_test_expe1(:,i) - MU(i))/SIGMA(i);
end

% [Z_S_all_expe,mu,sigma] = zscore(S_all_expe_clas,0,'all');
% % Z_S_test_expe = zscore(S_test_expe,0,1);
% Z_S_test_expe = zeros(concen_num*2,4004);
% for i=1:concen_num%test_num
% Z_S_test_expe(i,:) = (S_test_expe_clas(i,:) - mu)/sigma;
% end

% % Z_S_all_expe_devi169 = zscore(S_all_expe_devi169,0,1);
% % Z_S_test_expe_devi169 = zscore(S_test_expe_devi169,0,1);
% % concen_num = 10;
% [Z_S_voltage,mu,sigma] = zscore(S_voltage,0,'all');%_expe_clas
% Z_S_voltage_test = zeros(concen_num,4004*18);%*5
% % S_test_expe0 = S_test_expe(1:concen_num,:);
% for j=1%:5
%     for i=1:concen_num
%     Z_S_voltage_test(i+concen_num*(j-1),:) = (S_voltage_test(i+concen_num*(j-1),:) - mu)/sigma;
%     end
% end
% Z_S_test_newCole = zscore(S_test_newCole,0,1);
% Z_S_test_newCole_1 = zscore(S_test_newCole_1,0,1);
%% dataset for linear regression
% S_test_Thru_1(:,24025) = S_target;
% S_test_Thru(:,24025) = S_target;
% S_test_Thru_1 = S_test_Thru_1(:,1:24024);
train_S_all_regr_newCole = [S_all_regr_newCole S_target_newCole];
%% testing for linear regression & classification
% fit_Thru = LR_Thru_predictor_v.predictFcn(S_test_Thru_1);
% clas_newCole = LDR_clas_192_4004.predictFcn(S_test_clas_newCole_Thru);

% clas_newCole_gra = LDR_clas_192_4004_PCA2.predictFcn(S_test_clas_newCole);
% clas_newCole_gra_1 = LDR_clas_192_4004_PCA2.predictFcn(S_test_clas_newCole_1);
% Thru_newCole = LR_Thruwater_12012_PCA2.predictFcn(S_test_newCole_Thru);
% Thru_onefreq = LR_Thru_32_24_PCA2.predictFcn(S_test_newCole2_1);
% Thru_LR = Plas_LR_12012_PCA3.predictFcn(predictors_v);
% Thru_LR_1 = Plas_LR_12012_PCA3.predictFcn(predictors_t);
% aaaa = LR_Water_12012_zscore_PCA4.predictFcn(predictors_v);%LR_newCole_4004imag_1norm;S_test_newCole_norm
% bbbb = LR_Water_12012_zscore_PCA4.predictFcn(predictors_t);

% fit_expe = RobustLinearRegression_15concentrations_5features.predictFcn(predictors_v);
% fit_expe = RobustLR_12concentrations_3features.predictFcn(predictors_v);
% prediction = RLR_32simulation_PCA5.predictFcn(Z_S_test_newCole1);%classification
prediction = Tree.predictFcn(Features_test);%LinearDiscriminant_10_4004_PCA3

% regression
% yfitPLS_expe = zeros(concen_num,5);
% yfitPLS_expe = KNN28_5extracted.predictFcn(Features_test);

% fit_newCole = LR_32_24_imagreal_PCA2.predictFcn(S_test_newCole2_1);
% fit_newCole_1 = LR_32_24_imagreal_PCA2.predictFcn(S_test_newCole2_2);%LR_05to16.predictFcn
%% new training & testing data
S_mix_train = [S_all_regr(:,1:12012) S_test_gra001_1(:,12013:24024)];
S_mix_test = [S_all_regr(:,12013:24024) S_test_gra001_1(:,1:12012)];
%% ANN & CNN data preparation
% S_target_newCole_water = [0;S_target_newCole];
% [XTrain,~,YTrain] = digitTrain4DArrayData;
% XTrain = S_all_regr;
% YTrain = S_target;
% XValidation = S_test_gra001_1;
% XTesting = S_test_gra001_0;
samp_num = 12;%32
% convolution_num = 3969;%11881 for 1-3
% convolution_num_sqrt = 63;%109 for 1-3
% % XTrain = S_full_regr(:,1:(size(S_full_regr,2)-1));
XTrain = Z_S_all_expe;%S_all_newCole_norm'；Z_S_all_regr_newCole_Thru;S_all_newCole_norm
YTrain = S_target_expe;%_newCole
% XTesting = S_test_gra001_0_norm;
XValidation = Z_S_test_expe;
% XValidation = S_test_newCole_norm;%S_test_gra001_0;S_test_newCole_norm'
% XTesting = S_test_newCole_1_norm;%S_all_regr_newCole_Thru_t;S_test_newCole_Thru_1;S_test_newCole_1_norm'


% XTrain = S_train_newCole2';
% YTrain = S_target_newCole;
% XValidation = S_test_newCole2_1';
% XTesting = S_test_newCole2_2';

% PCA
[pcaCoefficients, pcaScores, ~, ~, explained, pcaCenters] = pca(...
    XTrain);
[pcaCoefficients_v, pcaScores_v, ~, ~, explained_v, pcaCenters_v] = pca(...
    XValidation);
% [pcaCoefficients_t, pcaScores_t, ~, ~, explained_t, pcaCenters_t] = pca(...
%     XTesting);
% Keep enough components to explain the desired amount of variance.
explainedVarianceToKeepAsFraction = 95/100;
numComponentsToKeep = find(cumsum(explained)/sum(explained) >= explainedVarianceToKeepAsFraction, 1);
numComponentsToKeep = 3;
pcaCoefficients = pcaCoefficients(:,1:numComponentsToKeep);
predictors = pcaScores(:,1:numComponentsToKeep);
% array2table(), predictors(:, isCategoricalPredictor)

numComponentsToKeep_v = find(cumsum(explained_v)/sum(explained_v) >= explainedVarianceToKeepAsFraction, 1);
pcaCoefficients_v = pcaCoefficients_v(:,(1:(numComponentsToKeep)));
predictors_v = pcaScores_v(:,1:(numComponentsToKeep));
% % 
% numComponentsToKeep_t = find(cumsum(explained_t)/sum(explained_t) >= explainedVarianceToKeepAsFraction, 1);
% pcaCoefficients_t = pcaCoefficients_t(:,1:(numComponentsToKeep_t-1));
% predictors_t = pcaScores_t(:,1:(numComponentsToKeep_t-1));

% % ads_t = arrayDatastore(XTrain(:,501:700)','ReadSize',200,'IterationDimension',2, ...
% %  'OutputType','cell');
% % test = read(ads_t);

% % original size for ANN (data 32*24)
% ads0_ori = arrayDatastore(XTrain,'ReadSize',samp_num,'IterationDimension',2, ...
%  'OutputType','cell');
% adsAngles0_ori = arrayDatastore(YTrain,'ReadSize',samp_num,'OutputType','cell');
% t0x_ori = read(ads0_ori);% to see whether data is in the right dimension
% t0y_ori = read(adsAngles0_ori);
% ads_ori = arrayDatastore(XValidation,'ReadSize',samp_num,'IterationDimension',2, ...
%  'OutputType','cell');
% adss_ori = arrayDatastore(XTesting,'ReadSize',samp_num,'IterationDimension',2, ...
%  'OutputType','cell');
% cds0_ori = combine(ads0_ori,adsAngles0_ori);
% cds_ori = combine(ads_ori,adsAngles0_ori);
% cdss_ori = combine(adss_ori, adsAngles0_ori);

% %original size for CNN
% xTrain = zeros(convolution_num_sqrt,convolution_num_sqrt,samp_num);
% xValidation = zeros(convolution_num_sqrt,convolution_num_sqrt,samp_num);
% xTesting = zeros(convolution_num_sqrt,convolution_num_sqrt,samp_num);
% for i = 1:samp_num
%     xtrain = XTrain(i,:);
%     xvalidation = XValidation(i,:);
%     xtesting = XTesting(i,:);
%     xtrain1 = xtrain(1:convolution_num);
%     xvalidation1 = xvalidation(1:convolution_num);
%     xtesting1 = xtesting(1:convolution_num);
%     xTrain(:,:,i) = reshape(xtrain1,convolution_num_sqrt,convolution_num_sqrt);
%     xValidation(:,:,i) = reshape(xvalidation1,convolution_num_sqrt,convolution_num_sqrt);
%     xTesting(:,:,i) = reshape(xtesting1,convolution_num_sqrt,convolution_num_sqrt);
% end
% ads0_ori = arrayDatastore(xTrain,'ReadSize',samp_num,'IterationDimension',3, ...
%  'OutputType','cell');
% adsAngles0_ori = arrayDatastore(YTrain,'ReadSize',samp_num,'OutputType','cell');
% t0x_ori = read(ads0_ori);% to see whether data is in the right dimension
% ads_ori = arrayDatastore(xValidation,'ReadSize',samp_num,'IterationDimension',3, ...
%  'OutputType','cell');
% adss_ori = arrayDatastore(xTesting,'ReadSize',samp_num,'IterationDimension',3, ...
%  'OutputType','cell');
% cds0_ori = combine(ads0_ori,adsAngles0_ori);
% cds_ori = combine(ads_ori,adsAngles0_ori);
% % cdss_ori = combine(adss_ori, adsAngles0_ori);
% % 
%PCA for ANN
ads0 = arrayDatastore(predictors','ReadSize',samp_num,'IterationDimension',2, ...
 'OutputType','cell');
adsAngles0 = arrayDatastore(YTrain,'ReadSize',samp_num,'OutputType','cell');

% ads = arrayDatastore(predictors_v','ReadSize',samp_num,'IterationDimension',2, ...
%  'OutputType','cell');

% adss = arrayDatastore(predictors_t','ReadSize',size(predictors_t,1),'IterationDimension',2, ...
%  'OutputType','cell');
% vx = read(ads0);
% vy = read(adsAngles0);
% 
% cds0 = combine(ads0,adsAngles0);%Training
% cds = combine(ads,adsAngles0);%Validation
% cdss = combine(adss, adsAngles0);%Testing
% read(cds);
% read(cdss);

% cat_S_target_clas_newCole = categorical(S_target_clas_newCole);
% ads0_clas = arrayDatastore(S_all_clas_newCole_Thru','ReadSize',192,'IterationDimension',2, ...
%  'OutputType','cell');
% adsAngles0_clas = arrayDatastore(cat_S_target_clas_newCole,'ReadSize',192,'OutputType','cell');
% 
% ads_clas = arrayDatastore(S_test_clas_newCole_Thru','ReadSize',192,'IterationDimension',2, ...
%  'OutputType','cell');

% cds0_clas = combine(ads0_clas,adsAngles0_clas);%Training
% cds_clas = combine(ads_clas,adsAngles0_clas);%Validation
% 
% vxclas = read(ads0_clas);
% vyclas = read(adsAngles0_clas);
%% neural net fitting data preparation
X_train_vali = [predictors;predictors_v];
Y_train_vali = [S_target;S_target];
%% ANN & CNN test
% % [Y,Xf,Af] = ANN(S_test_gra001);
% % [YY,~,~] = ANN5t2v(S_test_gra001);
% % [YYY,~,~] = PCA3_ANN10_6t1v(predictors_v);
% % [YYYY,~,~] = PCA3_ANN10_7t(predictors_v);
% % vx1 = subset(adss,1);
% % tx1 = XTrain(:,501:700);
% % Ynet_PCA_ANN310sigmid1_new = predict(ANN310sigmid1,predictors_v);
% vANN_result = predict(trainedNetwork_6,predictors_v);
% % tANN_result = predict(trainedNetwork_3,predictors_t);
% ANN_result = predict(trainedNetwork_6,predictors);

% tCNN_result = zeros(27,1);
% % cd 'C:\Users\93554\OneDrive - University of Cambridge\projects\project1\Images-test\0'
% % vCNN_result = zeros(32,1);
% for i=1:27
%     image_t = S_image_test(:,:,i);
%     image_t = uint8(255 * mat2gray(image_t));
% %     imwrite(image,[int2str(i),'.png']);
%     tCNN_result(i) = classify(trainedNetwork_1,image_t);
% %     vCNN_result(i) = predict(trainedNetwork_2,xValidation(:,:,i));
% end

image_t = S_image_test(:,:,27);
    image_t = uint8(255 * mat2gray(image_t));classify(trainedNetwork_1_7samples,image_t);
%% PLSR
[n,p] = size(Z_S_all_expe);%S_all_regr_newCole;S_all_newCole_norm;Z_S_all_expe;Z_S_voltage;Z_S_all_expe_devi169
y = S_target_expe;%S_voltage_Target;S_target_expe
% y = concen';
% y = S_target_nor';
% y = S_full_regr(:,size(S_full_regr,2));
%PLS components
k = 4;%27试一下小点值会不会好些 - 貌似并没有
[Xloadings,Yloadings,Xscores,Yscores,betaPLS10,PLSPctVar,MSE] = plsregress(...
	Z_S_all_expe,y,k);%Z_S_voltage;11
RMSE_tra = sqrt(MSE);
% hold on; for i=1:2 plot(Xloadings(:,i),'LineWidth',1); end %hold off
% plot(betaPLS10,'LineWidth',2);
%% PLSR vs PCR
[Xl,Yl,Xs,Ys,beta,pctVar,PLSmsep] = plsregress(Z_S_all_expe,S_target_expe,11,'CV',3);
PCRmsep = sum(crossval(@pcrssee,Z_S_all_expe,S_target_expe,'KFold',3),1) / concen_num;
plot(0:11,PLSmsep(2,:),'b-o',0:11,PCRmsep,'r-^');
xlabel('Number of components');
ylabel('Estimated Mean Squared Prediction Error');
legend({'PLSR' 'PCR'},'location','NE');
%% PLSR testing
% conclusion: target normalise/not normalise didn't affect the final regression result
% with features and symbols normalisation:
% yfitPLS_1 = mapminmax('reverse',yfitPLS_1_nor,ts);
% S_test_gra001_0_nor = mapminmax('apply',S_test_gra001_0,PS);
% yfitPLS_expe = [ones(n,1) Z_S_test_expe]*betaPLS10;%Z_S_test_expe;Z_S_voltage_test;(1:concen_num,:)

% yfitPLS_expe = zeros(concen_num,5);
% for i = 1%:5
    yfitPLS_expe = [ones(n,1) Z_S_test_expe]*betaPLS10;%concen_num
% end
%% calculate RMSE
% RMSE = zeros(5,1);
% R_squared1 = zeros(5,1);
% R_squared = zeros(5,1);
for i=1%:5
    ygra = yfitPLS_expe(:,i);%*sigma+mu
%    ygra = yhat;
    RMSE = sqrt(mean((S_target_expe - ygra).^2));  % Root Mean Squared Error; S_target_expe;S_voltage_Target
    % RMSE_new = sqrt(mean((S_all_regr_newCole_target - ygra).^2));
%     MSE = mean((S_target_expe - ygra).^2);
%     Var = mean((S_target_expe - mean(S_target_expe)).^2);
    % R_squared = 1 - MSE/Var;

    % adjusted R-square
    
    R_squared0 = corrcoef(S_target_expe,ygra);
    R_squared1 = R_squared0(2)^2;
    R_squared2 = 1-(n-1)/(n-k)*(1-R_squared1);
end
% result = yfitPLS';
% ANN_result = ANN3_10_20_1_result';
%%----------------------------------------------------------------------------------
%% k-fold validaiton
max_compo = 8;
rMSE = zeros(1,concen_num-test_num);
r_squared1 = zeros(1,concen_num-test_num);
RMSE = zeros(1,max_compo);
R_squared1 = zeros(1,max_compo);
for k = 1:max_compo
    for i = 1: concen_num-test_num
        Sk = Z_S_all_expe;
        Sk(i,:) = [];
        [n,p] = size(Sk);
        y = S_train_target;
        y(i,:) = [];
    [Xloadings,Yloadings,Xscores,Yscores,betaPLS10,PLSPctVar,MSE] = plsregress(...
        Sk,y,k);
    
        yfitPLS_expe = [1 Z_S_all_expe(i,:)]*betaPLS10;
        rMSE(i) = (S_train_target(i) - yfitPLS_expe).^2;
%         r_squared0 = corrcoef(S_target_expe(i),yfitPLS_expe);
%         r_squared1(i) = r_squared0^2;      
    end
    RMSE(k) = sum(rMSE)/(concen_num-test_num);
%     R_squared1(k) = mean(r_squared1);
end
plot(linspace(1,max_compo,max_compo),RMSE,'-o','LineWidth',2,'MarkerSize',5);
xlabel('PLS component number')
ylabel('Estimated MSE of leave-one-out cross-validation');
ax = gca; % current axes
ax.FontSize = 14;
%% Ridge regression
% Z_S_target_expe = zscore(S_target_expe,0,'all');
k = 0.2;
b = ridge(S_target_expe,Z_S_all_expe,k,0);
yhat = b(1)+Z_S_test_expe*b(2:end);%b(1) + 
%% predict one
% SS_test = zeros(1,100);
S_test = zeros(1,4004);
for i=1
    s_test = sparameters(['C:\Users\93554\Desktop\Project 1\CST\Temp\waveguides\coplanar_tapered_copycopy\SPGlu\chemical-0.6\16\scam_2_7\glucose16_gap_l_0-5_',int2str(i),'.s2p']);
    s11_t = rfparam(s_test,1,1);
    s11_t_re = real(s11_t);
    s11_t_imag = imag(s11_t);
%     s11_t_amp = abs(s11_t);
%     s11_t_phase = angle(s11_t);
    s21_t = rfparam(s_test,2,1);
    s21_t_re = real(s21_t);
    s21_t_imag = imag(s21_t);
%     s21_t_amp = abs(s21_t);
%     s21_t_phase = angle(s21_t);
    ss_t=cat(1,s11_t_imag,s11_t_re,s21_t_imag,s21_t_re);
%     ss_t_norm = mapminmax('apply',ss_t,PS);
%     [ss_t_norm, PPS] = mapminmax(ss_t,-1,1);%??为啥这里这个归一化函数不好使？
    S_test(1+4004*(i-1):4004*i) = ss_t;
end
S_test_norm = (S_test - min(S_test,[],'all'))/(max(S_test,[],'all')- min(S_test,[],'all'));%这里归一化用到了测试集合的信息，是不正确的
S_test_norm_column = zeros(1,4004);
for m = 1:4004
    S_test_norm_column(1,m) = (S_test(1,m) - min_eachcolumn(1,m))/(max_eachcolumn(1,m)- min_eachcolumn(1,m));
end
% for m = 1:4004
%     S_test_norm_column(1,m) = (S_test(1,m) - mean(S_test))/std(S_test);
% end
% RandIndex = randperm(1001);%只用到S11imag
% for i = 1:10
%     SS_test(:,1+25*(i-1):25*i) = S_test(:,4004*(i-1)+RandIndex(1:25));
% end

% SS_test = S_test(3504:3603);

% SS_test = zeros(1,400);
% for j = 1:4
%     SS_test(:,1+100*(j-1):100*j) = S_test_norm(:,4004+RandIndex(1:100)+1001*(j-1));%+4004-change device...too much influence?
% end
% SS_test = S_test_norm(3504:3603);
% f = Test.predictFcn(S_test_norm);

% fit = LR_4004Srownorm_PCA1.predictFcn(S_test_norm);
% fitt = LR_4004columnNorm_PCA3.predictFcn(S_test_norm);
% fittt = LR_4004_nonorm_PCA1.predictFcn(S_test);
% fitttt = LR_4004Smeanstd_PCA3.predictFcn(S_test_norm);
f = LR_nonorm4004_noPCA.predictFcn(S_test);
% y_fit = LR_PCA2_All_imag_re.predictFcn(S_test_norm);
% yy_fit = LR_PCA2_AllSnorm_imag_re.predictFcn(S_test_norm);
% yyy_fit = LR_PCA2_AllS_imag_re_trainedwithScam_2_7.predictFcn(S_test_norm);
% % y_fit = LR_PCA1_scam2_7_500600S21re.predictFcn(SS_test);
% % y_fit = LR_400AllS_imag_re1.predictFcn(SS_test);
% % y_fit = LR_allS_imag_re2.predictFcn(SS_test);
% yyyy_fit = Linear_PCA2_RMSE4.predictFcn(S_test_norm);
% yyyyy_fit = LR_AllSnorm_PCA1_35043603.predictFcn(SS_test);
%% adding symbols and normalisation of training
mater_num = 7;
sam_num = 3;
feature_num = 4004*sam_num;
S_full_regr = zeros(mater_num,feature_num+1);
% S_full = zeros(mater_num*sam_num,feature_num+1);
[S_all_regr_norm, PS] = mapminmax(S_all_regr,0,1);
S_full_regr(:,1:feature_num)=S_all_regr_norm;
% min_eachrow = min(S_all_regr,[],2);
% max_eachrow = max(S_all_regr,[],2);
% min_eachcolumn = min(S_all_regr,[],1);
% max_eachcolumn = max(S_all_regr,[],1);
% S_all_regr_norm = zeros(7,4004);
% for n =1:7
%     S_all_regr_norm(n,:) = (S_all_regr(n,:) - min_eachrow(n,1))/(max_eachrow(n,1)- min_eachrow(n,1));
% end
% for n =1:4004
%     S_all_regr_norm(:,n) = (S_all_regr(:,n) - min_eachcolumn(1,n))/(max_eachcolumn(1,n)- min_eachcolumn(1,n));
% end

% mu = mean(S_all_regr);
% sigma = std(S_all_regr);

% for n =1:4004%column_meanstd_norm
%     S_all_regr_meanstd(:,n) = (S_all_regr(:,n) - mu(n))/sigma(n);
% end
% S_all_regr_norm = (S_all_regr - min(S_all_regr,[],'all'))/(max(S_all_regr,[],'all')- min(S_all_regr,[],'all'));

% S_full_regr_Thru = zeros(mater_num,feature_num+1);
% S_full_regr_Thru(:,1:feature_num) = S_all_regr_Thru;

% S_full_regr_norm=zeros(mater_num,feature_num+1);
% S_full_regr_norm(:,1:feature_num)=S_all_regr_norm;
y = zeros(sam_num,1);
for i=1:mater_num
    switch(i)
        case 1
            yy=0.25;
        case 2
            yy=0.5;
        case 3
            yy=1;
        case 4
            yy=2;
        case 5
            yy=4;
        case 6
            yy=8;
        case 7
            yy=16;
    end
%     for j=1:sam_num
%         y(j) = yy;
%     end
%     S_full((sam_num*(i-1)+1):sam_num*i,feature_num+1)= y;
%     S_full_regr_Thru(i,size(S_full_regr_Thru,2)) = yy;
    S_full_regr(i,size(S_full_regr,2)) = yy;
%     S_full_regr_norm(i,size(S_full_regr_norm,2)) = yy;
end
S_target = S_full_regr(:,size(S_full_regr,2));
[S_target_nor,ts] = mapminmax(S_target',0,1);
% SS = zeros(6,251);

% SS = zeros(7,101);
% SS(:,1:100) = S_full_regr_norm(:,3504:3603);
% SS(:,101) = S_full_regr(:,size(S_full_regr,2));

% SS = zeros(7,401);
% RandIndex = randperm(1001);
% for j = 1:4
%     SS(:,1+100*(j-1):100*j) = S_full_regr_norm(:,RandIndex(1:100)+1001*(j-1));
% end
% SS(:,401) = S_full_regr_norm(:,size(S_full_regr,2));

% RandIndex = randperm(size(S_full_regr,2));
% RandIndex = randperm(1001);%只用到S11imag
% for i = 1:10
%     SS(:,1+25*(i-1):25*i) = S_full_regr(:,4004*(i-1)+RandIndex(1:25));
% end
% SS(:,251) = S_full_regr(:,size(S_full_regr,2));

% RandIndex = randperm(size(S_full,1));
% SS_full = S_full(RandIndex,:);
%% new cole-cole data functions
all_concentration = linspace(1000,16000,16);
Y_inf = zeros(1,16);
Y_0 = zeros(1,16);
tau = zeros(1,16);
for i = 1:16
    Y_inf(i) = Cole(-8.214*10^(-8),2.148*10^(-3),8.722,all_concentration(i));
    Y_0(i) = Cole(2.318*10^(-9),-2.793*10^(-4),81.015,all_concentration(i));
    tau(i) = Cole(-8.37*10^(-9),5.15*10^(-4),8.776,all_concentration(i));
end
%% correlation between concentration and s parameters
COR = cell(1,4);
COR_1 = cell(1,4);
target = S_target';
for i =1:4
    COR_1{1,i} = corrcoef(target,s_fixfre(i,:));
end

%% small data set
y0 = zeros(49,1);
y1 = ones(49,1);
y2 = zeros(49,1);
for i=1:49
    y2(i)=2;
end
S_0 = [S_0 y0];
S_1 = [S_1 y1];
S_2 = [S_2 y2];
S = [S_0; S_1; S_2];
RandIndex = randperm(size(S,1));
SS = S(RandIndex,:);
%% test
load ionosphere
ionosphere = array2table(X);
ionosphere.Group = Y;
%%
%import data
opts = delimitedTextImportOptions();
opts.Whitespace = '''';
opts.LeadingDelimitersRule = 'ignore';
opts.ConsecutiveDelimitersRule = 'join';
opts.DataLines = 3;
opts.VariableNames = {'Frequency / GHz','S1,1 (w=49)/abs,dB'};
opts.Delimiter = ' ';
% M = readtable('C:\Users\93554\Desktop\Project 1\CST\Temp\waveguides\coplanar\S11_0.txt',opts);
%N = readtable('C:\Users\93554\Desktop\Project 1\CST\Temp\waveguides\coplanar\S12_0.txt',opts);
% f = str2double(M.Frequency_GHz);
% s = str2double(M.S11);
f_ = str2double(N.Frequency_GHz);
s_ = str2double(N.S21);
plot(f_,s_,'r');
ylabel('S21 / dB');
xlabel('frequency / GHz');
% s = s.';
% s = abs(s);
% [load,scores,var] = pca(s,'Centered',false);
%ss = scores*load(1:10);
% s_origin = s.';
% plot(1:1001, s_origin,'bo');
% hold on;
%plot(1:10,ss,'r^');

%preview('airlinesmall_subset.xlsx',opts)
% x3=dlmread('C:\Users\93554\Desktop\Project 1\CST\Temp\waveguides\coplanar\S11_0.txt',',',2) ;
% [x1,x2] = textscan('C:\Users\93554\Desktop\Project 1\CST\Temp\waveguides\coplanar\S11_0.txt','%f%f',1001);
% [PCALoadings,PCAScores,PCAVar] = pca(X,'Economy',false);
% betaPCR = regress(y-mean(y), PCAScores(:,1:10));
% %To make the PCR results easier to interpret in terms of the original 
% %spectral data, transform to regression coefficients for the original, uncentered variables.
% betaPCR = PCALoadings(:,1:10)*betaPCR;
% betaPCR = [mean(y) - mean(X)*betaPCR; betaPCR];
% yfitPCR = [ones(n,1) X]*betaPCR;
% %plot
% plot(y,yfitPCR,'r^');
% xlabel('Frequency / GHz');
% ylabel('Fitted Response');
% legend({'PLSR with 10 Components' 'PCR with 10 Components'},  ...
% 	'location','NW');
