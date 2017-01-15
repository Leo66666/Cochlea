%Lei Guo Nov/23/2016
%This file is to culculate the numerical solution of F_vis and F_p
%with changing of frequency and G and find the relationship of
%F_vis with G & frequency and F_vis with G & frequency.
%Then compare the solution of numerical method and theoretical method

clc, clear, close all;
[num, txt, raw] = xlsread('comsol_data.xlsx');
data_raw = str2double(raw);
%change the size if change domain in comsol!
frequency_top = 2500; %Maximum of the frequency
frequency_bot = 100; %Minimum of the frequency
frequency_step = 600; %Frequency step
G_top = 6e-6; %Maximum of G
G_bot = 0; %Minimum of G
G_step = 0.2e-6; %G step

frequency_size = (frequency_top-frequency_bot)/frequency_step+1; %Number of the frequency
G_size = floor((G_top-G_bot)/G_step)+1; %Number of G
P_positive_raw = zeros(G_size , frequency_size); %Raw data of positive P
P_negative_raw = zeros(G_size , frequency_size); %Raw data of negative P
F_vis_raw = zeros(G_size , frequency_size); %Raw data of viscous force
dpdx = zeros(G_size , frequency_size); %Raw data of dpdx

%Relocate the data and culculate the simulation results
%For G(i) and Frequency(j), save positive P, negative P, viscous force and dpdx data to (i , j)
for i = 1 : G_size
    for j = 1 : frequency_size
        P_positive_raw(i , j) = data_raw((i-1)*frequency_size+j , 1);
        P_negative_raw(i , j) = abs(data_raw((i-1)*frequency_size+j , 2));
        F_vis_raw(i , j) = data_raw((i-1)*frequency_size+j , 3);
        dpdx(i , j) = data_raw((i-1)*frequency_size+j , 4);
        F_pre(i , j) = abs(P_positive_raw(i , j)-P_negative_raw(i , j)); %F_pre is the simulation result of pressure force
        F_vis(i , j) = abs(F_vis_raw(i , j)); %F_vis is the simulation result of viscous force
%       F_total(i , j) = abs(F_vis_raw(i , j)+P_positive_raw(i , j)-P_negative_raw(i , j)); %F_total is the simulation result of the total force
    end
end

figure
subplot(1,2,1)
plot([G_bot : G_step : G_top] , F_pre(: , 3))
title('Comsol Pressure torque')
subplot(1,2,2)
plot([G_bot : G_step : G_top] , F_vis(: , 3))
title('Comsol Viscous torque')


%This section is going to culculate the theoretical solution.
certain = 11;
for k = 1 : G_size

    %log in all the basic geometry variables which are same as comsol.
    G0 = 6e-6;
    G = G_bot+G_step*(k-1);
    H = 5.5e-6;
    alpha = (G0-G)/H;
    L = 50e-6;
    R = 0.2e-6;
    a = 1.3e-6;
    b = 1.5e-6;
    c = 1e-6;
    d = 0;%0.2e-6;
    utms = 0.05e-6;
    
    density = 1000;
    
    
    x_span = [-R : 0.01*R : R];
    %calculate g(x)
    g = G0-H-sqrt(R^2-x_span.^2)-d/2+d/2/R.*x_span-d/2/R*utms;
    %calculate Ii
    I1 = trapz(x_span,1./g);
    I2 = trapz(x_span,1./g.^2);
    I3 = trapz(x_span,1./g.^3);
    %calculate g0(y),ksi,ksi0,omega
    for index = 1:101
        ymin = 0.01*(1-alpha)*H*(index-1)+alpha*H;
        y = [ymin:0.01*(H-ymin):H];
        g0 = -(y-G0+d)*a/(G-d)-utms;
        ksi(index) = trapz(y,1./g0.^3);
    end
    ksi0 = ksi(1);
    y = [alpha*H:0.01*(1-alpha)*H:H];
    for index = 1:101
        htemp(index) = y(index)*ksi(index);
    end
    omega = trapz(y,htemp/H);
    
    miu=1.01*10^(-3);
    i = sqrt(-1);
    
        for j = 1 : frequency_size
            % This index(j) is for frequency.
            frequency(j) = frequency_bot+frequency_step*(j-1);
            gama(j) = sqrt(-i*frequency(j)*2*pi*density/miu);
            A(j) = 6*10^-6-2*(cosh(gama(j)*6*10^-6)-1)/(gama(j)*sinh(gama(j)*6*10^-6));
            B(j) = (cosh(gama(j)*6*10^-6)-1)/(gama(j)*sinh(gama(j)*6*10^-6));
            y = [alpha*H:0.01*(1-alpha)*H:H];
            for index = 1 : 101
                delta_P_HB(index) = -12*miu*A(j)/(i*2*pi*frequency(j)*density)*(I3+ksi(index))*dpdx(k,j)+(12*miu*i*2*pi*frequency(j)*B(j)*(I3+ksi(index))-6*miu*i*2*pi*frequency(j)*I2)*utms;
                dtemp(index) = y(index)*delta_P_HB(index);
            end
            inte = trapz(y,dtemp);
            tp(j) = -(alpha^2*H^2*delta_P_HB(1)/2+inte);
            tvis(j) = 6*miu*I2*H*A(j)/(i*2*pi*frequency(j)*density)*dpdx(k,j)+2*i*2*pi*frequency(j)*miu*H*(I1-3*I2*B(j))*utms;
            Q2(j) = 1/(i*2*pi*frequency(j)*density)*dpdx(k,j)*A(j)-i*2*pi*frequency(j)*B(j)*utms;
            %for certain frequency(certain) and G(k), save data.
            F_pre_theo(k,j) = abs(tp(j));
            F_vis_theo(k,j) = abs(tvis(j));
            Q2_theory(k,j) = (Q2(j));
            Q2_p(k,j)=1/(i*2*pi*frequency(j)*density)*dpdx(k,j)*A(j);
            Q2_v(k,j)=-i*2*pi*frequency(j)*B(j)*utms;
            A1(k,j)=A(j);
            B1(k,j)=B(j);
        end
end

figure
subplot(1,2,1)
plot([G_bot : G_step : G_top] , F_pre_theo(: , 3))
title('Theoretical Pressure torque')
subplot(1,2,2)
plot([G_bot : G_step : G_top] , F_vis_theo(: , 3))
title('Theoretical Viscous torque')

figure
plot([G_bot : G_step : G_top] , abs(Q2_v(: , 3)))
% 
% figure
% plot([G_bot : G_step : G_top] , F_certain_notunderline)
% 
% %Plot the F_total vs G at frequency at 2000 (frequency can be changed by
% %changing i)
% figure
% subplot(2,2,1)
% for i=certain:certain
%     G_interp = linspace(G_bot , G_top);
%     F_interp = interp1([G_bot : G_step : G_top] , F_total(: , i) , G_interp , 'pchip');
%     plot(G_interp , F_interp)
%     hold on
%     %plot([G_bot : G_step : G_top] , F_certain_notunderline,'r');
%     hold on
%     %plot([G_bot : G_step : G_top] , F_certain_underline,'k');
% end
% %Plot the F_pre vs G at frequency at 2000 (frequency can be changed by
% %changing i)
% subplot(2,2,2)
% for i=certain:certain
%     G_interp = linspace(G_bot , G_top);
%     F_interp = interp1([G_bot : G_step : G_top] , F_pre(: , i) , G_interp , 'pchip');
%     plot(G_interp , F_interp)
%     hold on
%     %plot([G_bot : G_step : G_top] , Fpre_certain_notunderline,'r');
%     hold on
%     %plot([G_bot : G_step : G_top] , Fpre_certain_underline,'k');
% end
% %Plot the F_pre vs G at frequency at 2000 (frequency can be changed by
% %changing i)
% subplot(2,2,3)
% for i=certain:certain
%     G_interp = linspace(G_bot , G_top);
%     F_interp = interp1([G_bot : G_step : G_top] , F_vis(: , i) , G_interp , 'pchip');
%     plot(G_interp , F_interp)
%     hold on
%     %plot([G_bot : G_step : G_top] , Fvis_certain_notunderline,'r');
%     hold on
%     %plot([G_bot : G_step : G_top] , Fvis_certain_underline,'k');
% end

%             Is_underlined = 1;
%             Fpre_underline(j) = 6*miu*A(j)/(i*2*pi*frequency(j)*density)*(I2+I3*H+Is_underlined*(ksi0*alpha^2*H+2*omega))*dpdx(k,j);
%             Fvis_underline(j) = i*2*pi*frequency(j)*miu*(2*I1+3*I2*H-B(j)*(6*I2+6*I3*H+Is_underlined*(6*ksi0*alpha^2*H+12*omega)))*utms;
%             F_underline(j) = Fpre_underline(j)+Fvis_underline(j);
%             Is_underlined = 0;
%             Fpre_notunderline(j) = 6*miu*A(j)/(i*2*pi*frequency(j)*density)*(I2+I3*H+Is_underlined*(ksi0*alpha^2*H+2*omega))*dpdx(k,j);
%             Fvis_notunderline(j) = i*2*pi*frequency(j)*miu*(2*I1+3*I2*H-B(j)*(6*I2+6*I3*H+Is_underlined*(6*ksi0*alpha^2*H+12*omega)))*utms;
%             F_notunderline(j) = Fpre_underline(j)+Fvis_underline(j);