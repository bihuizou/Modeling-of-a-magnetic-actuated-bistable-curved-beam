clc
clear all
close all

% ---This file demonstrates 
%    1. the multimodal deformation of the bistable curved beam under various directions of magnetic fields;
%    2. the calculation of magnetic torques;
%    3. the calculation of S-ratio;
%    4. the multimodal deformation in the seceond step actuation.
%
% ---The step size of thetaB is set to be 10 degrees (thetaB_step = pi/18) in this file.
%
% ---To analyze the critical angle of snapping or the range of S-Shape, we
%    need to use a smaller step size (eg. thetaB_step = pi/180)

%% Initialize

% Input Parameter 
L = 15*1e-3;        % span of the beam [m]
h = 5*1e-3;         % amplitude of the beam (apex height) [m]
L1 = 2*1e-3;        % position of the left magnet [m]
L2 = 13*1e-3;       % position of the right magnet [m]
b = 3.5*1e-3;       % width of the beam [m]
t = 0.6*1e-3;       % thickness of the beam [m]
II = b*t^3/12;      % area moment of inertia of the beam 
EE = 3*1e6;         % modulus of beam material [Pa]
B = 50*1e-3;        % intensity of magnetic field [T]
m1 = 0.1140;        % magnetic moment of the left magnet [A m^2]
m2 = m1/5;          % magnetic moment of the right magnet [A m^2]

thetaB_step = pi/18;  % step size of thetaB

% stationary condition with respect to a1 and a2
syms a1 a2 thetaB
pp = pi;  
Eq1 = (EE*II*(16*a1*pp^4 - (327184*a2*pp^2)/10449 - 8*h*pp^4 + (1144*a2*pp^2*cos((43*pp)/50))/43 + (1144*a2*pp^2*cos((243*pp)/50))/243 + 4*a1*pp^3*sin(4*pp) + (40898*a2*pp^3*sin((43*pp)/50))/1075 + (40898*a2*pp^3*sin((243*pp)/50))/6075 - 2*h*pp^3*sin(4*pp)))/(2*L^3) + (2*B*m1*pp*sin((2*L1*pp)/L)*sin(a2*((2*cos((143*L1*pp)/(50*L)))/L - 2/L + (143*pp*sin((143*L1*pp)/(50*L)))/(50*L)) - thetaB + (2*a1*pp*sin((2*L1*pp)/L))/L))/L + (2*B*m2*pp*sin((2*L2*pp)/L)*sin(pp - thetaB + a2*((2*cos((143*L2*pp)/(50*L)))/L - 2/L + (143*pp*sin((143*L2*pp)/(50*L)))/(50*L)) + (2*a1*pp*sin((2*L2*pp)/L))/L))/L - (2*EE*b*t*(((81796*a2)/10449 - 4*a2*cos(2*pp) - (200*a2*cos((43*pp)/50))/43 + (200*a2*cos((243*pp)/50))/243)/(2*L) - (2*a1*pp^2)/L + (pp*(a1*sin(4*pp) - (286*a2*sin((43*pp)/50))/43 + (286*a2*sin((243*pp)/50))/243))/(2*L))*((4*a2^2*cos((143*pp)/50) - a2^2*cos((143*pp)/25) - (81796*a1*a2)/10449 + 3*a2^2 + 4*a1*a2*cos(2*pp) + (200*a1*a2*cos((43*pp)/50))/43 - (200*a1*a2*cos((243*pp)/50))/243)/(2*L) + ((50*a2^2*sin((143*pp)/25))/143 - (400*a2^2*sin((143*pp)/50))/143)/(2*L*pp) - (pp*((a1^2*sin(4*pp))/2 + (143*a2^2*sin((143*pp)/25))/200 - (286*a1*a2*sin((43*pp)/50))/43 + (286*a1*a2*sin((243*pp)/50))/243))/(2*L) + (pp^2*(2*a1^2 + (20449*a2^2)/5000))/(2*L) - (h^2*pp*(4*pp - sin(4*pp)))/(16*L)))/((pp*(4*pp - sin(4*pp))*h^2)/(16*L) + L);
Eq2 = B*m1*sin(a2*((2*cos((143*L1*pp)/(50*L)))/L - 2/L + (143*pp*sin((143*L1*pp)/(50*L)))/(50*L)) - thetaB + (2*a1*pp*sin((2*L1*pp)/L))/L)*((2*cos((143*L1*pp)/(50*L)))/L - 2/L + (143*pp*sin((143*L1*pp)/(50*L)))/(50*L)) + B*m2*sin(pp - thetaB + a2*((2*cos((143*L2*pp)/(50*L)))/L - 2/L + (143*pp*sin((143*L2*pp)/(50*L)))/(50*L)) + (2*a1*pp*sin((2*L2*pp)/L))/L)*((2*cos((143*L2*pp)/(50*L)))/L - 2/L + (143*pp*sin((143*L2*pp)/(50*L)))/(50*L)) + (EE*II*((20449*a2*pp^2)/1250 - (327184*a1*pp^2)/10449 + (418161601*a2*pp^4)/6250000 + (163592*h*pp^2)/10449 - (143*a2*pp*sin((143*pp)/25))/25 + (1144*a1*pp^2*cos((43*pp)/50))/43 + (20449*a2*pp^2*cos((143*pp)/25))/1250 + (1144*a1*pp^2*cos((243*pp)/50))/243 - (572*h*pp^2*cos((43*pp)/50))/43 - (572*h*pp^2*cos((243*pp)/50))/243 + (40898*a1*pp^3*sin((43*pp)/50))/1075 + (2924207*a2*pp^3*sin((143*pp)/25))/250000 + (40898*a1*pp^3*sin((243*pp)/50))/6075 - (20449*h*pp^3*sin((43*pp)/50))/1075 - (20449*h*pp^3*sin((243*pp)/50))/6075))/(2*L^3) + (2*EE*b*t*((6*a2 - (81796*a1)/10449 + 4*a1*cos(2*pp) + (200*a1*cos((43*pp)/50))/43 - 2*a2*cos((143*pp)/25) + 8*a2*cos((143*pp)/50) - (200*a1*cos((243*pp)/50))/243)/(2*L) + (20449*a2*pp^2)/(5000*L) + ((100*a2*sin((143*pp)/25))/143 - (800*a2*sin((143*pp)/50))/143)/(2*L*pp) - (pp*((143*a2*sin((143*pp)/25))/100 - (286*a1*sin((43*pp)/50))/43 + (286*a1*sin((243*pp)/50))/243))/(2*L))*((4*a2^2*cos((143*pp)/50) - a2^2*cos((143*pp)/25) - (81796*a1*a2)/10449 + 3*a2^2 + 4*a1*a2*cos(2*pp) + (200*a1*a2*cos((43*pp)/50))/43 - (200*a1*a2*cos((243*pp)/50))/243)/(2*L) + ((50*a2^2*sin((143*pp)/25))/143 - (400*a2^2*sin((143*pp)/50))/143)/(2*L*pp) - (pp*((a1^2*sin(4*pp))/2 + (143*a2^2*sin((143*pp)/25))/200 - (286*a1*a2*sin((43*pp)/50))/43 + (286*a1*a2*sin((243*pp)/50))/243))/(2*L) + (pp^2*(2*a1^2 + (20449*a2^2)/5000))/(2*L) - (h^2*pp*(4*pp - sin(4*pp)))/(16*L)))/((pp*(4*pp - sin(4*pp))*h^2)/(16*L) + L);



%% magnetic field pointing to right hand side
thetaB_0 = pi/2;
thetaBrange1 = thetaB_0:(-thetaB_step):(-pi/2);  
d = 0;

for k = 1:length(thetaBrange1)
% specify angle
    thetaB_value = thetaBrange1(k);
    Eq11 = subs(Eq1,thetaB,thetaB_value);
    Eq21 = subs(Eq2,thetaB,thetaB_value);
% specify initial guess    
    if k==1                  % at the first step we specify the initial guess
        a1_iniguess = h/2;
        a2_iniguess = 0;
    else                     % we set the solution from the last step as the initial guess at current step
        a1_iniguess = a1_sol_r(k-1);
        a2_iniguess = a2_sol_r(k-1);
    end
% solve a1 and a2
    [sol_1 sol_2] = vpasolve([Eq11,Eq21],[a1,a2],[a1_iniguess;a2_iniguess]);
    if isempty(double(sol_1))==1     
        % snapping happens
        [sol_1 sol_2] = vpasolve([Eq11,Eq21],[a1,a2],[-a1_sol_r(1);0]); % we reset the initial guess at this step and recalculate the value of a1 and a2 
        a1_sol_r(k) = double(sol_1);
        a2_sol_r(k) = double(sol_2);
        % save the position of switching point
        d = d+1; 
        sw_r(d) = k;   
    else
        % general case
        a1_sol_r(k) = double(sol_1);
        a2_sol_r(k) = double(sol_2);
    end
end

%%   magnetic field pointing to left hand side
thetaBrange2 = thetaB_0:(thetaB_step):(pi*3/2);
d=0;

for k = 1:length(thetaBrange2)
% specify angle
    thetaB_value = thetaBrange2(k);
    Eq11 = subs(Eq1,thetaB,thetaB_value);
    Eq21 = subs(Eq2,thetaB,thetaB_value);
% specify initial guess    
    if k==1                  % at the first step we specify the initial guess
        a1_iniguess = h/2;
        a2_iniguess = 0;
    else                     % we set the solution from the last step as the initial guess at current step
        a1_iniguess = a1_sol_l(k-1);
        a2_iniguess = a2_sol_l(k-1);
    end
% solve a1 and a2
    [sol_1 sol_2] = vpasolve([Eq11,Eq21],[a1,a2],[a1_iniguess;a2_iniguess]);
    if isempty(double(sol_1))==1      
        % snapping happens
        [sol_1 sol_2] = vpasolve([Eq11,Eq21],[a1,a2],[-a1_sol_l(1);a2_iniguess]);  % reset the initial guess and recalculate 
        a1_sol_l(k) = double(sol_1);
        a2_sol_l(k) = double(sol_2);
        % save the position of switching point
        d = d+1;
        sw_l(d) = k;
    else
        a1_sol_l(k) = double(sol_1);
        a2_sol_l(k) = double(sol_2);
    end
end


%% switching angle when snapping happens
% % we run the section when the thetaB_step is set to be 1 degree 
% thetaBcrit_r = thetaBrange1(sw_r);        %  switching angle at right hand side  [rad]
% thetaBcrit_l = thetaBrange2(sw_l)-2*pi;   %  switching angle at left  hand side  [rad]
% 
% thetaBcrit_r_deg = thetaBcrit_r/pi*180    %  switching angle at right hand side  [degree]
% thetaBcrit_l_deg = thetaBcrit_l/pi*180    %  switching angle at left  hand side  [degree]


%% combine solution (left and right)
thetaBrange2(thetaBrange2>pi) = thetaBrange2(thetaBrange2>pi)-2*pi;
thetaB1color = 38-(36*(thetaBrange1+pi)/(2*pi)+1);
thetaB2color = 38-(36*(thetaBrange2+pi)/(2*pi)+1);
thetaBrange = [thetaBrange1 thetaBrange2];
thetaBcolor = int8([thetaB1color thetaB2color]);
a1_sol = [a1_sol_r a1_sol_l];
a2_sol = [a2_sol_r a2_sol_l];


%% plot all (2 halves) solutions
fontsize = 17;
xplot = 0:L/100:L;
colormap( parula (  length(a1_sol)  )) % use color to represent the direction
cmap = colormap;
figure(1)
yplot = h/2*( 1-cos(2*pi*xplot/L)  );
plot(xplot,yplot,'k--','LineWidth',1.2)   % plot the initial shape
hold on
for k = 1:length(a1_sol)
        figure(1)
        yplot = a1_sol(k)*( 1-cos(2*pi*xplot/L)  )+ a2_sol(k)*( 1-2*xplot/L - cos( 2.86*pi*xplot/L ) + 2/2.86/pi*sin(2.86*pi*xplot/L )  );
        co = (k-1)*1/length(a1_sol);
        plot(xplot,yplot,'Color',[cmap(thetaBcolor(k),1),cmap(thetaBcolor(k),2),cmap(thetaBcolor(k),3)],'LineWidth',1.2)
        hold on
        xlabel('x [m]','FontSize',fontsize,'FontWeight','bold','FontName','Calibri')
        ylabel('y [m]','FontSize',fontsize,'FontWeight','bold','FontName','Calibri')
        xlim([0 0.015])
        ylim([-0.006 0.006])
        title ('Deformatino under various directions of B','FontSize',fontsize,'FontWeight','bold','FontName','Calibri')
        axis equal
        set(gca,'FontName','Calibri','FontSize',fontsize,'FontWeight',...
            'bold','YTick',...
            [-0.006 -0.004 -0.002 0 0.002 0.004 0.006],'YTickLabel',...
            {'-0.006','-0.004','-0.002','0','0.002','0.004','0.006'});
        
        figure(2)
        quiver(0,0,cos(thetaBrange(k)),sin(thetaBrange(k)),1,'Color',[cmap(thetaBcolor(k),1),cmap(thetaBcolor(k),2),cmap(thetaBcolor(k),3)],'LineWidth',1.2)
        hold on
        xlim([-1.5 1.5])
        ylim([-1.5 1.5])
        axis equal
        title ('Direction of magnetic field','FontSize',fontsize,'FontWeight','bold','FontName','Calibri')
end



%%  calculate S-ratio (only analyze right hand side which has S-shape)
fontsize = 13
S_crit = 0.3;  % set the critical value of S-ratio
for k = 1:length(a1_sol_r)
    yplot = a1_sol_r(k)*( 1-cos(2*pi*xplot/L)  )+ a2_sol_r(k)*( 1-2*xplot/L - cos( 2.86*pi*xplot/L ) + 2/2.86/pi*sin(2.86*pi*xplot/L )  );
    ymax = abs(max(yplot));   
    ymin = abs(min(yplot));
    Sratio(k) = min([ymax,ymin])/max([ymax,ymin]);
    
    figure(3)
    plot(thetaBrange1(k)/pi*180,Sratio(k),'o','MarkerFaceColor',[0,0,0],'MarkerEdgeColor',[0,0,0],'markersize',5)
    hold on
    ylim([0 1.2])
    xlim([-90 90])
    xlabel('Direction of magnetic field [\circ]','FontSize',fontsize,'FontWeight','bold','FontName','Calibri')
    ylabel('S ratio','FontSize',fontsize,'FontWeight','bold','FontName','Calibri')
    set(gca,'FontName','Calibri','FontSize',fontsize,'FontWeight',...
        'bold','XTick',...
    [-90 -60 -30 0 30 60 90],'XTickLabel',...
    {'-90','-60','-30','0','30','60','90'})
set(gca,'box','on');
end

plot([-90 90],[S_crit S_crit],'k--','LineWidth',1)  


%% find the range of S-shape
% % we run the section when the thetaB_step is set to be 1 degree 
% thetaB_Srange = thetaBrange1(Sratio>=S_crit);
% Srange = thetaB_Srange(1)-thetaB_Srange(end);
% S_start_deg = thetaB_Srange(1)/pi*180
% S_end_deg = thetaB_Srange(end)/pi*180
% Srange_deg = S_start_deg-S_end_deg


%% calculate the magnetic torques

T1cal = B*m1*sin(a2_sol.*((2*cos((143*L1*pp)/(50*L)))/L - 2/L + (143*pp*sin((143*L1*pp)/(50*L)))/(50*L)) - thetaBrange + (2.*a1_sol.*pp*sin((2*L1*pp)/L))/L);
T2cal = B*m2.*sin(pp - thetaBrange + a2_sol.*((2*cos((143*L2*pp)/(50*L)))/L - 2/L + (143*pp*sin((143*L2*pp)/(50*L)))/(50*L)) + (2.*a1_sol.*pp*sin((2*L2*pp)/L))/L);

fontsize = 17;

figure(4)
hold on
for k = 1:length(a1_sol)
      plot(thetaBrange(k)/pi*180,T1cal(k),'o','MarkerFaceColor',[0,0,0],'MarkerEdgeColor',[0,0,0],'markersize',3)
      hold on
      xlim([-180 180])
      xlabel('Direction of magnetic field [\circ]','FontSize',fontsize,'FontWeight','bold','FontName','Calibri')
      ylabel('Magnetic torque T1 [N·m]','FontSize',fontsize,'FontWeight','bold','FontName','Calibri')
        set(gca,'box','on','FontName','Calibri','FontSize',fontsize,'FontWeight',...
    'bold','XTick',...
    [-180  -120  -60  0  60  120  180],'XTickLabel',...
    {'-180','-120','-60','0','60','120','180'});
end

figure(5)
for k = 1:length(a1_sol)
      plot(thetaBrange(k)/pi*180,T2cal(k),'o','MarkerFaceColor',[0,0,0],'MarkerEdgeColor',[0,0,0],'markersize',3)
      hold on 
      xlim([-180 180])
      xlabel('Direction of magnetic field [\circ]','FontSize',fontsize,'FontWeight','bold','FontName','Calibri')
      ylabel('Magnetic torque T2 [N·m]','FontSize',fontsize,'FontWeight','bold','FontName','Calibri')
      set(gca,'box','on','FontName','Calibri','FontSize',fontsize,'FontWeight',...
    'bold','XTick',...
    [-180  -120  -60  0  60  120  180],'XTickLabel',...
    {'-180','-120','-60','0','60','120','180'});
end


%% analyze the second step ( right half)
thetaBrange3 = -pi/2:(thetaB_step):(pi/2);

for k = 1:length(thetaBrange3)
% specify angle
    thetaB_value = thetaBrange3(k);
    Eq11 = subs(Eq1,thetaB,thetaB_value);
    Eq21 = subs(Eq2,thetaB,thetaB_value);
% specify initial guess    
    if k==1                 % the initial guess is set to be the solution of C2 mode obtained from the first step  
        a1_iniguess = a1_sol_r(end);
        a2_iniguess = a2_sol_r(end);
    else                    % the initial guess is set to be the solution from last iterative step
        a1_iniguess = a1_sol_rev_r(k-1);
        a2_iniguess = a2_sol_rev_r(k-1);
    end
% solve a1 and a2
    [sol_1 sol_2] = vpasolve([Eq11,Eq21],[a1,a2],[a1_iniguess;a2_iniguess]);
    if isempty(double(sol_1))==1   % snapping happens 
        [sol_1 sol_2] = vpasolve([Eq11,Eq21],[a1,a2],[-a1_sol_rev_r(1);a2_iniguess]);
        a1_sol_rev_r(k) = double(sol_1);
        a2_sol_rev_r(k) = double(sol_2);
    else
        a1_sol_rev_r(k) = double(sol_1);
        a2_sol_rev_r(k) = double(sol_2);
    end
end

%% analyze the second step (left half)
thetaBrange4 = -pi/2:(-thetaB_step):(-pi*3/2);

for k = 1:length(thetaBrange4)
% specify angle
    thetaB_value = thetaBrange4(k);
    Eq11 = subs(Eq1,thetaB,thetaB_value);
    Eq21 = subs(Eq2,thetaB,thetaB_value);
% specify initial guess    
    if k==1       % the initial guess is set to be the solution of C2 mode obtained from the first step 
        a1_iniguess = a1_sol_l(end);
        a2_iniguess = a2_sol_l(end);
    else          % the initial guess is set to be the solution from last iterative step
        a1_iniguess = a1_sol_rev_l(k-1);
        a2_iniguess = a2_sol_rev_l(k-1);
    end
% solve a1 and a2
    [sol_1 sol_2] = vpasolve([Eq11,Eq21],[a1,a2],[a1_iniguess;a2_iniguess]);
    if isempty(double(sol_1))==1    % snapping happens
        [sol_1 sol_2] = vpasolve([Eq11,Eq21],[a1,a2],[-a1_sol_rev_l(1);a2_iniguess]);
        a1_sol_rev_l(k) = double(sol_1);
        a2_sol_rev_l(k) = double(sol_2);
    else
        a1_sol_rev_l(k) = double(sol_1);
        a2_sol_rev_l(k) = double(sol_2);
    end
end



%% combine the 2 halves - second step
thetaBrange4(thetaBrange4<=-pi) = thetaBrange4(thetaBrange4<=-pi)+2*pi;
thetaB3color = 38-(36*(thetaBrange3+pi)/(2*pi)+1);
thetaB4color = 38-(36*(thetaBrange4+pi)/(2*pi)+1);
thetaBrange_rev = [thetaBrange3 thetaBrange4];
thetaBcolor_rev = int8([thetaB3color thetaB4color]);
a1_sol_rev = [a1_sol_rev_r a1_sol_rev_l];
a2_sol_rev = [a2_sol_rev_r a2_sol_rev_l];

%% plot all solution - second step
xplot = 0:L/100:L;
colormap( parula (  length(a1_sol_rev)  ))
cmap = colormap;
figure(6)
yplot = -h/2*( 1-cos(2*pi*xplot/L)  );
plot(xplot,yplot,'k--','LineWidth',1.2)
hold on
for k = 1:length(a1_sol_rev)
        figure(6)
        yplot = a1_sol_rev(k)*( 1-cos(2*pi*xplot/L)  )+ a2_sol_rev(k)*( 1-2*xplot/L - cos( 2.86*pi*xplot/L ) + 2/2.86/pi*sin(2.86*pi*xplot/L )  );
        co = (k-1)*1/length(a1_sol_rev);
        plot(xplot,yplot,'Color',[cmap(thetaBcolor_rev(k),1),cmap(thetaBcolor_rev(k),2),cmap(thetaBcolor_rev(k),3)],'LineWidth',1.2)
        hold on
        xlabel('x [m]','FontSize',fontsize,'FontWeight','bold','FontName','Calibri')
        ylabel('y [m]','FontSize',fontsize,'FontWeight','bold','FontName','Calibri')
        xlim([0 0.015])
        ylim([-0.006 0.006])
        title ('Deformatino under various directions of B','FontSize',fontsize,'FontWeight','bold','FontName','Calibri')
        axis equal
        set(gca,'FontName','Calibri','FontSize',fontsize,'FontWeight',...
            'bold','YTick',...
            [-0.006 -0.004 -0.002 0 0.002 0.004 0.006],'YTickLabel',...
            {'-0.006','-0.004','-0.002','0','0.002','0.004','0.006'});
        
        figure(7)
        quiver(0,0,cos(thetaBrange_rev(k)),sin(thetaBrange_rev(k)),1,'Color',[cmap(thetaBcolor_rev(k),1),cmap(thetaBcolor_rev(k),2),cmap(thetaBcolor_rev(k),3)],'LineWidth',1.2)
        hold on
        xlim([-1.5 1.5])
        ylim([-1.5 1.5])
        axis equal
        title ('Direction of magnetic field','FontSize',fontsize,'FontWeight','bold','FontName','Calibri')
end
