clc
clear all
close all

% ---This file calculates
%    1. the critical angle for snapping;
%    2. the range of thetaB for S shape.
%
% ---The step size of thetaB is set to be 1 degrees (thetaB_step = pi/180).


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

thetaB_step = pi/180;  % step size of thetaB

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
% we run the section when the thetaB_step is set to be 1 degree 
thetaBcrit_r = thetaBrange1(sw_r);        %  switching angle at right hand side  [rad]
thetaBcrit_l = thetaBrange2(sw_l)-2*pi;   %  switching angle at left  hand side  [rad]

thetaBcrit_r_deg = thetaBcrit_r/pi*180;    %  switching angle at right hand side  [degree]
thetaBcrit_l_deg = thetaBcrit_l/pi*180;    %  switching angle at left  hand side  [degree]

disp(['The critical angles for snapping is ',num2str(thetaBcrit_r_deg),' degree and ',num2str(thetaBcrit_l_deg),' degree.'])

%%  calculate S-ratio (only analyze right hand side which has S-shape)
fontsize = 13;
S_crit = 0.3;  % set the critical value of S-ratio
xplot = 0:L/100:L;
for k = 1:length(a1_sol_r)
    yplot = a1_sol_r(k)*( 1-cos(2*pi*xplot/L)  )+ a2_sol_r(k)*( 1-2*xplot/L - cos( 2.86*pi*xplot/L ) + 2/2.86/pi*sin(2.86*pi*xplot/L )  );
    ymax = abs(max(yplot));   
    ymin = abs(min(yplot));
    Sratio(k) = min([ymax,ymin])/max([ymax,ymin]);
    
    figure(1)
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
% we run the section when the thetaB_step is set to be 1 degree 
thetaB_Srange = thetaBrange1(Sratio>=S_crit);
Srange = thetaB_Srange(1)-thetaB_Srange(end);
S_start_deg = thetaB_Srange(1)/pi*180;
S_end_deg = thetaB_Srange(end)/pi*180;
Srange_deg = S_start_deg-S_end_deg;

disp(['The range of thetaB for S-Shape is from ',num2str(S_start_deg),' degree to ',num2str(S_end_deg),' degree.'])

