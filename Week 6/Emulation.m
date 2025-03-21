%% EGH445 - Modern Control
%--------------------------------------------------------------------------
% Modern Control
% **** Support Functions
%     **** Emulation Design
%
% Author: Aaron McFadyen
%--------------------------------------------------------------------------

% Design or emulate discrete controllers from continous controllers.

% G(s): transmembrane voltage potential from rest, 
% t: Vector of sampling times 
% Methods: Cell array of strings defining the emulation approach (see c2d.m) 

%% Emulation Function

function [F] = Emulation(Gc,Gs,Fz,t,methods,colors)

T = diff(t(1:2));
F = cell(numel(T),1);

% Get Feedback Transfer Functions (unity feedback, discrete)
Hz.tf = tf([0 1],[0 1],T); 
Hs.tf = tf([0 1],[0 1]);      

% Continuous Feedback Control
Gs.tfclosed     = feedback(series(Gc.tf,Gs.tf),Hs.tf);
Gs.polesclosed  = roots(cell2mat(Gs.tfclosed.den));
Gs.stepinfo     = stepinfo(Gs.tfclosed);

for i=1:length(methods)

    % Discrete Feedback Control
    Fc.tf           = c2d(Gc.tf,T,methods{i});
    Fz.tfclosed     = feedback(series(Fc.tf,Fz.tf),Hz.tf);
    Fz.polesclosed  = roots(cell2mat(Fz.tfclosed.den));
    Fz.zerosclosed  = roots(cell2mat(Fz.tfclosed.num));
    Fz.stepinfo     = stepinfo(Fz.tfclosed);
    %[Wn,zeta,P] = damp(Fz.tfclosed)

    F{i} = Fz;
    
    % Closed Loop Response (discrete - control)
    figure(101);hold on
    [x2,t2] = step(Fz.tfclosed,t,horzcat(colors.full{i},'-'));       % State 2 (Output State)
    stairs(t2,x2,horzcat(colors.full{i},'-'));
    %x1 = cumtrapz(t2,x2);                                      % State 1 (Non-Output State)
    %stairs(t2,x1,'-','Color', lightcolors(i,:),'LineWidth',2);
    
    if i == length(methods)       
        % Closed Loop Response (continuous - control)
        figure(101); hold on
        [x2,t2] = step(Gs.tfclosed, t);
        plot(t2,x2,'k-');
        figure(101);title('Step Response (Simulink and MATLAB) - Closed Loop', 'interpreter','latex');
        legend([methods, "continuous"])
        %x1 = cumtrapz(t2,x2);
        %plot(t2,x1,'-','Color', [0.8 0.8 0.8],'LineWidth',2);
        
        figure(303); hold on
        rlocus(series(Gc.tf,Gs.tf),[1 2 3 4 5]);
        figure(303);title('Root Locus - Continuous', 'interpreter','latex')
        
        figure(404); hold on
        rlocus(series(Fz.tf,Fc.tf),[1 2 3 4 5]);
        figure(404);title('Root Locus - Discrete', 'interpreter','latex');
        
        % Open Loop Response (continuous - no control)
        figure(202); hold on
        [x2,t2] = step(Gs.tf, t);                   % (output 2)
        plot(t2,x2,'k-');   
        
        %x1 = cumtrapz(t2,x2);                      % (output 1)
        %plot(t2,x1,'-','Color', [0.8 0.8 0.8 1.0],'LineWidth',2);
        
        % Open Loop Response (discrete - no control)
        [x2,t2] = step(Fz.tf,t);                    % (output 2)
        stairs(t2,x2,'k.-');   
        %x1 = cumtrapz(t2,x2);                      % (output 1)
        %stairs(t2,x1,'.-','Color', [0.8 0.8 0.8],'LineWidth',2);
        figure(202);title('Step Response - Open Loop', 'interpreter','latex');   
        xlabel('Time (s)'); ylabel('Amplitude')
    end
end

%% Plotting/Checking Resources
%
%         Tsim = 1/100;
%         sim('Equivalence')
%         figure(505); hold on;
%         plot(ycont.Time,ycont.Data,'k-');
%         plot(ydisc.Time,ydisc.Data,horzcat(colors{i},'--'));
%         plot(ydiscdisc.Time,ydiscdisc.Data,'-','Color',horzcat(lightcolors(i,:),1.0));
%         figure(505);title('Step Response (Simulink) - Closed Loop', 'interpreter','latex'); 
%         xlabel('Time (s)', 'interpreter','latex');
%         ylabel('Output (m)', 'interpreter','latex')
%         %ExportFigJpg(horzcat('tustin'));
%         %ExportFigJpg(horzcat('zoh'));
%         %ExportFigJpg(horzcat('mpz'));
%         
%         
%         figure(101)
%         plot(ydisc.Time,ydisc.Data,horzcat(colors{i},'--'));
%         plot(ydiscdisc.Time,ydiscdisc.Data,'-','Color',horzcat(lightcolors(i,:),1.0));
%         kk=1;
   