%% Initialize MATLAB
clear, clc;

%% Define Variables
w_0 = 1;
Resolution = 1e-5;
AllQ = [10 100];
Terms = 10;
w = -Terms*w_0:Resolution:Terms*w_0;
k = -Terms:1:Terms;

figure;
for m = 1:length(AllQ)
    Q = AllQ(m);
    a = w_0/(2*Q);
    % Coeffs
    for j = 1:length(k)
        Eta(j).Fcn = (exp(-1i.*pi.*w./w_0).*sin(pi*w/w_0)) ./ ((pi.*w./w_0) .* (a + (w - w_0.*k(j))*1i));
    end
    
    % Add
    V_SH = zeros(size(w));
    for j=1:length(k)
        V_SH = V_SH + Eta(j).Fcn;
    end
    
    % Plot
    subplot(1,2,m), hold on;
    plot(w,abs(V_SH)/(2*Q),'b');
    plot(w,abs(1./(a+1i.*w))/(2*Q),'k--');
    axis tight;
    ylabel('$|V_{out}(j\omega)|$ (a.u.)','Interpreter','latex');
    xlim([-2*w_0,2*w_0]);
    ylim([0 1])
    xlabel('\omega/{\omega_0}');
    set(gcf, 'Units','centimeters', 'Position',[8 8 11 5.5])
end
annotation('textbox',[0.13, 0.81, 0.1, 0.1],'EdgeColor','none', 'String', ['Q = 10']);
annotation('textbox',[0.57, 0.81, 0.1, 0.1],'EdgeColor','none', 'String', ['Q = 100']);