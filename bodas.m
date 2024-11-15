% Drawing asymptotic bode from 0.1 to 1000 Hz, by default
% Does not work with complex poles

% Auralius Manurung
% manurunga@yandex.com

% A system defined as:
%
%          (s+z1)(s+z2) ... (s+zm)
% G(s) = K -----------------------
%          (s+p1)(s+z2) ... (p+zn)     
%
% To draw the asymptotic Bode plots, we use the following script:
%
% bodas([z1 z2 ... zm], [p1 p2 ... pn], K)
%
% We can also define the frequency range as follows:       
%
% bodas([z1 z2 ... zm], [p1 p2 ... pn], K, [-2 3])
% Here, at the last argument, -2 corresponds to 10^-2 and 3 corresponds to 10^3.
%
% Please, check the PDF file for more detailed examples.


function [G, w] = bodas(Z, P, K)
Z = sort(Z);
P = sort(P);

W = sort(unique(nonzeros([Z, P])));
wmin = W(1)/100;
wmax = W(end)*100;
w = {wmin, wmax};
omega = logspace(log10(wmin), log10(wmax), 1000);

s = tf('s');

% -------------------------------------------------------------------------
% Gain plot
% -------------------------------------------------------------------------

figure;
set(gcf,'units','normalized','outerposition',[0 0 1 1])

[ha, pos] = tight_subplot(2, 1, [0.05, 0.05]);

axes(ha(1));
hold on;
grid on;

legend_text = cell(length(Z) + length(P) + 3, 1);
ctr = 1;

A = zeros(1,length(omega));
B = zeros(1,length(omega));
C = zeros(1,length(omega));
G = 1;

for i = 1:length(Z)
    if Z(i) == 0 % Zero at the origin
        offset = log10(1/omega(1))*20;
        for j = 1:length(omega)
            A(i,j)=20*log10(omega(j)/omega(max(j-1,1)))+ A(i,max(j-1,1));
        end
        A(i,:)=A(i,:) - offset;
    else
        for j = 1:length(omega)
            if omega(j) >= Z(i)
                A(i,j)=20*log10(omega(j)/omega(max(j-1,1)))+ A(i,max(j-1,1));
            else
                A(i,j)=20*log10(Z(i));
            end
        end
    end
    
    plot(omega,A(i,:), 'LineWidth', 2);
    if Z(i) == 0
        legend_text{ctr} = "G(s)=s";
        G = G * s;
    else
        legend_text{ctr} = "G(s)=s+"+ num2str(Z(i));
        G = G * s+Z(i);
    end
    ctr = ctr + 1;
end

for i = 1:length(P)
    
    if P(i) == 0  % Pole at the origin
        offset = log10(1/omega(1))*20;
        for j = 1:length(omega)
            B(i,j)=-20*log10(omega(j)/omega(max(j-1,1)))+ B(i,max(j-1,1));
        end
        B(i,:)=B(i,:)+offset;
    else
        for j = 1:length(omega)
            if omega(j) >= P(i)
                B(i,j)=-20*log10(omega(j)/omega(max(j-1,1)))+ B(i,max(j-1,1));
            else
                B(i,j)=-20*log10(P(i));
            end
        end
    end
    
    
    plot(omega,B(i,:), 'LineWidth', 2);
    if P(i) == 0
        legend_text{ctr} = "G(s)=1/s";
        G = G * 1/s;
    else
        legend_text{ctr} = "G(s)=1/(s+"+ num2str(P(i))+")" ;
        G = G * 1/(s+P(i));
    end
    ctr = ctr + 1;
end

% The gain
if K ~= 0 && nargin > 2
    kDb = 20*log10(abs(K));
    for j = 1:length(omega)
        C(j)=kDb;
    end
    
    plot(omega, C, 'LineWidth', 2);
    legend_text{ctr} = "G(s)="+ num2str(K);
    G = G * K;

[mag,phase, wout] = bode(G, {wmin, wmax});       

% Combinations
plot(omega, sum([sum(A,1);sum(B,1);C],1), 'LineWidth',6, 'LineStyle', '-');
plot(wout, 20*log10(squeeze(mag)), 'LineWidth', 3, 'LineStyle', '-')

legend_text{end-1} = "Asymptotic Bode";
legend_text{end}= "Actual Bode";
legend(legend_text, 'Location', 'best');
ylabel('Magnitude (dB)');
set(gca, 'XScale', 'log');

% -------------------------------------------------------------------------
% Phase plot
% -------------------------------------------------------------------------

axes(ha(2));
hold on;
grid on;

A = zeros(1,length(omega));
B = zeros(1,length(omega));
C = zeros(1,length(omega));

for i = 1:length(Z)
    if Z(i) == 0
        A(i,:)=ones(1,length(omega))*90;
    else
        for j = 1:length(omega)
            if omega(j) > Z(i)/10 && omega(j) <= Z(i)*10
                A(i,j)=min(45*log10(omega(j)/omega(j-1))+ A(i,j-1),90);
            elseif omega(j) <= Z(i)/10
                A(i,j)=0;
            elseif omega(j) > Z(i)*10
                A(i,j)=90;
            end
        end
    end
    plot(omega,A(i,:), 'LineWidth', 2);
end

for i = 1:length(P)
     if P(i) == 0
        B(i,:)=ones(1,length(omega))*-90;
    else
        for j = 1:length(omega)
            if omega(j) > P(i)/10 && omega(j) <= P(i)*10
                B(i,j)=max(-45*log10(omega(j)/omega(j-1))+ B(i,j-1), -90);
            elseif omega(j) <= P(i)/10
                B(i,j)=0;
            elseif omega(j) > P(i)*10
                B(i,j)=-90;
            end

        end
     end
    
    plot(omega,B(i,:), 'LineWidth', 2);
end

% The gain
if K ~= 0 && nargin > 2
    for j = 1:length(omega)
        if (K > 0)
            C(j)=0;
        else
            C(j)=-180;
        end
    end
    plot(omega, C, 'LineWidth', 2);
end

% Combinations
plot(omega, sum([sum(A,1);sum(B,1);C],1), 'LineWidth',6, 'LineStyle', '-');
ylabel('Phase (degrees)')
xlabel('Frequency (log scale)')
set(gca, 'XScale', 'log');

end