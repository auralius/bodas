% Drawing asymptotic bode from 0.1 to 1000 Hz
% Does not work with imaginary poles

% Auralius Manurung
% manurung.auralius@gmail.com

% How to use:
%
% A system with gain = [0.1], and zeros = [10 100], poles = [1]
%   bodas([10 100], [1], 0.1)
%
% A system with gain = -10, zeros = [0], poles = [1 1 10]. The bode is
% drawn in frequency 0.002 to 1000 rad/s
%   bodas([0],[1 1 10], -10, [-2 3]
% -2 corresponds to 10^-2 and 3 corresponds to 10^3.

function bodas(Z, P, K, range)

if nargin > 3
    % Frequency = from 10e-1 to 10e3
    omega = logspace(range(1), range(2), 100);
else
    % Frequency = from 10e-1 to 10e3
    omega = logspace(-1, 3, 100);
end

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Gain plot
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

figure;
subplot(2,1,1)
hold on;
grid on;

legend_text = cell(length(Z) + length(P ), 1);
ctr = 1;

A = zeros(1,length(omega));
B = zeros(1,length(omega));
C = zeros(1,length(omega));
K = 0;

for i = 1:length(Z)
    if Z(i) == 0 % Zero at the origin
        offset = log10(1/omega(1))*20;
        for j = 1:length(omega)
            A(i,j)=20*log10(omega(j)/omega(max(j-1,1)))+ A(i,max(j-1,1));
        end
        A(i,:)=A(i,:) - offset;
    else
        for j = 1:length(omega)
            if omega(j) > Z(i)
                A(i,j)=20*log10(omega(j)/omega(max(j-1,1)))+ A(i,max(j-1,1));
            else
                A(i,j)=0;
            end
        end
    end
    
    plot(omega,A(i,:));
    if Z(i) == 0
        legend_text{ctr} = "G(s)= s";
    else
        legend_text{ctr} = "G(s)=" + "(s+"+ num2str(Z(i))+")/"+num2str(Z(i));
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
            if omega(j) > P(i)
                B(i,j)=-20*log10(omega(j)/omega(j-1))+ B(i,j-1);
            else
                B(i,j)=0;
            end
        end
    end
    
    
    plot(omega,B(i,:));
    if P(i) == 0
        legend_text{ctr} = "G(s)= 1/s";
    else
        legend_text{ctr} = "G(s)=" + num2str(P(i)) + "/(s+"+ num2str(P(i))+")" ;
    end
    ctr = ctr + 1;
end

% The gain
if K ~= 0 && nargin > 2
    kDb = 20*log10(abs(K));
    for j = 1:length(omega)
        C(j)=kDb;
    end
    
    plot(omega,C);
    legend_text{ctr} = "G(s)="+ num2str(K);
end

% Combinations
plot(omega, sum([sum(A,1);sum(B,1);C],1), 'LineWidth',3, 'LineStyle', '-.');
set(gca, 'XScale', 'log');
legend(legend_text, 'Location', 'best');
ylabel('Magnitude (dB)');

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Phase plot
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

subplot(2,1,2)
hold on;
grid on;

A = zeros(1,length(omega));
B = zeros(1,length(omega));
C = zeros(1,length(omega));
K = 0;

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
    plot(omega,A(i,:));
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
    
    plot(omega,B(i,:));
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
    plot(omega, C);
end

% Combinations
plot(omega, sum([sum(A,1);sum(B,1);C],1), 'LineWidth',3, 'LineStyle', '-.');
set(gca, 'XScale', 'log');
ylabel('Phase (degrees)')
xlabel('Frequency (log scale)')

end