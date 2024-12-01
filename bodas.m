% Drawing asymptotic bode

% Auralius Manurung
% manurung.auralius@gmail.com

% A system defined as:
%
%          (s+z1)(s+z2) ... (s+zm)
% G(s) = K -----------------------
%          (s+p1)(s+p2) ... (s+pn)     
%
% To draw the asymptotic Bode plots, we use the following script:
%
% bodas([z1 z2 ... zm], [p1 p2 ... pn], K)
%
% We can also define the frequency range as follows:       
%
% bodas([z1 z2 ... zm], [p1 p2 ... pn], K)
%
% Please, check the PDF file for more detailed examples.


function [G, w] = bodas(sys)
q = zpk(sys); 
z = q.Z{1};
p = q.P{1};
k = q.K;

i = 1;
while(~isempty(z))
    if isreal(z(i)) == false
        z(find(z == conj(z(i)))) = [];
    end
    z(i) = -z(i);
    i = i + 1;
    if i > length(z)
        break;
    end
end

i = 1;
while(~isempty(p))
    if isreal(p(i)) == false
        p(find(p == conj(p(i)))) = [];
    end
        p(i) = -p(i);

    i = i + 1;
    if i > length(p)
        break;
    end
end

z = sort(z);
p = sort(p);

W = unique(nonzeros([z; p]));
for j  = 1 : length(W)
    if isreal(W(j)) == false
        W(j) = ceil(log10(sqrt(W(j) * conj(W(j)))));
    else
        W(j) = ceil(log10(W(j)));
    end
end

W = sort(W);

if isempty(W)
    W = 1;
end

wmin = W(1)-2;
wmax = W(end)+2;
omega = logspace(wmin, wmax, 1000);

s = tf('s');

% -------------------------------------------------------------------------
% Gain plot
% -------------------------------------------------------------------------

figure;
set(gcf,'units','normalized','outerposition',[0 0 1 1])

[ha, ] = tight_subplot(2, 1, [0.05, 0.05]);

axes(ha(1));
hold on;
grid on;

legend_text = cell(length(z) + length(p) + length(k) + 2, 1);
ctr = 1;

A = zeros(length(z), length(omega));
B = zeros(length(p), length(omega));
C = zeros(1, length(omega));

G = 1;

for i = 1:length(z)
    if z(i) == 0 % Zero at the origin
        offset = log10(1/omega(1))*20;
        for j = 1:length(omega)
            A(i,j)=20*log10(omega(j)/omega(max(j-1,1)))+ A(i,max(j-1,1));
        end
        A(i,:)=A(i,:) - offset;
    elseif isreal(z(i)) == true
        for j = 1:length(omega)
            if omega(j) >= z(i)
                A(i,j)=20*log10(omega(j)/omega(max(j-1,1)))+ A(i,max(j-1,1));
            else
                A(i,j)=20*log10(z(i));
            end
        end
    elseif isreal(z(i)) == false
        peak = false;
        wn = sqrt(z(i)*conj(z(i)));
        zeta = (z(i)+conj(z(i)))/2/wn;
        for j = 1:length(omega)
            if omega(j) >= wn
                A(i,j)=40*log10(omega(j)/omega(max(j-1,1)))+ A(i,max(j-1,1));
                if peak == false
                    peak_idx = j;
                    peak = true;
                end
            else
                A(i,j)=40*log10(wn);
            end
        end
        if zeta < 0.5
            A(i,peak_idx)=A(i,peak_idx)-20*log10(1/(2*zeta));
        end
        
    end
    
    plot(omega,A(i,:), 'LineWidth', 2);
    if z(i) == 0
        legend_text{ctr} = "s";
        G = G * s;
    else
        if isreal(z(i)) == true
            legend_text{ctr} = "s+"+ num2str(z(i));
            G = G * (s+z(i));
        elseif isreal(z(i)) == false
            legend_text{ctr} = "(s+"+ num2str(z(i))+") (s+"+ num2str(conj(z(i)))+")" ;
            G = G * (s+z(i)) * (s+conj(z(i)));
        end
    end
    ctr = ctr + 1;
end

for i = 1:length(p)
    if p(i) == 0  % Pole at the origin
        offset = log10(1/omega(1))*20;
        for j = 1:length(omega)
            B(i,j)=-20*log10(omega(j)/omega(max(j-1,1)))+ B(i,max(j-1,1));
        end
        B(i,:)=B(i,:)+offset;
    elseif isreal(p(i)) == true
        for j = 1:length(omega)
            if omega(j) >= p(i)
                B(i,j)=-20*log10(omega(j)/omega(max(j-1,1)))+ B(i,max(j-1,1));
            else
                B(i,j)=-20*log10(p(i));
            end
        end
    elseif isreal(p(i)) == false
        peak = false;
        wn = sqrt(p(i)*conj(p(i)));
        zeta = (p(i)+conj(p(i)))/2/wn;
        for j = 1:length(omega)
            if omega(j) >= wn
                B(i,j)=-40*log10(omega(j)/omega(max(j-1,1)))+ B(i,max(j-1,1));
                if peak == false
                    peak_idx = j;
                    peak = true;
                end
            else
                B(i,j)=-40*log10(wn);
            end
        end
        if zeta < 0.5
            B(i,peak_idx)=B(i,peak_idx)+20*log10(1/(2*zeta));
        end
    end
    
    
    plot(omega,B(i,:), 'LineWidth', 2);
    if p(i) == 0
        legend_text{ctr} = "1/s";
        G = G * 1/s;
    elseif isreal(p(i)) == false
        legend_text{ctr} = "1/( (s+"+ num2str(p(i))+") (s+"+ num2str(conj(p(i)))+")" +" )" ;
        G = G * 1/( (s+p(i)) * (s+conj(p(i))) );
    elseif isreal(p(i)) == true
        legend_text{ctr} = "1/(s+"+ num2str(p(i))+")" ;
        G = G * 1/(s+p(i));
    end
    ctr = ctr + 1;
end

% The gain
if k ~= 0
    kDb = 20*log10(abs(k));
    for j = 1:length(omega)
        C(j)=kDb;
    end
    
    plot(omega, C, 'LineWidth', 2);
    legend_text{ctr} = num2str(k);
    G = G * k;
end

[mag,phase, wout] = bode(G, {10^wmin, 10^wmax});       

% Combinations
M = sum([sum(A,1);sum(B,1);C],1);
plot(omega, M, 'LineWidth',3, 'LineStyle', '-', 'Color', [.7, .7,.7]);
plot(wout, 20*log10(squeeze(mag)), 'LineWidth', 2, 'LineStyle', ':', 'Color', 'k');

legend_text{end-1} = "Asymptotic Bode";
legend_text{end}= "Exact Bode";
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

for i = 1:length(z)
    if z(i) == 0
        A(i,:)=ones(1,length(omega))*90;
    elseif isreal(z(i)) == true
        for j = 1:length(omega)
            if omega(j) > z(i)/10 && omega(j) <= z(i)*10
                A(i,j)=min(45*log10(omega(j)/omega(j-1))+ A(i,j-1),90);
            elseif omega(j) <= z(i)/10
                A(i,j)=0;
            elseif omega(j) > z(i)*10
                A(i,j)=90;
            end
        end
    elseif isreal(z(i)) == false
        wn = sqrt(z(i)*conj(z(i))); 
        zeta = (z(i)+conj(z(i)))/2/wn;
        delta = 10^zeta;

        for j = 1:length(omega)
            if (omega(j) > wn/delta) && (omega(j) <= wn*delta)
                A(i,j)=90/zeta*log10(omega(j)/omega(j-1))+ A(i,j-1);
            elseif omega(j) <= wn/delta
                A(i,j)=0;
            elseif omega(j) > wn*delta
                A(i,j)=180;
            end
        end
    end
    plot(omega,A(i,:), 'LineWidth', 2);
end

for i = 1:length(p)
    if p(i) == 0
        B(i,:)=ones(1,length(omega))*-90;
    elseif isreal(p(i)) == true
        for j = 1:length(omega)
            if omega(j) > p(i)/10 && omega(j) <= p(i)*10
                B(i,j)=max(-45*log10(omega(j)/omega(j-1))+ B(i,j-1), -90);
            elseif omega(j) <= p(i)/10
                B(i,j)=0;
            elseif omega(j) > p(i)*10
                B(i,j)=-90;
            end
        end
     elseif isreal(p(i)) == false
        wn = sqrt(p(i)*conj(p(i))); 
        zeta = (p(i)+conj(p(i)))/2/wn;
        delta = 10^zeta;

        for j = 1:length(omega)
            if (omega(j) > wn/delta) && (omega(j) <= wn*delta)
                B(i,j)=-90/zeta*log10(omega(j)/omega(j-1))+ B(i,j-1);
            elseif omega(j) <= wn/delta
                B(i,j)=0;
            elseif omega(j) > wn*delta
                B(i,j)=-180;
            end
        end
    end    
    plot(omega,B(i,:), 'LineWidth', 2);
end

% The gain
if k ~= 0 && nargin > 2
    for j = 1:length(omega)
        if (k > 0)
            C(j)=0;
        else
            C(j)=-180;
        end
    end
    plot(omega, C, 'LineWidth', 2);
end

% Combinations
S = sum([sum(A,1);sum(B,1);C],1);
if min(S) < 0
    S = S + 180;
end
plot(omega, S, 'LineWidth',3, 'LineStyle', '-', 'Color', [.7, .7,.7]);
plot(wout, squeeze(phase), 'LineWidth', 2, 'LineStyle', ':', 'Color', 'k');
ylabel('Phase (degrees)')
xlabel('Frequency (log scale)')
set(gca, 'XScale', 'log');

figure;
[hb, ] = tight_subplot(2, 1, [0.05, 0.05]);

axes(hb(1));
hold on;
grid on;
plot(omega, M, 'LineWidth',3, 'LineStyle', '-', 'Color', [.7, .7,.7]);
plot(wout, 20*log10(squeeze(mag)), 'LineWidth', 2, 'LineStyle', ':', 'Color', 'k');
legend({'Asymptotic Bode', 'Exact Bode'}, 'Location', 'best');
ylabel('Magnitude (dB)');
set(gca, 'XScale', 'log');

axes(hb(2));
hold on;
grid on;
plot(omega, S, 'LineWidth',3, 'LineStyle', '-', 'Color', [.7, .7,.7]);
plot(wout, squeeze(phase), 'LineWidth', 2, 'LineStyle', ':', 'Color', 'k');
legend({'Asymptotic Bode', 'Exact Bode'}, 'Location', 'best');
ylabel('Phase (degrees)')
xlabel('Frequency (log scale)')
set(gca, 'XScale', 'log');

w ={10^wmin 10^wmax};
end