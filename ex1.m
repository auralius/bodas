%
%                  -100 s
%  G(s) = ------------------------
%         s^3 + 12 s^2 + 21 s + 10
%
%                         s+0
%  G(s) = -100 * ----------------------
%                (s+1) * (s+1) * (s+10)

clc
close all

[G, w] = bodas([0], [1 1 10], -100);

