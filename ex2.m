clc
close all
s = tf('s');
G = -1 / (s^3 + 0.2*s^2 + 4.01*s)
[G, w] = bodas(G);

