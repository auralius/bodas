clc
close all
s = tf('s');
G = (-s^2-0.2*s-4.01) / s
[G, w] = bodas(G);
