s = tf('s');
sys = 20*(s+4)^2 / (s * (s^2+6*s+100));
bodas(sys)
