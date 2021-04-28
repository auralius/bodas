# bodas
The asymptotic bode diagram in MATLAB

## Formalization

In order to use bodas, we must first decompose our transfer function as a product of:

G1(s)=1/(s+a)  
G2(s)=(s+b)  
G3(s)=1/s  
G4(s)=s  
G5(s)=K 

Please note that bodas CAN NOT handle complex poles/zeros. 

## A screenshot of how the result looks like

![Screenshot](sshot.png)


## How to use

A system with gain = 0.1, and zeros = [10 100], poles = [1]
```
% G(s) =  (s+10) * (s+100) * 0.1
%     
bodas([10 100], [1], 0.1)
```

A system with gain = -10, zeros = [0], poles = [1 1 10]. The bode is drawn in frequency 0.002 to 1000 rad/s
```
%              1       1       1
% G(s) = s * ----- * ----- * ------ * (-10)
%            (s+1)   (s+1)   (s+10)
%        
bodas([0],[1 1 10], -10, [-2 3]) 
```
Here, at the last argument, -2 corresponds to 10^-2 and 3 corresponds to 10^3.
