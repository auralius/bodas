# bodas
The asymptotic bode diagram in MATLAB

## Formalization

An open loop transfer function is decomposed into the following forms:

G1(s)=a/(s+a)  
G2(s)=(s+b)/s  
G3(s)=1/s  
G4(s)=s  
G5(s)=K 

It is a must to formalize the system first. Also, bodas CAN NOT handle complex poles/zeros. 

## A screenshot of how the result looks like

![Screenshot](sshot.png)


## How to use

A system with gain = 0.1, and zeros = [10 100], poles = [1]
```
%         (s+10)   (s+100)     1
% G(s) =  ------ * ------- * ----- * 0.1
%           10       100     (s+1)
%     
bodas([10 100], [1], 0.1)
```

A system with gain = -10, zeros = [0], poles = [1 1 10]. The bode is drawn in frequency 0.002 to 1000 rad/s
```
%              1       1       10
% G(s) = s * ----- * ----- * ------ * (-10)
%            (s+1)   (s+1)   (s+10)
%        
bodas([0],[1 1 10], -10, [-2 3]) 
```
Here, in the last argument, -2 corresponds to 10^-2 and 3 corresponds to 10^3.
