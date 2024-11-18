# bodas
The asymptotic bode diagram in MATLAB

## Formalization

In order to use bodas, we must first decompose our transfer function as a product of:

Real pole                  : ```G(s) = 1/(s+a)```  
Real zero                  : ```G(s) = (s+b)```  
Pole at the origin         : ```G(s) = 1/s```  
Zero at the origin         : ```G(s) = s```  
Constant gain      	       : ```G(s) = K```  
Complex conjugate zeros    : ```G(s) = (s+a+bi)*(s+a-bi)```   
Complex conjugate poles    : ```G(s) = 1 / ( (s+a+bi)*(s+a-bi) )```  

## A screenshot of how the result looks like

### Ex.1 with real poles/zeros

![Screenshot](sshot1.png)

![Screenshot](sshot2.png)

### Ex2. with complex conjugate poles/zeros

![Screenshot](sshot3.png)

![Screenshot](sshot4.png)

## How to use

Please read: bodas-manual.pdf
