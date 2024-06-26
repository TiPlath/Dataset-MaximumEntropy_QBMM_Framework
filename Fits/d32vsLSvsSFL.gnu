# guessed parameters to be used as initial values
a= 1 # zero is not allowed!
b= 1
c= 1
#
ffit2(t,h)=(a*t**2) * (h/(1550000))**(0.33333333333)
fit ffit2(x, y) 'd32vsLSvsSFL.txt' using 1:2:3 via a

# Plot the data and the fitted function
splot 'd32vsLSvsSFL.txt' using 1:2:3 with points, ffit2(x,y)
pause -1
