function da_dt = muscleActODEfun(t,a)


da_dt = (1.7*a*sin(pi/2*t)-a)/(.1*(.5+1.5*a)); %**Fill in your differential equation from part a. Don't forget
            % to delete the brackets first, or the code won't work**.