function makefile( filename )   % Takes the initial position, and velocities of the particles, and the name of the file to save the data to.

    fid = fopen(filename, 'w'); % Create and open a file for writing.
    
    N = 36;              % Number of particles
	
    boltzmann = 1.381e-23;		%Gas Constant in J/kgK
	    
    x(N,1) = 0;
    y(N,1) = 0;
    v(N,1) = 0;
    u(N,1) = 0;
    m(N,1) = 0;
    
    d=2^(1/6);				%Lowest energy interatomic distance in nm
	Vrms = sqrt(3*boltzmann*300/(1.661e-27));	%calculate initial RMS velocity in nm/ns, assuming the mass of each particle is 1 AMU
    
    for i=1:N
        switch mod(i,4)
            case 1
                x(i,1) = floor(i/4)*d;
                y(i,1) = 0;
            case 2
                x(i,1) = (d/2)+floor(i/4)*d;
                y(i,1) = tand(60)*d/2;
            case 3
                x(i,1) = floor(i/4)*d;
                y(i,1) = 2*tand(60)*d/2;
            case 0
                x(i,1) = (d/2) + (floor(i/4)-1)*d;
                y(i,1) = 3*tand(60)*d/2;
        end
        
		%difference = rand(1);
		
		%%% Randomize the initial velocity distribution %%		
        %v(i,1) = Vrms*sqrt(1-difference);
        %u(i,1) = Vrms*sqrt(difference);
        v(i,1) = 0;
        u(i,1) = 0;
        m(i,1) = 1;
        
    end
    
    for i=1:N                   % Print off all x, y, u, v values.
        fprintf(fid, '%7.3f %7.3f %7.3f %7.3f %7.3f \n', x(i), y(i), v(i), u(i), m(i));
    end
    
    fclose(fid);

end