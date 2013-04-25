function makefile( filename )   % Takes the initial position, and velocities of the particles, and the name of the file to save the data to.

    fid = fopen(filename, 'w'); % Create and open a file for writing.
    
    N = 36;              % Number of particles
    
    % The lowest energy interatomic distance is (2^(1/6))
    
    
    x(N,1) = 0;
    y(N,1) = 0;
    v(N,1) = 0;
    u(N,1) = 0;
    m(N,1) = 0;
    
    d=2^(1/6);
    
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
        
        v(i,1) = 0;
        u(i,1) = 0;
        m(i,1) = 1;
        
    end
    
    for i=1:N                   % Print off all x, y, u, v values.
        fprintf(fid, '%7.3f %7.3f %7.3f %7.3f %7.3f \n', x(i), y(i), v(i), u(i), m(i));
    end
    
    fclose(fid);

end