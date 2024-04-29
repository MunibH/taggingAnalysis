function time = getTimeVector(tmin,tmax,dt)

    edges = tmin:dt:tmax;
    time = edges + dt/2;
    time = time(1:end-1);
end