function [] = saveSimResults()   
    filename = 'simResults.mat';

    yout = evalin('base','yout');
    tout = evalin('base','tout');
    quadModel = evalin('base','quadModel');
    IC = evalin('base','IC');
    pathExist = evalin('base', "exist('path','var')");
    
    if pathExist
        path = evalin('base','path');
        save(filename, 'yout', 'tout', 'quadModel', 'IC', 'path');
    else
        save(filename, 'yout', 'tout', 'quadModel', 'IC');
    end