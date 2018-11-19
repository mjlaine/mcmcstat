% publish options
popts = struct();
popts.format = 'html';
popts.outputDir = '../docs/ex/';
popts.iameFormat = 'png';

pub = @(file) publish(file,popts);
