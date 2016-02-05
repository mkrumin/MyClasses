function data = loadTMazeData(obj)

filenames = dat.expFilePath(obj.expRef, 'TMaze');
try
    data = load(filenames{1});
catch
    try
        data = load(filenames{2});
    catch
        fprintf('No TMaze data found for ExpRef %s\n', ExpRef);
        return;
    end
end

