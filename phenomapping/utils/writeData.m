function writeData(filename,data,typeColumns,heading,typeHeading)

if (nargin < 4)
    heading = [];
end

if exist(filename,'file') == 2
    warning('the filename provided already exists. It will be overwritten')
end

fid = fopen(filename,'wt');
if fid>0
    if ~isempty(heading)
        fprintf(fid,strcat(typeHeading,'\n'),heading{1,:});
    end
    for k = 1:size(data,1)
        fprintf(fid,strcat(typeColumns,'\n'),data{k,:});
    end
    fclose(fid);
end