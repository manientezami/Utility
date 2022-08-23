function [ ] = csvAsColumnVectors2( fileName, headerRows, columnNames, formatSpec, ws )
%CSVASCOLUMNVECTORS Loads a CSV with a single header line as column vectors
    %fileName = full file path
    %headerRows = number of header rows (NOT including the column names)
    %columnNames = array of column names - if empty, will use first row after headerRows
    %formatSpec = column types (e.g. '%d %s %f %f %f %f' NOTE spaces not commas)
    %ws = workspace (default 'base')
    
    % Version 2 allows string columns
    
    if nargin < 5   
        ws = 'base';
    end
    
    F = strsplit(formatSpec,' ');
    nc = length(F);
    
    if (~isempty(columnNames)) && (nc ~= length(columnNames))
       error('Number of column names does not match formatSpec');
    end
    
    if headerRows > 0
        Headers = cell(headerRows, 1);
    end
    
    fid = fopen(fileName); %Open for text read
    for i=1:headerRows
        Headers{i} = fgets(fid);
    end
    
    if isempty(columnNames)
        fspec = '%s ';
        for i=2:nc
            fspec = [fspec '%s '];
        end;
        fspec = fspec(1:(end-1));
        firstline = textscan(fid,fspec,1,'HeaderLines',0,'Delimiter',','); %Get the first line, after headers
        
        vn = genvarname(firstline);
    else
        vn = genvarname(columnNames);
    end;
    
    Vectors = cell(0,nc);
    while ~feof(fid)
        cols = textscan(fid,formatSpec,1,'HeaderLines',0,'Delimiter',',');
        Vectors = [Vectors; cols];
    end;
    fclose(fid);
        
    for i=1:nc
        assignin(ws, vn{i}, [Vectors{:,i}]');
    end;
    
    if(headerRows > 0)
        assignin(ws, 'Headers', Headers);
    end
    
    %disp('Done.');

end

