function [ ] = csvAsColumnVectors( fileName, headerRows, ws )
%CSVASCOLUMNVECTORS Loads a CSV with a single header line as column vectors
    %fileName = full file path, headerRows = number of header rows (NOT
    % including the column names (default 0), ws = workspace (default
    % 'base')
    
    if nargin < 3   
        ws = 'base';
        if nargin < 2
            headerRows = 0;
        end
    end
    
    disp('Loading csv data...');
    data = csvread(fileName,headerRows); %Load in the data, ignoring the header rows and column headers
    
    disp('Finished loading data. Now loading vector names...');
    
    if headerRows > 0
        Headers = cell(headerRows, 1);
    end
    
    fid = fopen(fileName); %Open for text read
    %firstline = textscan(fid,'%s',1,'HeaderLines',headerRows); %Get the first line
    for i=1:headerRows
        Headers{i} = fgets(fid);
    end
    firstline = textscan(fid,'%s',1,'HeaderLines',0); %Get the first line, after headers
    fclose(fid);
    
    paramStr = firstline{1}{1};
    p = 1;
    param = [];
    for c=1:length(paramStr)
        cc = paramStr(c);
        if cc == ','
            if ~isempty(param)
                paramIds{p} = param;
                param = [];
            else
                paramIds{p} = ['Column' num2str(p)];
            end
            p = p + 1;
        else
            param = [param cc];
        end
    end
    if ~isempty(param)
        paramIds{p} = param;
    else
        paramIds{p} = ['Column' num2str(p)];
    end
        
    
    disp('Now assigning vectors...');

    %Generate a column vector for each CSV column, using paramIds as 
    % variable names
    for i=1:size(data,2)
        if i <= length(paramIds)
            vn = genvarname(paramIds{i});
        else
            vn = genvarname(['Col' num2str(i)]);
        end

        assignin(ws, vn, data(:,i));
    end
    
    if(headerRows > 0)
        assignin(ws, 'Headers', Headers);
    end
    
    disp('Done.');

end

