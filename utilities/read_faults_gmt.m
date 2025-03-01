% 读取 .gmt 文件中的断层数据
function faults = read_faults_gmt(file_path)
    faults = {};
    fid = fopen(file_path, 'r');
    current_fault = [];
    while ~feof(fid)
        line = fgetl(fid);
        if startsWith(line, '>')
            if ~isempty(current_fault)
                faults{end+1} = current_fault;
                current_fault = [];
            end
        elseif ~isempty(line) && (isstrprop(line(1), 'digit') || line(1) == '-')
            data = sscanf(line, '%f %f');
            lon = data(1);
            lat = data(2);
            current_fault = [current_fault; lon, lat];
        end
    end
    if ~isempty(current_fault)
        faults{end+1} = current_fault;  % Append the last fault
    end
    fclose(fid);
end
