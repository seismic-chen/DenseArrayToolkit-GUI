function ver_out = FMI_ver;

% Return the current version number of FMI

ver_number = 0.91;

if nargout == 1
    ver_out = ver_number;
else
    fprintf('FMI version: %.2f\n\n', ver_number);
end;
