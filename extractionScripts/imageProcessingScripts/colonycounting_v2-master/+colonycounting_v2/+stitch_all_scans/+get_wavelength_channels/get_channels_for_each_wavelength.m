function channels = get_channels_for_each_wavelength(path_folder, name_scan, images)

    % get list of wavelengths:
    list_wavelengths = extractfield(images, 'wavelength');

    % ask user to assign a channel to each wavelength:
    prompt = [{path_folder}, list_wavelengths];
    default = cell(numel(list_wavelengths) + 1, 1);
    default(:,1) = {''};
    default{1} = name_scan;
    channels = inputdlg(prompt, 'Enter Channel for Each Wavelength', [1 200], default);
    
    % remove first answer:
    channels(1) = [];

end