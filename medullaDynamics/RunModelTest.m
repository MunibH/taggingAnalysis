% Define filepaths and parameters
feats_filepath = 'Z:\OE\Data - All Neural and Kinematic Activity\Kinematics Data\YH35_240828_Kinematics_All.csv';
hidden_dim = 128;
num_layers = 3;
Animal = 0;
trial_len = 1480;
model_params_filepath = 'trained_model_parameters_20241112.pth';

% Construct the command to call the batch file with arguments
command = sprintf('RunModel.bat "%s" %d %d %d %d "%s"', ...
    feats_filepath, hidden_dim, num_layers, Animal, trial_len, model_params_filepath);

% Execute the command
[status, result] = system(command);

% Check for errors
if status == 0
    disp('Python script executed successfully');
    disp(result);

    % Load the predicted PCs
    Predicted_PCs = csvread('Predicted_PCs.csv');  % Read the array into MATLAB

    % Reshape if needed (based on trial length and number of PCs)
    Predicted_PCs = reshape(Predicted_PCs, trial_len, [], size(Predicted_PCs,2)); % [time, trials, PCs]
else
    disp('Error in executing Python script');
    disp(result);
end

%%
f = figure;
ax = gca;
hold on;
for i = 1:8
    cla(ax)
    hold on
    imagesc(Predicted_PCs(:,:,i)')
    ax.YDir = 'reverse';
    pause
end

















