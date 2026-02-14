%NOTE: real time plot/1 drug/1 patient/1 compartment
%patient weight specified through V1 (blood volume) parameter
%default window parameters == nitroglycerin
function two_compartment_model_drug_specifics
help prompt
    prompt = {'Drug name:', 'Half-life (min):', 'Molar Mass (g/mol):', ...
              'Concentration in diluent (mg/mL):', 'V1 (L):', 'k12:', 'k21:', ...
              'Initial dose (mg):', 'Dose rate (mg/h):', 'Time step (hours):'};
    dlgtitle = 'Drug and its parameters';
    definput = {'DrugName', '2', '227.09', '0.4', '5', '0.1', '0.05', '100', '10', '0.1'};
    answer = inputdlg(prompt, dlgtitle, [1 50], definput);
    [drug_name, half_life, molar_mass, drug_concentration, V1, k12, k21, initial_dose, dose_rate, dt] = ...
        convert_inputs(answer);
    ke = log(2) / half_life; 
    initial_concentration = (initial_dose * 1e3) / (molar_mass * V1);  
    time_points = 0; central_array = initial_concentration;peripheral_array = 0;
    [hFig, hCentralLine, hPeripheralLine] = setup_figure(drug_name, time_points, central_array, peripheral_array);
    [hSliderK12, hSliderK21, hSliderKe, hSliderDoseRate] = setup_sliders(k12, k21, ke, dose_rate);
    isRunning = true;
    tic;
    while isRunning
        if toc >= dt
            tic;
            [k12, k21, ke, dose_rate] = get_slider_values(hSliderK12, hSliderK21, hSliderKe, hSliderDoseRate);
            [central, peripheral] = update_concentrations(central_array(end), peripheral_array(end), ...
                k12, k21, ke, dose_rate, molar_mass, V1, dt);            
            time_points(end+1) = time_points(end) + dt;
            central_array(end+1) = central;
            peripheral_array(end+1) = peripheral;

            set(hCentralLine, 'XData', time_points, 'YData', central_array);
            set(hPeripheralLine, 'XData', time_points, 'YData', peripheral_array);
            drawnow;
        end
    end
    function [drug_name, half_life, molar_mass, drug_concentration, V1, k12, k21, initial_dose, dose_rate, dt] = convert_inputs(answer)
        drug_name = answer{1};
        half_life = str2double(answer{2});
        molar_mass = str2double(answer{3});
        drug_concentration = str2double(answer{4});
        V1 = str2double(answer{5});
        k12 = str2double(answer{6});
        k21 = str2double(answer{7});
        initial_dose = str2double(answer{8});
        dose_rate = str2double(answer{9});
        dt = str2double(answer{10});
    end
    function [hFig, hCentralLine, hPeripheralLine] = setup_figure(drug_name, time_points, central_array, peripheral_array)
        hFig = figure('Name', ['Two-Compartment Model for ', drug_name], 'NumberTitle', 'off');
        hCentralLine = plot(time_points, central_array, 'r', 'DisplayName', 'Central Compartment');
        hold on;
        hPeripheralLine = plot(time_points, peripheral_array, 'b', 'DisplayName', 'Peripheral Compartment');
        xlabel('Time (hours)'); ylabel('Concentration (nmol/L)');
        title(['Two-Compartment Model for ', drug_name]);
        legend; grid on;
        hold off;
    end
    function [hSliderK12, hSliderK21, hSliderKe, hSliderDoseRate] = setup_sliders(k12, k21, ke, dose_rate)
        hSliderK12 = create_slider('k12', k12, [100, 130]);
        hSliderK21 = create_slider('k21', k21, [100, 100]);
        hSliderKe = create_slider('ke', ke, [100, 70]);
        hSliderDoseRate = create_slider('Dose Rate (mg/h)', dose_rate, [100, 40], 100);
    end
    function hSlider = create_slider(label, value, position, maxVal)
        if nargin < 4, maxVal = 1; end
        uicontrol('Style', 'text', 'Position', [20, position(2), 80, 20], 'String', label);
        hSlider = uicontrol('Style', 'slider', 'Min', 0, 'Max', maxVal, 'Value', value, ...
                            'Position', [position(1), position(2), 200, 20]);
        uicontrol('Style', 'text', 'Position', [310, position(2), 50, 20], 'String', num2str(value));
    end
    function [k12, k21, ke, dose_rate] = get_slider_values(hSliderK12, hSliderK21, hSliderKe, hSliderDoseRate)
        k12 = get(hSliderK12, 'Value');
        k21 = get(hSliderK21, 'Value');
        ke = get(hSliderKe, 'Value');
        dose_rate = get(hSliderDoseRate, 'Value');
    end
    function [central, peripheral] = update_concentrations(central, peripheral, k12, k21, ke, dose_rate, molar_mass, V1, dt)
        dose_input_nmoles = (dose_rate * dt * 1e3) / molar_mass;
        dC_central = (dose_input_nmoles / V1) - (ke + k12) * central + k21 * peripheral;
        dC_peripheral = k12 * central - k21 * peripheral;
        central = central + dC_central * dt;
        peripheral = peripheral + dC_peripheral * dt;
    end
end
