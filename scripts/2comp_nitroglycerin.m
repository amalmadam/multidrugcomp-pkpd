%NOTE:this script generates one static graph over set time interval with traces overlaid
%between patients/within drug comparison
function simulates_two_compartment_nitroglycerin()
    prompts = {'What is the drug name?', 'k12 (central --> peripheral rate):', 'k21 (peripheral --> central rate):', ...
               'Half-life (minutes):', 'Molar mass (g/mol):', 'Dose (mg/mL):', 'Number of patients:'};
    dlg_title = 'Simulation parameters';
    default = {'Nitroglycerin', '0.1', '0.05', '120', '227.09', '0.4', '3'}; 
    inputs = inputdlg(prompts, dlg_title, [1 50], default);
                    drug = inputs{1};
                    k12 = str2double(inputs{2});
                    k21 = str2double(inputs{3});
                    half_life = str2double(inputs{4});
                    mol_mass = str2double(inputs{5});
                    dose_mg_per_mL = str2double(inputs{6});
                    num_patients = str2double(inputs{7});
                        ke = log(2) / half_life;
                        conc_nmol_per_L = (dose_mg_per_mL * 1e3 / mol_mass) * 1e6;
    patient_weights = zeros(1, num_patients);
    for i = 1:num_patients
        prompt_weight = ['Enter weight in kilograms...'];
        weight_input = inputdlg(prompt_weight, 'What is the patient weight?', [1 50], {'70'});
        patient_weights(i) = str2double(weight_input{1});
    end
    V1 = patient_weights * 75 / 1000; %assuming avg. 70 ml blood per kg across all pts
    time_steps = [0, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60]; 
    infusion_rates = [0.5, 1, 2, 4, 8, 16, 32, 16, 8, 4, 2, 0.5];  
            sim_time = linspace(0, 60, 601);  %clinically relevant one hour window
    infusion_rates_interp = interp1(time_steps, infusion_rates, sim_time, 'previous');

    figure;
    hold on;
    traces = lines(num_patients);    
    for i = 1:num_patients
        Cc = 0; Cp = 0; 
        Cc_vals = zeros(1, length(sim_time)); Cp_vals = zeros(1, length(sim_time));
        matrixA = [- (ke + k12), k21; k12, -k21];
        for j = 1:length(sim_time)
                             t = sim_time(j);
                                    infusion_rate = infusion_rates_interp(j);
                                    infusion_nmol_per_min = (infusion_rate / 60) * conc_nmol_per_L;
                                 if j == 1
                                     dt = sim_time(j);
                                 else
                                     dt = sim_time(j) - sim_time(j - 1);
                                 end
            expA_dt = expm(matrixA * dt);
            B = [infusion_nmol_per_min / V1(i); 0];
            concentrations = expA_dt * [Cc; Cp] + B * dt;
            Cc = concentrations(1); Cp = concentrations(2);
            Cc_vals(j) = Cc; Cp_vals(j) = Cp;
        end
        plot(sim_time, Cc_vals, 'Color', traces(i, :), 'LineWidth', 2, 'DisplayName', ...
            ['Patient ', num2str(i), ' (', num2str(patient_weights(i)), ' kg) - Central']);
        plot(sim_time, Cp_vals, '--', 'Color', traces(i, :), 'LineWidth', 1.5, 'DisplayName', ...
            ['Patient ', num2str(i), ' (', num2str(patient_weights(i)), ' kg) - Peripheral']);
    end
    xlabel('Time (min)');
    ylabel('Concentration (nmol/L)');
    yyaxis right;
    plot(sim_time, infusion_rates_interp, 'r--', 'LineWidth', 2, 'DisplayName', 'Set flowrate (mL/hr)');
    ylabel('Flow rate (mL/hr)');
    title([drug, ' pharmacokinetics for patients of varying weight']);    
    legend('show');    
    grid on;
    hold off;
end
