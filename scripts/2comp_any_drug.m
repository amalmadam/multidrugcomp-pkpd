%NOTE: this script is not interactive/will generate static graphs according to written parameters 
%between drugs/within patients comparison
function sim2comp()
    ke = log(2)/1; k12 = 0.1; k21 = 0.05;molmass = 456.32;  dose_mgMl = 0.5;
    concentration_nm_per_L = ((dose_mgMl * 1e3)/molmass)*1e6;
    patients = struct('Male_73kg', 5.475, 'Female_16kg', 1.040, 'Infant_4kg', 0.00366);
    step_times = [0, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60];  
    infusion_rates = [0.5, 1, 2, 4, 8, 16, 32, 16, 8, 4, 2, 0.5];     
    timeinterval = linspace(0, 60, 601); 
    patient_fields = fieldnames(patients);
    for i = 1:numel(patient_fields)
        patient_name = patient_fields{i};
        V1 = patients.(patient_name);  
        Cc = 0; Cp = 0;  
        Cc_values = zeros(1, length(timeinterval));
        Cp_values = zeros(1, length(timeinterval));
        matrixA = [- (ke + k12), k21; k12, -k21];
        for j = 1:length(timeinterval)
            t = timeinterval(j);
            infusion_rate_mL_hr = interp1(step_times, infusion_rates, t, 'previous'); %interpolation for concentration at time t
            infusion_rate_nm_per_min = (infusion_rate_mL_hr / 60) * concentration_nm_per_L;  
            if j == 1
                dt = timeinterval(j);  
            else
                dt = timeinterval(j) - timeinterval(j - 1);  
            end
            expA_dt = expm(matrixA * dt);  
            B = [infusion_rate_nm_per_min / V1; 0]; 
            concentrations = expA_dt * [Cc; Cp] + B * dt; Cc = concentrations(1); Cp = concentrations(2);
            Cc_values(j) = Cc; Cp_values(j) = Cp;
        end
        figure;
        plot(timeinterval,Cc_values,'b-','DisplayName','Central compartment');
        hold on;
        plot(timeinterval,Cp_values,'g-','DisplayName','Peripheral Compartment');
        xlabel('Time (min)');
        ylabel('Concentration (nmol/L)');
        title(['Clevidipine pharmacokinetics for ', strrep(patient_name, '_', ' ')]);
        legend;
        yyaxis right;
        plot(step_times, infusion_rates, 'r--','DisplayName','Flowrate (mL/hr)');
        ylabel('Set flowate (mL/hr)');
        grid on;
        hold off;
    end
end
