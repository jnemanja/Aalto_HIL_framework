f = figure;
set(f, 'Position', [0 0 1000 400]);
dm = plot(magnetic_dipole.time, magnetic_dipole.signals.values);
xlabel('time, $s$', 'Interpreter', 'LaTex');
ylabel('dipole moment, $\frac{Nm}{T}$', 'Interpreter', 'LaTex');
export_fig './tests/' -eps -png -m4 dm
close

f = figure;
set(f, 'Position', [0 0 1000 400]);
qe = plot(quaternion_error.time, quaternion_error.signals.values);
xlabel('time, $s$', 'Interpreter', 'LaTex');
ylabel('quaternion error $q_{ob}$', 'Interpreter', 'LaTex');
export_fig './tests/' -eps -png -m4 qe
close

f = figure;
set(f, 'Position', [0 0 1000 400]);
ae = plot(angle_error.time, angle_error.signals.values);
xlabel('time, $s$', 'Interpreter', 'LaTex');
ylabel('angle error $\alpha$, ${deg}$', 'Interpreter', 'LaTex');
export_fig './tests/' -eps -png -m4 ae
close

f = figure;
set(f, 'Position', [0 0 1000 400]);
ae1 = plot(angle_error1.time, angle_error1.signals.values);
xlabel('time, $s$', 'Interpreter', 'LaTex');
ylabel('angle error $\alpha$, ${deg}$', 'Interpreter', 'LaTex');
export_fig './tests/' -eps -png -m4 ae2
close